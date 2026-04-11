using System.Diagnostics;
using Testbed;
using Testbed.Nudge;
using Testbed.Bepu;
using Testbed.Jolt;

var adapters = new List<Func<IPhysicsAdapter>>
{
	() => new NudgeAdapter(),
	() => new BepuAdapter(),
	() => new JoltAdapter(),
};

var scenarios = new (string name, Func<IPhysicsAdapter, int> setup)[]
{
	("StackBoxes_100",   a => Scenarios.StackBoxes(a, 100)),
	("StackBoxes_500",   a => Scenarios.StackBoxes(a, 500)),
	("SphereDrop_10x10", a => Scenarios.SphereDrop(a, 10)),
	("PyramidBoxes_20",  a => Scenarios.PyramidBoxes(a, 20)),
};

const float dt = 1.0f / 60.0f;
const int warmupSteps = 300;
const int measureSteps = 600;

Console.WriteLine("=".PadRight(110, '='));
Console.WriteLine($"  PHYSICS ENGINE BENCHMARK  (warmup={warmupSteps} steps, measure={measureSteps} steps)");
Console.WriteLine("=".PadRight(110, '='));
Console.WriteLine();

foreach (var (scenarioName, setup) in scenarios)
{
	Console.WriteLine($"--- {scenarioName} ---");
	Console.WriteLine($"{"Engine",-12} {"Bodies",6} {"Awake",6} {"Avg ms",8} {"Min ms",8} {"Max ms",8} {"P50 ms",8} {"P95 ms",8} {"Total s",8}");
	Console.WriteLine(new string('-', 90));

	foreach (var createAdapter in adapters)
	{
		try
		{
			using var adapter = createAdapter();
			adapter.CreateWorld(0, -9.81f, 0);
			int bodyCount = setup(adapter);

			for (int i = 0; i < warmupSteps; i++)
				adapter.Step(dt);

			var times = new double[measureSteps];
			bool isNudge = adapter is NudgeAdapter;
			var nudgeTimers = new NudgePerfTimers();
			int awakeSampleInterval = 120;
			var awakeLog = new List<(int step, int awake)>();

			var sw = Stopwatch.StartNew();
			for (int i = 0; i < measureSteps; i++)
			{
				adapter.Step(dt);
				times[i] = adapter.GetLastStepTimeMs();
				if (i % awakeSampleInterval == 0 || i == measureSteps - 1)
					awakeLog.Add((warmupSteps + i, adapter.GetActiveBodyCount()));

				if (isNudge)
				{
					var p = ((NudgeAdapter)adapter).GetPerfTimers();
					nudgeTimers.Broadphase += p.Broadphase;
					nudgeTimers.PreSolve += p.PreSolve;
					nudgeTimers.PgsSolve += p.PgsSolve;
					nudgeTimers.PositionCorrect += p.PositionCorrect;
					nudgeTimers.Integrate += p.Integrate;
					nudgeTimers.Islands += p.Islands;
					nudgeTimers.Total += p.Total;
					nudgeTimers.PgsPreSolve += p.PgsPreSolve;
					nudgeTimers.PgsWarmStart += p.PgsWarmStart;
					nudgeTimers.PgsGraphColor += p.PgsGraphColor;
					nudgeTimers.PgsIterations += p.PgsIterations;
					nudgeTimers.PgsJointLimits += p.PgsJointLimits;
					nudgeTimers.PgsLdl += p.PgsLdl;
					nudgeTimers.PgsRelax += p.PgsRelax;
					nudgeTimers.PgsPosContacts += p.PgsPosContacts;
					nudgeTimers.PgsPosJoints += p.PgsPosJoints;
					nudgeTimers.PgsPostSolve += p.PgsPostSolve;
				}
			}
			sw.Stop();

			// Debug: dump sleep diagnostics for Nudge
			if (isNudge) ((NudgeAdapter)adapter).DebugSleep();

			Array.Sort(times);
			double avg = times.Average();
			double min = times[0];
			double max = times[^1];
			double p50 = times[measureSteps / 2];
			double p95 = times[(int)(measureSteps * 0.95)];
			int awake = adapter.GetActiveBodyCount();

			Console.WriteLine($"{adapter.Name,-12} {bodyCount,6} {awake,6} {avg,8:F3} {min,8:F3} {max,8:F3} {p50,8:F3} {p95,8:F3} {sw.Elapsed.TotalSeconds,8:F3}");
			Console.Write($"  awake: ");
			foreach (var (step, aw) in awakeLog) Console.Write($"t{step/60.0:F0}s={aw} ");
			Console.WriteLine();

			if (isNudge)
			{
				double n = measureSteps;
				Console.WriteLine($"  phases: bp={nudgeTimers.Broadphase/n:F3} pre={nudgeTimers.PreSolve/n:F3} pgs={nudgeTimers.PgsSolve/n:F3} pos={nudgeTimers.PositionCorrect/n:F3} int={nudgeTimers.Integrate/n:F3} isl={nudgeTimers.Islands/n:F3}");
				Console.WriteLine($"  pgs:    pre={nudgeTimers.PgsPreSolve/n:F3} iter={nudgeTimers.PgsIterations/n:F3} relax={nudgeTimers.PgsRelax/n:F3} post={nudgeTimers.PgsPostSolve/n:F3}");
			}
		}
		catch (Exception ex)
		{
			Console.WriteLine($"{"ERROR",-12} {ex.GetType().Name}: {ex.Message}");
		}
	}
	Console.WriteLine();
}
