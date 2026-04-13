using System.Diagnostics;
using System.Runtime.InteropServices;
using Testbed;

namespace Testbed.Nudge;

internal static partial class Native
{
	const string Lib = "nudge";

	[LibraryImport(Lib, EntryPoint = "nudge_create_world")]
	public static partial ulong CreateWorld(float gx, float gy, float gz, int solverType, int subSteps, int velocityIters);

	[LibraryImport(Lib, EntryPoint = "nudge_destroy_world")]
	public static partial void DestroyWorld(ulong world);

	[LibraryImport(Lib, EntryPoint = "nudge_step")]
	public static partial void Step(ulong world, float dt);

	[LibraryImport(Lib, EntryPoint = "nudge_get_step_time")]
	public static partial double GetStepTime(ulong world);

	[LibraryImport(Lib, EntryPoint = "nudge_create_body")]
	public static partial ulong CreateBody(ulong world, float px, float py, float pz, float mass, float friction, float restitution);

	[LibraryImport(Lib, EntryPoint = "nudge_create_body_rotated")]
	public static partial ulong CreateBodyRotated(ulong world, float px, float py, float pz, float qx, float qy, float qz, float qw, float mass, float friction, float restitution);

	[LibraryImport(Lib, EntryPoint = "nudge_body_add_box")]
	public static partial void BodyAddBox(ulong world, ulong body, float hx, float hy, float hz);

	[LibraryImport(Lib, EntryPoint = "nudge_body_add_sphere")]
	public static partial void BodyAddSphere(ulong world, ulong body, float radius);

	[LibraryImport(Lib, EntryPoint = "nudge_body_add_capsule")]
	public static partial void BodyAddCapsule(ulong world, ulong body, float halfHeight, float radius);

	[LibraryImport(Lib, EntryPoint = "nudge_get_position")]
	public static unsafe partial void GetPosition(ulong world, ulong body, float* outXyz);

	[LibraryImport(Lib, EntryPoint = "nudge_get_rotation")]
	public static unsafe partial void GetRotation(ulong world, ulong body, float* outXyzw);

	[LibraryImport(Lib, EntryPoint = "nudge_body_is_asleep")]
	public static partial int BodyIsAsleep(ulong world, ulong body);

	[LibraryImport(Lib, EntryPoint = "nudge_get_perf")]
	public static unsafe partial void GetPerf(ulong world, double* out17);

	[LibraryImport(Lib, EntryPoint = "nudge_debug_sleep")]
	public static partial void DebugSleep(ulong world);

	[LibraryImport(Lib, EntryPoint = "nudge_body_set_velocity")]
	public static partial void SetVelocity(ulong world, ulong body, float vx, float vy, float vz);

	[LibraryImport(Lib, EntryPoint = "nudge_create_ball_socket")]
	public static partial ulong CreateBallSocket(ulong world, ulong bodyA, ulong bodyB, float ax, float ay, float az, float bx, float by, float bz, float freq, float dampingRatio);

	[LibraryImport(Lib, EntryPoint = "nudge_body_set_position")]
	public static partial void SetPosition(ulong world, ulong body, float px, float py, float pz);

	[LibraryImport(Lib, EntryPoint = "nudge_destroy_body")]
	public static partial void DestroyBody(ulong world, ulong body);

	[LibraryImport(Lib, EntryPoint = "nudge_destroy_joint")]
	public static partial void DestroyJoint(ulong world, ulong joint);

	[LibraryImport(Lib, EntryPoint = "nudge_create_distance_joint")]
	public static partial ulong CreateDistanceJoint(ulong world, ulong bodyA, ulong bodyB, float ax, float ay, float az, float bx, float by, float bz, float restLength);
}

public struct NudgePerfTimers
{
	// Top-level phases (ms)
	public double Broadphase, PreSolve, PgsSolve, PositionCorrect, Integrate, Islands, Total;
	// PGS sub-timers (ms)
	public double PgsPreSolve, PgsWarmStart, PgsGraphColor, PgsIterations;
	public double PgsJointLimits, PgsLdl, PgsRelax, PgsPosContacts, PgsPosJoints, PgsPostSolve;
}

public class NudgeAdapter : IPhysicsAdapter
{
	public string Name => "Nudge";
	public ulong WorldId => _world;

	ulong _world;
	readonly List<ulong> _bodies = new();
	double _lastStepMs;

	static (float x, float y, float z, float w) NormalizeRot(BodyDesc d) =>
		(d.RotX == 0 && d.RotY == 0 && d.RotZ == 0 && d.RotW == 0) ? (0, 0, 0, 1f) : (d.RotX, d.RotY, d.RotZ, d.RotW);

	public void CreateWorld(float gravityX, float gravityY, float gravityZ)
	{
		_world = Native.CreateWorld(gravityX, gravityY, gravityZ, 0, 0, 0);
	}

	public int AddBody(BodyDesc desc)
	{
		var (qx, qy, qz, qw) = NormalizeRot(desc);
		ulong body = Native.CreateBodyRotated(_world, desc.PosX, desc.PosY, desc.PosZ, qx, qy, qz, qw, desc.Mass, desc.Friction, desc.Restitution);

		switch (desc.Shape)
		{
			case ShapeType.Box:
				Native.BodyAddBox(_world, body, desc.HalfExtentX, desc.HalfExtentY, desc.HalfExtentZ);
				break;
			case ShapeType.Sphere:
				Native.BodyAddSphere(_world, body, desc.Radius);
				break;
			case ShapeType.Capsule:
				Native.BodyAddCapsule(_world, body, desc.HalfHeight, desc.Radius);
				break;
		}

		int index = _bodies.Count;
		_bodies.Add(body);
		return index;
	}

	public void Step(float dt)
	{
		Native.Step(_world, dt);
		_lastStepMs = Native.GetStepTime(_world) * 1000.0;
	}

	public unsafe (float x, float y, float z) GetPosition(int bodyIndex)
	{
		float* pos = stackalloc float[3];
		Native.GetPosition(_world, _bodies[bodyIndex], pos);
		return (pos[0], pos[1], pos[2]);
	}

	public unsafe (float x, float y, float z, float w) GetRotation(int bodyIndex)
	{
		float* rot = stackalloc float[4];
		Native.GetRotation(_world, _bodies[bodyIndex], rot);
		return (rot[0], rot[1], rot[2], rot[3]);
	}

	public double GetLastStepTimeMs() => _lastStepMs;

	public int GetActiveBodyCount()
	{
		int count = 0;
		foreach (var b in _bodies)
			if (Native.BodyIsAsleep(_world, b) == 0) count++;
		return count;
	}

	public bool IsBodyActive(int bodyIndex) => Native.BodyIsAsleep(_world, _bodies[bodyIndex]) == 0;

	public void DebugSleep() => Native.DebugSleep(_world);

	public void SetVelocity(int bodyIndex, float vx, float vy, float vz) => Native.SetVelocity(_world, _bodies[bodyIndex], vx, vy, vz);

	struct DragState { public ulong AnchorBody, Joint; }
	readonly Dictionary<int, DragState> _drags = new();
	int _nextDragId = 1;

	public int BeginDrag(int bodyIndex, float hitX, float hitY, float hitZ)
	{
		// Create static anchor body at hit point
		ulong anchor = Native.CreateBody(_world, hitX, hitY, hitZ, 0, 0.5f, 0);
		Native.BodyAddSphere(_world, anchor, 0.01f);

		// Compute body-local offset: inverse-rotate (hit - bodyPos)
		var (bx, by, bz) = GetPosition(bodyIndex);
		var (qx, qy, qz, qw) = GetRotation(bodyIndex);
		float dx = hitX - bx, dy = hitY - by, dz = hitZ - bz;
		// q^-1 * d: conjugate quaternion rotation
		float tx = 2 * (qy * dz - qz * dy), ty = 2 * (qz * dx - qx * dz), tz = 2 * (qx * dy - qy * dx);
		float lx = dx - qw * tx - (qy * tz - qz * ty);
		float ly = dy - qw * ty - (qz * tx - qx * tz);
		float lz = dz - qw * tz - (qx * ty - qy * tx);

		ulong joint = Native.CreateBallSocket(_world, anchor, _bodies[bodyIndex], 0, 0, 0, lx, ly, lz, 5.0f, 0.7f);

		int id = _nextDragId++;
		_drags[id] = new DragState { AnchorBody = anchor, Joint = joint };
		return id;
	}

	public void UpdateDrag(int dragHandle, float targetX, float targetY, float targetZ)
	{
		if (!_drags.TryGetValue(dragHandle, out var state)) return;
		Native.SetPosition(_world, state.AnchorBody, targetX, targetY, targetZ);
	}

	public void EndDrag(int dragHandle)
	{
		if (!_drags.Remove(dragHandle, out var state)) return;
		Native.DestroyJoint(_world, state.Joint);
		Native.DestroyBody(_world, state.AnchorBody);
	}

	public void AddDistanceJoint(int bodyA, int bodyB, float localAx, float localAy, float localAz, float localBx, float localBy, float localBz, float restLength)
	{
		Native.CreateDistanceJoint(_world, _bodies[bodyA], _bodies[bodyB], localAx, localAy, localAz, localBx, localBy, localBz, restLength);
	}

	public unsafe NudgePerfTimers GetPerfTimers()
	{
		double* t = stackalloc double[17];
		Native.GetPerf(_world, t);
		return new NudgePerfTimers
		{
			Broadphase = t[0]*1000, PreSolve = t[1]*1000, PgsSolve = t[2]*1000,
			PositionCorrect = t[3]*1000, Integrate = t[4]*1000, Islands = t[5]*1000, Total = t[6]*1000,
			PgsPreSolve = t[7]*1000, PgsWarmStart = t[8]*1000, PgsGraphColor = t[9]*1000,
			PgsIterations = t[10]*1000, PgsJointLimits = t[11]*1000, PgsLdl = t[12]*1000,
			PgsRelax = t[13]*1000, PgsPosContacts = t[14]*1000, PgsPosJoints = t[15]*1000,
			PgsPostSolve = t[16]*1000,
		};
	}

	public unsafe string GetPerfBreakdown()
	{
		var p = GetPerfTimers();
		return $"bp={p.Broadphase:F2} pre={p.PreSolve:F2} pgs={p.PgsSolve:F2} pos={p.PositionCorrect:F2} int={p.Integrate:F2} isl={p.Islands:F2}";
	}

	public void Dispose()
	{
		if (_world != 0)
		{
			foreach (var state in _drags.Values)
			{
				Native.DestroyJoint(_world, state.Joint);
				Native.DestroyBody(_world, state.AnchorBody);
			}
			_drags.Clear();
			Native.DestroyWorld(_world);
			_world = 0;
		}
	}
}
