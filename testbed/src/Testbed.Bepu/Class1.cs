using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using BepuPhysics;
using BepuPhysics.Collidables;
using BepuPhysics.CollisionDetection;
using BepuPhysics.Constraints;
using BepuUtilities;
using BepuUtilities.Memory;
using Testbed;

namespace Testbed.Bepu;

struct NarrowPhaseCallbacks : INarrowPhaseCallbacks
{
	public SpringSettings ContactSpringiness;

	public void Initialize(Simulation simulation)
	{
		if (ContactSpringiness.AngularFrequency == 0)
			ContactSpringiness = new SpringSettings(30, 1);
	}

	public bool AllowContactGeneration(int workerIndex, CollidableReference a, CollidableReference b, ref float speculativeMargin) => true;
	public bool AllowContactGeneration(int workerIndex, CollidablePair pair, int childIndexA, int childIndexB) => true;

	public bool ConfigureContactManifold<TManifold>(int workerIndex, CollidablePair pair, ref TManifold manifold, out PairMaterialProperties pairMaterial) where TManifold : unmanaged, IContactManifold<TManifold>
	{
		pairMaterial.FrictionCoefficient = 0.5f;
		pairMaterial.MaximumRecoveryVelocity = 2f;
		pairMaterial.SpringSettings = ContactSpringiness;
		return true;
	}

	public bool ConfigureContactManifold(int workerIndex, CollidablePair pair, int childIndexA, int childIndexB, ref ConvexContactManifold manifold) => true;
	public void Dispose() { }
}

struct PoseIntegratorCallbacks : IPoseIntegratorCallbacks
{
	public Vector3 Gravity;
	Vector3Wide gravityWide;

	public readonly AngularIntegrationMode AngularIntegrationMode => AngularIntegrationMode.Nonconserving;
	public readonly bool AllowSubstepsForUnconstrainedBodies => false;
	public readonly bool IntegrateVelocityForKinematics => false;

	public void Initialize(Simulation simulation) { }

	public void PrepareForIntegration(float dt)
	{
		gravityWide = Vector3Wide.Broadcast(Gravity * dt);
	}

	public void IntegrateVelocity(Vector<int> bodyIndices, Vector3Wide position, QuaternionWide orientation, BodyInertiaWide localInertia, Vector<int> integrationMask, int workerIndex, Vector<float> dt, ref BodyVelocityWide velocity)
	{
		velocity.Linear += gravityWide;
	}
}

public class BepuAdapter : IPhysicsAdapter
{
	public string Name => "Bepu";

	BufferPool? _pool;
	Simulation? _sim;
	readonly List<BodyHandle> _dynamicHandles = new();
	readonly List<StaticHandle> _staticHandles = new();
	readonly List<bool> _isStatic = new();
	double _lastStepMs;

	public void CreateWorld(float gravityX, float gravityY, float gravityZ)
	{
		_pool = new BufferPool();
		_sim = Simulation.Create(
			_pool,
			new NarrowPhaseCallbacks(),
			new PoseIntegratorCallbacks { Gravity = new Vector3(gravityX, gravityY, gravityZ) },
			new SolveDescription(8, 1));
	}

	public int AddBody(BodyDesc desc)
	{
		int index = _isStatic.Count;

		if (desc.Mass <= 0)
		{
			// Static body
			TypedIndex shapeIndex = AddShape(desc);
			var staticDesc = new StaticDescription(
				new Vector3(desc.PosX, desc.PosY, desc.PosZ),
				shapeIndex);
			_staticHandles.Add(_sim!.Statics.Add(staticDesc));
			_dynamicHandles.Add(default);
			_isStatic.Add(true);
		}
		else
		{
			// Dynamic body
			switch (desc.Shape)
			{
				case ShapeType.Box:
				{
					var shape = new Box(desc.HalfExtentX * 2, desc.HalfExtentY * 2, desc.HalfExtentZ * 2);
					var bodyDesc = BodyDescription.CreateConvexDynamic(
						new RigidPose(new Vector3(desc.PosX, desc.PosY, desc.PosZ)),
						desc.Mass, _sim!.Shapes, shape);
					_dynamicHandles.Add(_sim!.Bodies.Add(bodyDesc));
					break;
				}
				case ShapeType.Sphere:
				{
					var shape = new BepuPhysics.Collidables.Sphere(desc.Radius);
					var bodyDesc = BodyDescription.CreateConvexDynamic(
						new RigidPose(new Vector3(desc.PosX, desc.PosY, desc.PosZ)),
						desc.Mass, _sim!.Shapes, shape);
					_dynamicHandles.Add(_sim!.Bodies.Add(bodyDesc));
					break;
				}
				case ShapeType.Capsule:
				{
					var shape = new Capsule(desc.Radius, desc.HalfHeight * 2);
					var bodyDesc = BodyDescription.CreateConvexDynamic(
						new RigidPose(new Vector3(desc.PosX, desc.PosY, desc.PosZ)),
						desc.Mass, _sim!.Shapes, shape);
					_dynamicHandles.Add(_sim!.Bodies.Add(bodyDesc));
					break;
				}
			}
			_staticHandles.Add(default);
			_isStatic.Add(false);
		}

		return index;
	}

	TypedIndex AddShape(BodyDesc desc)
	{
		return desc.Shape switch
		{
			ShapeType.Box => _sim!.Shapes.Add(new Box(desc.HalfExtentX * 2, desc.HalfExtentY * 2, desc.HalfExtentZ * 2)),
			ShapeType.Sphere => _sim!.Shapes.Add(new BepuPhysics.Collidables.Sphere(desc.Radius)),
			ShapeType.Capsule => _sim!.Shapes.Add(new Capsule(desc.Radius, desc.HalfHeight * 2)),
			_ => throw new ArgumentException("Unknown shape type"),
		};
	}

	public void Step(float dt)
	{
		var sw = Stopwatch.StartNew();
		_sim!.Timestep(dt);
		sw.Stop();
		_lastStepMs = sw.Elapsed.TotalMilliseconds;
	}

	public (float x, float y, float z) GetPosition(int bodyIndex)
	{
		if (_isStatic[bodyIndex])
		{
			var pose = _sim!.Statics[_staticHandles[bodyIndex]].Pose;
			return (pose.Position.X, pose.Position.Y, pose.Position.Z);
		}
		else
		{
			var pose = _sim!.Bodies[_dynamicHandles[bodyIndex]].Pose;
			return (pose.Position.X, pose.Position.Y, pose.Position.Z);
		}
	}

	public (float x, float y, float z, float w) GetRotation(int bodyIndex)
	{
		if (_isStatic[bodyIndex])
		{
			var q = _sim!.Statics[_staticHandles[bodyIndex]].Pose.Orientation;
			return (q.X, q.Y, q.Z, q.W);
		}
		else
		{
			var q = _sim!.Bodies[_dynamicHandles[bodyIndex]].Pose.Orientation;
			return (q.X, q.Y, q.Z, q.W);
		}
	}

	public void SetVelocity(int bodyIndex, float vx, float vy, float vz)
	{
		if (!_isStatic[bodyIndex])
			_sim!.Bodies[_dynamicHandles[bodyIndex]].Velocity.Linear = new Vector3(vx, vy, vz);
	}

	public double GetLastStepTimeMs() => _lastStepMs;

	public int GetActiveBodyCount()
	{
		int count = 0;
		for (int i = 0; i < _isStatic.Count; i++)
		{
			if (!_isStatic[i] && _sim!.Bodies[_dynamicHandles[i]].Awake) count++;
		}
		return count;
	}

	public bool IsBodyActive(int bodyIndex)
	{
		if (_isStatic[bodyIndex]) return false;
		return _sim!.Bodies[_dynamicHandles[bodyIndex]].Awake;
	}

	public void Dispose()
	{
		_sim?.Dispose();
		_pool?.Clear();
		_sim = null;
		_pool = null;
	}
}
