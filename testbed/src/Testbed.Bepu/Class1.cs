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

class MaterialLookup
{
	public readonly Dictionary<BodyHandle, float> DynamicFriction = new();
	public readonly Dictionary<StaticHandle, float> StaticFriction = new();

	public float GetFriction(CollidableReference r)
	{
		if (r.Mobility == CollidableMobility.Static)
		{
			if (StaticFriction.TryGetValue(r.StaticHandle, out float f)) return f;
		}
		else
		{
			if (DynamicFriction.TryGetValue(r.BodyHandle, out float f)) return f;
		}
		return 0.5f;
	}
}

struct NarrowPhaseCallbacks : INarrowPhaseCallbacks
{
	public SpringSettings ContactSpringiness;
	public MaterialLookup Materials;

	public void Initialize(Simulation simulation)
	{
		if (ContactSpringiness.AngularFrequency == 0)
			ContactSpringiness = new SpringSettings(30, 1);
	}

	public bool AllowContactGeneration(int workerIndex, CollidableReference a, CollidableReference b, ref float speculativeMargin) => true;
	public bool AllowContactGeneration(int workerIndex, CollidablePair pair, int childIndexA, int childIndexB) => true;

	public bool ConfigureContactManifold<TManifold>(int workerIndex, CollidablePair pair, ref TManifold manifold, out PairMaterialProperties pairMaterial) where TManifold : unmanaged, IContactManifold<TManifold>
	{
		float f1 = Materials?.GetFriction(pair.A) ?? 0.5f;
		float f2 = Materials?.GetFriction(pair.B) ?? 0.5f;
		pairMaterial.FrictionCoefficient = MathF.Sqrt(f1 * f2);
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
	readonly MaterialLookup _materials = new();
	double _lastStepMs;

	static Quaternion NormalizeRot(BodyDesc d) =>
		(d.RotX == 0 && d.RotY == 0 && d.RotZ == 0 && d.RotW == 0) ? Quaternion.Identity : new Quaternion(d.RotX, d.RotY, d.RotZ, d.RotW);

	public void CreateWorld(float gravityX, float gravityY, float gravityZ)
	{
		_pool = new BufferPool();
		_sim = Simulation.Create(
			_pool,
			new NarrowPhaseCallbacks { Materials = _materials },
			new PoseIntegratorCallbacks { Gravity = new Vector3(gravityX, gravityY, gravityZ) },
			new SolveDescription(8, 1));
	}

	public int AddBody(BodyDesc desc)
	{
		int index = _isStatic.Count;
		var rot = NormalizeRot(desc);
		var pose = new RigidPose(new Vector3(desc.PosX, desc.PosY, desc.PosZ), rot);
		float friction = desc.Friction > 0 ? desc.Friction : 0.5f;

		if (desc.Mass <= 0)
		{
			// Static body
			TypedIndex shapeIndex = AddShape(desc);
			var handle = _sim!.Statics.Add(new StaticDescription(pose, shapeIndex));
			_staticHandles.Add(handle);
			_dynamicHandles.Add(default);
			_isStatic.Add(true);
			_materials.StaticFriction[handle] = friction;
		}
		else
		{
			// Dynamic body
			switch (desc.Shape)
			{
				case ShapeType.Box:
				{
					var shape = new Box(desc.HalfExtentX * 2, desc.HalfExtentY * 2, desc.HalfExtentZ * 2);
					var bodyDesc = BodyDescription.CreateConvexDynamic(pose, desc.Mass, _sim!.Shapes, shape);
					var handle = _sim!.Bodies.Add(bodyDesc);
					_dynamicHandles.Add(handle);
					_materials.DynamicFriction[handle] = friction;
					break;
				}
				case ShapeType.Sphere:
				{
					var shape = new BepuPhysics.Collidables.Sphere(desc.Radius);
					var bodyDesc = BodyDescription.CreateConvexDynamic(pose, desc.Mass, _sim!.Shapes, shape);
					var handle = _sim!.Bodies.Add(bodyDesc);
					_dynamicHandles.Add(handle);
					_materials.DynamicFriction[handle] = friction;
					break;
				}
				case ShapeType.Capsule:
				{
					var shape = new Capsule(desc.Radius, desc.HalfHeight * 2);
					var bodyDesc = BodyDescription.CreateConvexDynamic(pose, desc.Mass, _sim!.Shapes, shape);
					var handle = _sim!.Bodies.Add(bodyDesc);
					_dynamicHandles.Add(handle);
					_materials.DynamicFriction[handle] = friction;
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

	BodyHandle EnsureDynamic(int bodyIndex)
	{
		if (!_isStatic[bodyIndex])
			return _dynamicHandles[bodyIndex];

		// Convert static to kinematic so it can participate in constraints
		var staticRef = _sim!.Statics[_staticHandles[bodyIndex]];
		var pose = staticRef.Pose;
		var shape = staticRef.Shape;
		var oldHandle = _staticHandles[bodyIndex];
		_sim.Statics.Remove(oldHandle);

		var desc = BodyDescription.CreateKinematic(pose, new CollidableDescription(shape), new BodyActivityDescription(-1));
		var handle = _sim.Bodies.Add(desc);
		_dynamicHandles[bodyIndex] = handle;
		_isStatic[bodyIndex] = false;

		// Migrate friction data
		if (_materials.StaticFriction.Remove(oldHandle, out float f))
			_materials.DynamicFriction[handle] = f;

		return handle;
	}

	struct DragState { public BodyHandle Anchor; public ConstraintHandle Constraint; }
	readonly Dictionary<int, DragState> _drags = new();
	int _nextDragId = 1;

	public int BeginDrag(int bodyIndex, float hitX, float hitY, float hitZ)
	{
		var anchorDesc = BodyDescription.CreateKinematic(
			new RigidPose(new Vector3(hitX, hitY, hitZ)),
			new CollidableDescription(_sim!.Shapes.Add(new BepuPhysics.Collidables.Sphere(0.01f))),
			new BodyActivityDescription(-1));
		var anchor = _sim!.Bodies.Add(anchorDesc);

		var dynHandle = EnsureDynamic(bodyIndex);

		// Compute body-local offset
		var bodyPose = _sim.Bodies[dynHandle].Pose;
		var worldOffset = new Vector3(hitX, hitY, hitZ) - bodyPose.Position;
		var invRot = Quaternion.Conjugate(bodyPose.Orientation);
		var localOffset = Vector3.Transform(worldOffset, invRot);

		var constraint = _sim.Solver.Add(anchor, dynHandle, new BallSocket
		{
			LocalOffsetA = Vector3.Zero,
			LocalOffsetB = localOffset,
			SpringSettings = new SpringSettings(5, 0.7f),
		});

		int id = _nextDragId++;
		_drags[id] = new DragState { Anchor = anchor, Constraint = constraint };
		return id;
	}

	public void UpdateDrag(int dragHandle, float targetX, float targetY, float targetZ)
	{
		if (!_drags.TryGetValue(dragHandle, out var state)) return;
		var bodyRef = _sim!.Bodies.GetBodyReference(state.Anchor);
		bodyRef.Pose.Position = new Vector3(targetX, targetY, targetZ);
		bodyRef.Awake = true;
	}

	public void EndDrag(int dragHandle)
	{
		if (!_drags.Remove(dragHandle, out var state)) return;
		_sim!.Solver.Remove(state.Constraint);
		_sim.Bodies.Remove(state.Anchor);
	}

	public void AddDistanceJoint(int bodyA, int bodyB, float localAx, float localAy, float localAz, float localBx, float localBy, float localBz, float restLength)
	{
		var hA = EnsureDynamic(bodyA);
		var hB = EnsureDynamic(bodyB);
		_sim!.Solver.Add(hA, hB, new DistanceServo
		{
			LocalOffsetA = new Vector3(localAx, localAy, localAz),
			LocalOffsetB = new Vector3(localBx, localBy, localBz),
			TargetDistance = restLength,
			SpringSettings = new SpringSettings(120, 1),
		});
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
