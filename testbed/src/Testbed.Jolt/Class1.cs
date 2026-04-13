using Testbed;
using System.Runtime.InteropServices;

namespace Testbed.Jolt;

internal static partial class Native
{
	const string Lib = "jolt";

	[LibraryImport(Lib, EntryPoint = "jolt_create_world")]
	public static partial nint CreateWorld(float gx, float gy, float gz, int maxBodies);

	[LibraryImport(Lib, EntryPoint = "jolt_destroy_world")]
	public static partial void DestroyWorld(nint world);

	// shape_type: 0=box, 1=sphere, 2=capsule
	// s0/s1/s2: box={hx,hy,hz}, sphere={radius,_,_}, capsule={half_height,radius,_}
	[LibraryImport(Lib, EntryPoint = "jolt_create_body")]
	public static partial int CreateBody(nint world, int shapeType, float s0, float s1, float s2, float px, float py, float pz, float mass, float friction, float restitution);

	[LibraryImport(Lib, EntryPoint = "jolt_create_body_rotated")]
	public static partial int CreateBodyRotated(nint world, int shapeType, float s0, float s1, float s2, float px, float py, float pz, float qx, float qy, float qz, float qw, float mass, float friction, float restitution);

	[LibraryImport(Lib, EntryPoint = "jolt_step")]
	public static partial void Step(nint world, float dt);

	[LibraryImport(Lib, EntryPoint = "jolt_get_step_time")]
	public static partial double GetStepTime(nint world);

	[LibraryImport(Lib, EntryPoint = "jolt_get_position")]
	public static unsafe partial void GetPosition(nint world, int bodyIndex, float* outXyz);

	[LibraryImport(Lib, EntryPoint = "jolt_get_rotation")]
	public static unsafe partial void GetRotation(nint world, int bodyIndex, float* outXyzw);

	[LibraryImport(Lib, EntryPoint = "jolt_is_body_active")]
	public static partial int IsBodyActive(nint world, int bodyIndex);

	[LibraryImport(Lib, EntryPoint = "jolt_get_active_count")]
	public static partial int GetActiveCount(nint world);

	[LibraryImport(Lib, EntryPoint = "jolt_set_velocity")]
	public static partial void SetVelocity(nint world, int bodyIndex, float vx, float vy, float vz);

	[LibraryImport(Lib, EntryPoint = "jolt_optimize_broadphase")]
	public static partial void OptimizeBroadphase(nint world);

	[LibraryImport(Lib, EntryPoint = "jolt_create_distance_joint")]
	public static partial void CreateDistanceJoint(nint world, int bodyA, int bodyB, float ax, float ay, float az, float bx, float by, float bz, float restLength);

	[LibraryImport(Lib, EntryPoint = "jolt_create_kinematic_body")]
	public static partial int CreateKinematicBody(nint world, float px, float py, float pz);

	[LibraryImport(Lib, EntryPoint = "jolt_set_position")]
	public static partial void SetPosition(nint world, int bodyIndex, float px, float py, float pz);

	[LibraryImport(Lib, EntryPoint = "jolt_remove_body")]
	public static partial void RemoveBody(nint world, int bodyIndex);

	[LibraryImport(Lib, EntryPoint = "jolt_create_spring_constraint")]
	public static partial int CreateSpringConstraint(nint world, int bodyA, int bodyB, float ax, float ay, float az, float bx, float by, float bz, float freq, float damping);

	[LibraryImport(Lib, EntryPoint = "jolt_remove_constraint")]
	public static partial void RemoveConstraint(nint world, int constraintId);
}

public class JoltAdapter : IPhysicsAdapter
{
	public string Name => "Jolt";

	nint _world;
	int _bodyCount;
	bool _broadphaseOptimized;

	static (float x, float y, float z, float w) NormalizeRot(BodyDesc d) =>
		(d.RotX == 0 && d.RotY == 0 && d.RotZ == 0 && d.RotW == 0) ? (0, 0, 0, 1f) : (d.RotX, d.RotY, d.RotZ, d.RotW);

	public void CreateWorld(float gravityX, float gravityY, float gravityZ)
	{
		_world = Native.CreateWorld(gravityX, gravityY, gravityZ, 65536);
		_broadphaseOptimized = false;
	}

	public int AddBody(BodyDesc desc)
	{
		int shapeType;
		float s0, s1, s2;

		switch (desc.Shape)
		{
			case ShapeType.Box:
				shapeType = 0;
				s0 = desc.HalfExtentX; s1 = desc.HalfExtentY; s2 = desc.HalfExtentZ;
				break;
			case ShapeType.Sphere:
				shapeType = 1;
				s0 = desc.Radius; s1 = 0; s2 = 0;
				break;
			case ShapeType.Capsule:
				shapeType = 2;
				s0 = desc.HalfHeight; s1 = desc.Radius; s2 = 0;
				break;
			default:
				throw new ArgumentException("Unknown shape type");
		}

		var (qx, qy, qz, qw) = NormalizeRot(desc);
		int index = Native.CreateBodyRotated(_world, shapeType, s0, s1, s2, desc.PosX, desc.PosY, desc.PosZ, qx, qy, qz, qw, desc.Mass, desc.Friction, desc.Restitution);
		_bodyCount++;
		return index;
	}

	public void Step(float dt)
	{
		if (!_broadphaseOptimized)
		{
			Native.OptimizeBroadphase(_world);
			_broadphaseOptimized = true;
		}
		Native.Step(_world, dt);
	}

	public unsafe (float x, float y, float z) GetPosition(int bodyIndex)
	{
		float* pos = stackalloc float[3];
		Native.GetPosition(_world, bodyIndex, pos);
		return (pos[0], pos[1], pos[2]);
	}

	public unsafe (float x, float y, float z, float w) GetRotation(int bodyIndex)
	{
		float* rot = stackalloc float[4];
		Native.GetRotation(_world, bodyIndex, rot);
		return (rot[0], rot[1], rot[2], rot[3]);
	}

	public double GetLastStepTimeMs() => Native.GetStepTime(_world) * 1000.0;

	public int GetActiveBodyCount() => Native.GetActiveCount(_world);

	public void SetVelocity(int bodyIndex, float vx, float vy, float vz) => Native.SetVelocity(_world, bodyIndex, vx, vy, vz);

	public bool IsBodyActive(int bodyIndex) => Native.IsBodyActive(_world, bodyIndex) != 0;

	struct DragState { public int AnchorBody, ConstraintId; }
	readonly Dictionary<int, DragState> _drags = new();
	int _nextDragId = 1;

	public unsafe int BeginDrag(int bodyIndex, float hitX, float hitY, float hitZ)
	{
		int anchor = Native.CreateKinematicBody(_world, hitX, hitY, hitZ);

		// Compute body-local offset: inverse-rotate (hit - bodyPos)
		float* pos = stackalloc float[3];
		Native.GetPosition(_world, bodyIndex, pos);
		float* rot = stackalloc float[4];
		Native.GetRotation(_world, bodyIndex, rot);
		float qx = rot[0], qy = rot[1], qz = rot[2], qw = rot[3];
		float dx = hitX - pos[0], dy = hitY - pos[1], dz = hitZ - pos[2];
		// q^-1 * d: conjugate quaternion rotation
		float tx = 2 * (qy * dz - qz * dy), ty = 2 * (qz * dx - qx * dz), tz = 2 * (qx * dy - qy * dx);
		float lx = dx - qw * tx - (qy * tz - qz * ty);
		float ly = dy - qw * ty - (qz * tx - qx * tz);
		float lz = dz - qw * tz - (qx * ty - qy * tx);

		int constraint = Native.CreateSpringConstraint(_world, anchor, bodyIndex, 0, 0, 0, lx, ly, lz, 5.0f, 0.7f);

		int id = _nextDragId++;
		_drags[id] = new DragState { AnchorBody = anchor, ConstraintId = constraint };
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
		Native.RemoveConstraint(_world, state.ConstraintId);
		Native.RemoveBody(_world, state.AnchorBody);
	}

	public void AddDistanceJoint(int bodyA, int bodyB, float localAx, float localAy, float localAz, float localBx, float localBy, float localBz, float restLength)
	{
		Native.CreateDistanceJoint(_world, bodyA, bodyB, localAx, localAy, localAz, localBx, localBy, localBz, restLength);
	}

	public void Dispose()
	{
		if (_world != 0)
		{
			// Clean up any active drags before destroying world
			foreach (var state in _drags.Values)
			{
				Native.RemoveConstraint(_world, state.ConstraintId);
				Native.RemoveBody(_world, state.AnchorBody);
			}
			_drags.Clear();
			Native.DestroyWorld(_world);
			_world = 0;
		}
	}
}
