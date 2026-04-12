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

	public void AddDistanceJoint(int bodyA, int bodyB, float localAx, float localAy, float localAz, float localBx, float localBy, float localBz, float restLength)
	{
		Native.CreateDistanceJoint(_world, bodyA, bodyB, localAx, localAy, localAz, localBx, localBy, localBz, restLength);
	}

	public void Dispose()
	{
		if (_world != 0)
		{
			Native.DestroyWorld(_world);
			_world = 0;
		}
	}
}
