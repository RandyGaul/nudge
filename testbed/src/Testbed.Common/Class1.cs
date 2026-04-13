namespace Testbed;

public enum ShapeType { Box, Sphere, Capsule }

public struct BodyDesc
{
	public ShapeType Shape;
	public float PosX, PosY, PosZ;
	public float RotX, RotY, RotZ, RotW; // quaternion; all-zero = identity
	public float HalfExtentX, HalfExtentY, HalfExtentZ; // box
	public float Radius;     // sphere/capsule
	public float HalfHeight; // capsule
	public float Mass;       // 0 = static
	public float Friction;   // default 0.5
	public float Restitution;
}

public interface IPhysicsAdapter : IDisposable
{
	string Name { get; }
	void CreateWorld(float gravityX, float gravityY, float gravityZ);
	int AddBody(BodyDesc desc);
	void Step(float dt);
	(float x, float y, float z) GetPosition(int bodyIndex);
	(float x, float y, float z, float w) GetRotation(int bodyIndex);
	double GetLastStepTimeMs();
	int GetActiveBodyCount();
	bool IsBodyActive(int bodyIndex);
	void SetVelocity(int bodyIndex, float vx, float vy, float vz);
	string GetPerfBreakdown() => ""; // optional, override for details
	void AddDistanceJoint(int bodyA, int bodyB, float localAx, float localAy, float localAz, float localBx, float localBy, float localBz, float restLength) { } // optional

	// Mouse picking: create spring constraint to drag bodies.
	// Returns a drag handle (opaque int). -1 = failed.
	int BeginDrag(int bodyIndex, float hitX, float hitY, float hitZ) => -1;
	void UpdateDrag(int dragHandle, float targetX, float targetY, float targetZ) { }
	void EndDrag(int dragHandle) { }
}
