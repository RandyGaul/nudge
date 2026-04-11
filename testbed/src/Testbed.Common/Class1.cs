namespace Testbed;

public enum ShapeType { Box, Sphere, Capsule }

public struct BodyDesc
{
	public ShapeType Shape;
	public float PosX, PosY, PosZ;
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
}
