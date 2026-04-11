namespace Testbed;

public static class Scenarios
{
	public static int StackBoxes(IPhysicsAdapter adapter, int count)
	{
		adapter.AddBody(new BodyDesc { Shape = ShapeType.Box, PosX = 0, PosY = -0.5f, PosZ = 0, HalfExtentX = 50, HalfExtentY = 0.5f, HalfExtentZ = 50, Mass = 0, Friction = 0.5f });
		float half = 0.5f;
		for (int i = 0; i < count; i++)
			adapter.AddBody(new BodyDesc { Shape = ShapeType.Box, PosX = 0, PosY = half + i, PosZ = 0, HalfExtentX = half, HalfExtentY = half, HalfExtentZ = half, Mass = 1.0f, Friction = 0.5f });
		return count + 1;
	}

	public static int SphereDrop(IPhysicsAdapter adapter, int countPerSide)
	{
		adapter.AddBody(new BodyDesc { Shape = ShapeType.Box, PosX = 0, PosY = -0.5f, PosZ = 0, HalfExtentX = 50, HalfExtentY = 0.5f, HalfExtentZ = 50, Mass = 0, Friction = 0.5f });
		int total = 0;
		float spacing = 2.0f;
		float start = -(countPerSide - 1) * spacing * 0.5f;
		for (int x = 0; x < countPerSide; x++)
			for (int z = 0; z < countPerSide; z++)
			{
				adapter.AddBody(new BodyDesc { Shape = ShapeType.Sphere, PosX = start + x * spacing, PosY = 5.0f + (x * countPerSide + z) * 0.1f, PosZ = start + z * spacing, Radius = 0.5f, Mass = 1.0f, Friction = 0.3f, Restitution = 0.5f });
				total++;
			}
		return total + 1;
	}

	public static int PyramidBoxes(IPhysicsAdapter adapter, int baseSize)
	{
		adapter.AddBody(new BodyDesc { Shape = ShapeType.Box, PosX = 0, PosY = -0.5f, PosZ = 0, HalfExtentX = 50, HalfExtentY = 0.5f, HalfExtentZ = 50, Mass = 0, Friction = 0.6f });
		int total = 0;
		float half = 0.5f;
		for (int row = 0; row < baseSize; row++)
		{
			int count = baseSize - row;
			float startX = -(count - 1) * 0.5f;
			for (int i = 0; i < count; i++)
			{
				adapter.AddBody(new BodyDesc { Shape = ShapeType.Box, PosX = startX + i, PosY = half + row, PosZ = 0, HalfExtentX = half, HalfExtentY = half, HalfExtentZ = half, Mass = 1.0f, Friction = 0.6f });
				total++;
			}
		}
		return total + 1;
	}
}
