using Testbed;

namespace Testbed.Visual;

public struct BodyInfo
{
	public int Index;
	public ShapeType Shape;
	public float HX, HY, HZ;
	public float Radius;
	public float HalfHeight;
	public bool IsStatic;
}

public class EngineSlot : IDisposable
{
	public IPhysicsAdapter Adapter;
	public string Name;
	public float OffsetX;
	public Raylib_cs.Color DynColor;
	public List<BodyInfo> Bodies = new();

	public EngineSlot(IPhysicsAdapter adapter, string name, float offsetX, Raylib_cs.Color dynColor)
	{
		Adapter = adapter;
		Name = name;
		OffsetX = offsetX;
		DynColor = dynColor;
	}

	public int AddBody(BodyDesc desc)
	{
		int idx = Adapter.AddBody(desc);
		Bodies.Add(new BodyInfo
		{
			Index = idx,
			Shape = desc.Shape,
			HX = desc.HalfExtentX, HY = desc.HalfExtentY, HZ = desc.HalfExtentZ,
			Radius = desc.Radius,
			HalfHeight = desc.HalfHeight,
			IsStatic = desc.Mass <= 0,
		});
		return idx;
	}

	public void Dispose() => Adapter.Dispose();
}

public delegate void SceneSetup(EngineSlot slot);

public static class VisualScenes
{
	static void AddGround(EngineSlot slot)
	{
		slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = 0, PosY = -0.5f, PosZ = 0,
			HalfExtentX = 20, HalfExtentY = 0.5f, HalfExtentZ = 20,
			Mass = 0, Friction = 0.5f,
		});
	}

	public static void BoxStack50(EngineSlot slot)
	{
		AddGround(slot);
		for (int i = 0; i < 50; i++)
		{
			slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Box,
				PosX = 0, PosY = 0.5f + i, PosZ = 0,
				HalfExtentX = 0.5f, HalfExtentY = 0.5f, HalfExtentZ = 0.5f,
				Mass = 1, Friction = 0.5f,
			});
		}
	}

	public static void Pyramid15(EngineSlot slot)
	{
		AddGround(slot);
		int baseSize = 15;
		for (int row = 0; row < baseSize; row++)
		{
			int count = baseSize - row;
			float startX = -(count - 1) * 0.5f;
			for (int i = 0; i < count; i++)
			{
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Box,
					PosX = startX + i, PosY = 0.5f + row, PosZ = 0,
					HalfExtentX = 0.5f, HalfExtentY = 0.5f, HalfExtentZ = 0.5f,
					Mass = 1, Friction = 0.6f,
				});
			}
		}
	}

	public static void SphereDrop(EngineSlot slot)
	{
		AddGround(slot);
		int n = 8;
		float spacing = 1.5f;
		float startX = -(n - 1) * spacing * 0.5f;
		for (int x = 0; x < n; x++)
		{
			for (int z = 0; z < n; z++)
			{
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Sphere,
					PosX = startX + x * spacing,
					PosY = 5 + (x * n + z) * 0.2f,
					PosZ = startX + z * spacing,
					Radius = 0.5f,
					Mass = 1, Friction = 0.3f, Restitution = 0.6f,
				});
			}
		}
	}

	public static void Dominos(EngineSlot slot)
	{
		AddGround(slot);
		int count = 20;
		float spacing = 1.2f;
		float startX = -(count - 1) * spacing * 0.5f;
		for (int i = 0; i < count; i++)
		{
			slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Box,
				PosX = startX + i * spacing, PosY = 1.0f, PosZ = 0,
				HalfExtentX = 0.1f, HalfExtentY = 1.0f, HalfExtentZ = 0.5f,
				Mass = 1, Friction = 0.4f,
			});
		}
		// Pusher ball — launch into first domino
		int ball = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Sphere,
			PosX = startX - 2, PosY = 1.0f, PosZ = 0,
			Radius = 0.5f,
			Mass = 5, Friction = 0.3f,
		});
		slot.Adapter.SetVelocity(ball, 8, 0, 0);
	}

	public static void BoxWall(EngineSlot slot)
	{
		AddGround(slot);
		int cols = 10, rows = 10;
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				float offset = (r % 2 == 0) ? 0 : 0.5f;
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Box,
					PosX = (c - cols * 0.5f) + offset, PosY = 0.5f + r, PosZ = 0,
					HalfExtentX = 0.5f, HalfExtentY = 0.5f, HalfExtentZ = 0.5f,
					Mass = 1, Friction = 0.5f,
				});
			}
		}
		// Wrecking sphere — launch into wall
		int wrecker = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Sphere,
			PosX = -10, PosY = 5, PosZ = 0,
			Radius = 1.5f,
			Mass = 50, Friction = 0.3f,
		});
		slot.Adapter.SetVelocity(wrecker, 15, 0, 0);
	}
}
