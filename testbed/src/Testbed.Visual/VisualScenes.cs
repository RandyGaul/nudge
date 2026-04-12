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

	// ---- New scenes ----

	public static void FrictionRamp(EngineSlot slot)
	{
		AddGround(slot);

		// Ramp tilted 25 degrees around Z (left side high, right side low)
		float halfAngle = 12.5f * MathF.PI / 180f;
		slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = 0, PosY = 3.5f, PosZ = 0,
			RotZ = -MathF.Sin(halfAngle), RotW = MathF.Cos(halfAngle),
			HalfExtentX = 8, HalfExtentY = 0.25f, HalfExtentZ = 5,
			Mass = 0, Friction = 0.5f,
		});

		// 5 objects at top of ramp with different friction values
		float[] frictions = { 0.0f, 0.1f, 0.3f, 0.5f, 1.0f };
		for (int i = 0; i < frictions.Length; i++)
		{
			float z = -3.0f + i * 1.5f;
			slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Box,
				PosX = -5, PosY = 6.5f, PosZ = z,
				HalfExtentX = 0.4f, HalfExtentY = 0.4f, HalfExtentZ = 0.4f,
				Mass = 1, Friction = frictions[i],
			});
		}

		// A few spheres too (different shape, same friction spread)
		for (int i = 0; i < frictions.Length; i++)
		{
			float z = -3.0f + i * 1.5f;
			slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Sphere,
				PosX = -3.5f, PosY = 6.0f, PosZ = z,
				Radius = 0.35f,
				Mass = 1, Friction = frictions[i],
			});
		}
	}

	public static void PendulumChain(EngineSlot slot)
	{
		AddGround(slot);

		// Static anchor at the top
		int anchor = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = 0, PosY = 12, PosZ = 0,
			HalfExtentX = 0.3f, HalfExtentY = 0.3f, HalfExtentZ = 0.3f,
			Mass = 0, Friction = 0.5f,
		});

		// Chain of 8 small spheres
		float spacing = 1.0f;
		int linkCount = 8;
		int[] links = new int[linkCount];
		for (int i = 0; i < linkCount; i++)
		{
			links[i] = slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Sphere,
				PosX = 0, PosY = 12 - (i + 1) * spacing, PosZ = 0,
				Radius = 0.2f,
				Mass = 0.5f, Friction = 0.3f,
			});
		}

		// Heavy wrecking ball at the end
		int ball = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Sphere,
			PosX = 0, PosY = 12 - (linkCount + 1) * spacing, PosZ = 0,
			Radius = 0.7f,
			Mass = 12, Friction = 0.3f,
		});

		// Distance joints: anchor -> link0 -> link1 -> ... -> ball
		slot.Adapter.AddDistanceJoint(anchor, links[0], 0, 0, 0, 0, 0, 0, spacing);
		for (int i = 0; i < linkCount - 1; i++)
			slot.Adapter.AddDistanceJoint(links[i], links[i + 1], 0, 0, 0, 0, 0, 0, spacing);
		slot.Adapter.AddDistanceJoint(links[linkCount - 1], ball, 0, 0, 0, 0, 0, 0, spacing);

		// Swing the whole chain to the right
		foreach (int link in links)
			slot.Adapter.SetVelocity(link, 12, 0, 0);
		slot.Adapter.SetVelocity(ball, 12, 0, 0);

		// Target wall for the pendulum to smash
		for (int r = 0; r < 5; r++)
			for (int c = 0; c < 3; c++)
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Box,
					PosX = 7, PosY = 0.5f + r, PosZ = -1 + c,
					HalfExtentX = 0.5f, HalfExtentY = 0.5f, HalfExtentZ = 0.5f,
					Mass = 1, Friction = 0.5f,
				});
	}

	public static void ShapeAvalanche(EngineSlot slot)
	{
		AddGround(slot);

		// Staircase of platforms at descending heights
		slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = -7, PosY = 6, PosZ = 0,
			HalfExtentX = 3, HalfExtentY = 0.25f, HalfExtentZ = 5,
			Mass = 0, Friction = 0.4f,
		});
		slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = -2, PosY = 4, PosZ = 0,
			HalfExtentX = 3, HalfExtentY = 0.25f, HalfExtentZ = 5,
			Mass = 0, Friction = 0.4f,
		});
		slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = 3, PosY = 2, PosZ = 0,
			HalfExtentX = 3, HalfExtentY = 0.25f, HalfExtentZ = 5,
			Mass = 0, Friction = 0.4f,
		});

		// Drop mixed shapes onto the top platform
		int idx = 0;
		for (int x = 0; x < 3; x++)
		{
			for (int z = 0; z < 3; z++)
			{
				float px = -8.0f + x * 1.0f;
				float pz = -2.0f + z * 2.0f;

				// Boxes
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Box,
					PosX = px, PosY = 8 + idx * 0.6f, PosZ = pz,
					HalfExtentX = 0.3f, HalfExtentY = 0.3f, HalfExtentZ = 0.3f,
					Mass = 1, Friction = 0.4f,
				});

				// Spheres (offset slightly)
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Sphere,
					PosX = px + 0.5f, PosY = 9 + idx * 0.6f, PosZ = pz,
					Radius = 0.3f,
					Mass = 1, Friction = 0.3f, Restitution = 0.3f,
				});

				// Capsules
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Capsule,
					PosX = px + 0.25f, PosY = 10 + idx * 0.6f, PosZ = pz,
					HalfHeight = 0.3f, Radius = 0.15f,
					Mass = 1, Friction = 0.3f,
				});

				idx++;
			}
		}
	}

	public static void Bowling(EngineSlot slot)
	{
		AddGround(slot);

		// 10 pins in standard triangle formation (capsules standing upright)
		float pinHH = 0.35f, pinR = 0.1f;
		float pinY = pinHH + pinR;
		float rowSpacing = 1.1f;
		float pinSpacing = 0.7f;
		float baseX = 6;
		int[][] rows = { new[] { 0 }, new[] { -1, 1 }, new[] { -2, 0, 2 }, new[] { -3, -1, 1, 3 } };
		for (int r = 0; r < rows.Length; r++)
		{
			foreach (int col in rows[r])
			{
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Capsule,
					PosX = baseX + r * rowSpacing, PosY = pinY, PosZ = col * pinSpacing * 0.5f,
					HalfHeight = pinHH, Radius = pinR,
					Mass = 0.8f, Friction = 0.4f,
				});
			}
		}

		// Bowling ball
		int ball = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Sphere,
			PosX = -6, PosY = 0.5f, PosZ = 0,
			Radius = 0.5f,
			Mass = 7, Friction = 0.2f,
		});
		slot.Adapter.SetVelocity(ball, 14, 0, 0);
	}

	public static void BouncyPit(EngineSlot slot)
	{
		// Bouncy ground
		slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = 0, PosY = -0.5f, PosZ = 0,
			HalfExtentX = 20, HalfExtentY = 0.5f, HalfExtentZ = 20,
			Mass = 0, Friction = 0.3f, Restitution = 0.95f,
		});

		// Pit walls
		float wallH = 3, wallHalf = 0.25f, pitSize = 5;
		slot.AddBody(new BodyDesc { Shape = ShapeType.Box, PosX = 0, PosY = wallH, PosZ = -pitSize, HalfExtentX = pitSize, HalfExtentY = wallH, HalfExtentZ = wallHalf, Mass = 0, Friction = 0.3f, Restitution = 0.95f });
		slot.AddBody(new BodyDesc { Shape = ShapeType.Box, PosX = 0, PosY = wallH, PosZ = pitSize, HalfExtentX = pitSize, HalfExtentY = wallH, HalfExtentZ = wallHalf, Mass = 0, Friction = 0.3f, Restitution = 0.95f });
		slot.AddBody(new BodyDesc { Shape = ShapeType.Box, PosX = -pitSize, PosY = wallH, PosZ = 0, HalfExtentX = wallHalf, HalfExtentY = wallH, HalfExtentZ = pitSize, Mass = 0, Friction = 0.3f, Restitution = 0.95f });
		slot.AddBody(new BodyDesc { Shape = ShapeType.Box, PosX = pitSize, PosY = wallH, PosZ = 0, HalfExtentX = wallHalf, HalfExtentY = wallH, HalfExtentZ = pitSize, Mass = 0, Friction = 0.3f, Restitution = 0.95f });

		// Drop 25 bouncy spheres from varying heights
		for (int i = 0; i < 25; i++)
		{
			float angle = i * 2.4f; // golden angle spread
			float r = 2.0f * (i % 5) / 5.0f;
			slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Sphere,
				PosX = r * MathF.Cos(angle),
				PosY = 6 + i * 0.5f,
				PosZ = r * MathF.Sin(angle),
				Radius = 0.25f + (i % 3) * 0.15f,
				Mass = 1, Friction = 0.1f, Restitution = 0.95f,
			});
		}
	}

	public static void Funnel(EngineSlot slot)
	{
		AddGround(slot);

		// Two angled walls forming a V-funnel
		float wallAngle = 30.0f * MathF.PI / 180f;
		float halfAngle = wallAngle * 0.5f;

		// Left wall: tilted so top leans right (toward center)
		slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = -3.5f, PosY = 5, PosZ = 0,
			RotZ = -MathF.Sin(halfAngle), RotW = MathF.Cos(halfAngle),
			HalfExtentX = 0.2f, HalfExtentY = 4, HalfExtentZ = 5,
			Mass = 0, Friction = 0.3f,
		});
		// Right wall: tilted so top leans left (toward center)
		slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = 3.5f, PosY = 5, PosZ = 0,
			RotZ = MathF.Sin(halfAngle), RotW = MathF.Cos(halfAngle),
			HalfExtentX = 0.2f, HalfExtentY = 4, HalfExtentZ = 5,
			Mass = 0, Friction = 0.3f,
		});

		// Drop mixed objects from above the funnel
		for (int i = 0; i < 20; i++)
		{
			float px = -2.0f + (i % 5) * 1.0f;
			float pz = -3.0f + (i / 5) * 2.0f;
			float py = 12 + i * 0.4f;

			if (i % 3 == 0)
			{
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Box,
					PosX = px, PosY = py, PosZ = pz,
					HalfExtentX = 0.3f, HalfExtentY = 0.3f, HalfExtentZ = 0.3f,
					Mass = 1, Friction = 0.4f,
				});
			}
			else if (i % 3 == 1)
			{
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Sphere,
					PosX = px, PosY = py, PosZ = pz,
					Radius = 0.3f,
					Mass = 1, Friction = 0.3f, Restitution = 0.3f,
				});
			}
			else
			{
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Capsule,
					PosX = px, PosY = py, PosZ = pz,
					HalfHeight = 0.25f, Radius = 0.15f,
					Mass = 1, Friction = 0.3f,
				});
			}
		}
	}
}
