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

	// Classic desk-toy Newton's cradle: 5 heavy spheres hanging from a top bar
	// via paired cables (V-suspension so they swing only in X-Y plane). Leftmost
	// ball is pulled back to start the sequential momentum transfer.
	public static void NewtonsCradle(EngineSlot slot)
	{
		AddGround(slot);

		const int nBalls = 5;
		const float ballR = 0.45f;
		const float spacing = ballR * 2.02f;
		const float stringLen = 4.5f;
		const float zOffset = 1.0f;
		const float topBarY = 10f;

		int topBar = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = 0, PosY = topBarY, PosZ = 0,
			HalfExtentX = (nBalls - 1) * spacing * 0.5f + 0.5f,
			HalfExtentY = 0.2f,
			HalfExtentZ = zOffset + 0.2f,
			Mass = 0, Friction = 0.5f,
		});

		float cableLen = MathF.Sqrt(stringLen * stringLen + zOffset * zOffset);

		for (int i = 0; i < nBalls; i++)
		{
			float anchorX = (i - (nBalls - 1) * 0.5f) * spacing;
			float px, py;
			if (i == 0)
			{
				// Pull leftmost ball back ~60 degrees from vertical to start the cascade.
				float theta = -1.05f;
				px = anchorX + stringLen * MathF.Sin(theta);
				py = topBarY - stringLen * MathF.Cos(theta);
			}
			else
			{
				px = anchorX;
				py = topBarY - stringLen;
			}

			int ball = slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Sphere,
				PosX = px, PosY = py, PosZ = 0,
				Radius = ballR,
				Mass = 2, Friction = 0.05f, Restitution = 0.95f,
			});

			// V-shaped paired cables from top bar to ball (two anchors at ±zOffset).
			slot.Adapter.AddDistanceJoint(topBar, ball, anchorX, 0, +zOffset, 0, ballR, 0, cableLen);
			slot.Adapter.AddDistanceJoint(topBar, ball, anchorX, 0, -zOffset, 0, ballR, 0, cableLen);
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

	// Pendulum chain with a heavy bob, hurled at a tall box wall it actually
	// demolishes. Starts displaced AND with tangential velocity so the bob
	// reaches the wall with real momentum.
	public static void PendulumChain(EngineSlot slot)
	{
		AddGround(slot);

		// Geometry check: chain total length = (linkCount+1)*spacing = 9.
		// Anchor at (-1, 11) puts the bottom of the arc at (-1, 2) (above ground)
		// and the +θ reach touches the wall's front face (x=5.5) at y≈4.8.
		const float anchorY = 11f;
		const float anchorX = -1f;
		const int linkCount = 8;
		const float spacing = 1.0f;

		int anchor = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = anchorX, PosY = anchorY, PosZ = 0,
			HalfExtentX = 0.3f, HalfExtentY = 0.3f, HalfExtentZ = 0.3f,
			Mass = 0, Friction = 0.5f,
		});

		// Chain begins displaced ~45 degrees to the LEFT (so it swings down and right into the wall).
		const float launchAngle = -0.8f; // radians, negative = displaced -x side
		int[] links = new int[linkCount];
		for (int i = 0; i < linkCount; i++)
		{
			float t = (i + 1) * spacing;
			links[i] = slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Sphere,
				PosX = anchorX + t * MathF.Sin(launchAngle),
				PosY = anchorY - t * MathF.Cos(launchAngle),
				PosZ = 0,
				Radius = 0.22f,
				Mass = 1.0f, Friction = 0.3f,
			});
		}

		float tBall = (linkCount + 1) * spacing;
		int ball = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Sphere,
			PosX = anchorX + tBall * MathF.Sin(launchAngle),
			PosY = anchorY - tBall * MathF.Cos(launchAngle),
			PosZ = 0,
			Radius = 0.9f,
			Mass = 60f, Friction = 0.4f,
		});

		slot.Adapter.AddDistanceJoint(anchor, links[0], 0, 0, 0, 0, 0, 0, spacing);
		for (int i = 0; i < linkCount - 1; i++)
			slot.Adapter.AddDistanceJoint(links[i], links[i + 1], 0, 0, 0, 0, 0, 0, spacing);
		slot.Adapter.AddDistanceJoint(links[linkCount - 1], ball, 0, 0, 0, 0, 0, 0, spacing);

		// Target wall: 8 tall × 3 deep × 5 wide — dense enough to feel the hit.
		const float wallX = 6f;
		for (int r = 0; r < 8; r++)
			for (int cz = 0; cz < 3; cz++)
				for (int cy = 0; cy < 2; cy++)
					slot.AddBody(new BodyDesc
					{
						Shape = ShapeType.Box,
						PosX = wallX + cy * 1.0f, PosY = 0.5f + r, PosZ = (cz - 1) * 1.0f,
						HalfExtentX = 0.5f, HalfExtentY = 0.5f, HalfExtentZ = 0.5f,
						Mass = 1, Friction = 0.5f,
					});
	}

	// Plank suspension bridge between two towers. Each pair of adjacent planks
	// is connected by 2 short cables at the front and back edges (so the planks
	// can't twist apart). Heavy balls dropped onto it stress the chain.
	public static void SuspensionBridge(EngineSlot slot)
	{
		AddGround(slot);

		// Bigger gap + lighter planks keep the distance joints in their stable
		// regime. Tiny rest lengths with heavy masses cause wild oscillation.
		const int nPlanks = 10;
		const float plankHX = 0.5f, plankHY = 0.15f, plankHZ = 1.5f;
		const float gap = 0.2f;
		const float plankSpacing = plankHX * 2 + gap;
		const float bridgeY = 6f;
		const float towerHalfH = bridgeY * 0.5f;

		float bridgeHalfLen = nPlanks * plankSpacing * 0.5f;
		float leftTowerX = -bridgeHalfLen - 0.4f;
		float rightTowerX = +bridgeHalfLen + 0.4f;

		int leftTower = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = leftTowerX, PosY = towerHalfH, PosZ = 0,
			HalfExtentX = 0.4f, HalfExtentY = towerHalfH, HalfExtentZ = plankHZ + 0.3f,
			Mass = 0, Friction = 0.5f,
		});
		int rightTower = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = rightTowerX, PosY = towerHalfH, PosZ = 0,
			HalfExtentX = 0.4f, HalfExtentY = towerHalfH, HalfExtentZ = plankHZ + 0.3f,
			Mass = 0, Friction = 0.5f,
		});

		int[] planks = new int[nPlanks];
		for (int i = 0; i < nPlanks; i++)
		{
			float x = -bridgeHalfLen + plankSpacing * (i + 0.5f);
			planks[i] = slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Box,
				PosX = x, PosY = bridgeY, PosZ = 0,
				HalfExtentX = plankHX, HalfExtentY = plankHY, HalfExtentZ = plankHZ,
				Mass = 1f, Friction = 0.6f,
			});
		}

		// Adjacent plank edges joined at ±Z corners.
		for (int i = 0; i < nPlanks - 1; i++)
		{
			slot.Adapter.AddDistanceJoint(planks[i], planks[i + 1],
				+plankHX, 0, +plankHZ,
				-plankHX, 0, +plankHZ,
				gap);
			slot.Adapter.AddDistanceJoint(planks[i], planks[i + 1],
				+plankHX, 0, -plankHZ,
				-plankHX, 0, -plankHZ,
				gap);
		}

		// End planks tethered to towers (top-inside corner of tower to outer edge of plank).
		slot.Adapter.AddDistanceJoint(leftTower, planks[0],
			+0.4f, +towerHalfH, +plankHZ,
			-plankHX, 0, +plankHZ,
			gap);
		slot.Adapter.AddDistanceJoint(leftTower, planks[0],
			+0.4f, +towerHalfH, -plankHZ,
			-plankHX, 0, -plankHZ,
			gap);
		slot.Adapter.AddDistanceJoint(rightTower, planks[nPlanks - 1],
			-0.4f, +towerHalfH, +plankHZ,
			+plankHX, 0, +plankHZ,
			gap);
		slot.Adapter.AddDistanceJoint(rightTower, planks[nPlanks - 1],
			-0.4f, +towerHalfH, -plankHZ,
			+plankHX, 0, -plankHZ,
			gap);

		// Drop a few light balls along the bridge (heavier balls excite too much
		// oscillation in the plank chain).
		for (int i = 0; i < 3; i++)
		{
			slot.AddBody(new BodyDesc
			{
				Shape = ShapeType.Sphere,
				PosX = -2f + i * 2f, PosY = 9 + i * 0.3f, PosZ = 0,
				Radius = 0.4f,
				Mass = 1.5f, Friction = 0.4f, Restitution = 0.05f,
			});
		}
	}

	// Heavy ball tethered to a central pole, launched tangentially so it orbits
	// as a conical pendulum and scatters tall pillars arranged in two rings.
	// Pillar height matches the ball's steady-state orbit height.
	public static void Tetherball(EngineSlot slot)
	{
		AddGround(slot);

		const float poleH = 5f;
		const float poleHalfH = poleH * 0.5f;
		const float poleR = 0.25f;
		int pole = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Box,
			PosX = 0, PosY = poleHalfH, PosZ = 0,
			HalfExtentX = poleR, HalfExtentY = poleHalfH, HalfExtentZ = poleR,
			Mass = 0, Friction = 0.5f,
		});

		const float tetherLen = 3.5f;
		int ball = slot.AddBody(new BodyDesc
		{
			Shape = ShapeType.Sphere,
			PosX = tetherLen, PosY = poleH - 0.5f, PosZ = 0,
			Radius = 0.7f,
			Mass = 25f, Friction = 0.4f,
		});

		slot.Adapter.AddDistanceJoint(pole, ball, 0, +poleHalfH, 0, 0, 0, 0, tetherLen);

		// Tangential launch (Z axis) -- orbits clockwise viewed from above.
		// With tetherLen=3.5 and v=11, steady-state conical orbit sits at ~y=4
		// with horizontal radius ~3.4 -- right through the ring of pillars.
		slot.Adapter.SetVelocity(ball, 0, 0, 11f);

		// Two concentric rings of tall pillars. Heights span roughly y=0..5 so
		// the orbiting ball actually collides with them.
		int[] ringCounts = { 10, 16 };
		float[] ringR = { 2.7f, 4.2f };
		for (int ring = 0; ring < 2; ring++)
		{
			int count = ringCounts[ring];
			float r = ringR[ring];
			for (int i = 0; i < count; i++)
			{
				float angle = i * 2f * MathF.PI / count + ring * 0.15f;
				slot.AddBody(new BodyDesc
				{
					Shape = ShapeType.Box,
					PosX = r * MathF.Cos(angle),
					PosY = 2.5f,
					PosZ = r * MathF.Sin(angle),
					HalfExtentX = 0.35f, HalfExtentY = 2.5f, HalfExtentZ = 0.35f,
					Mass = 1.5f, Friction = 0.5f,
				});
			}
		}
	}

}
