using System.Numerics;
using Raylib_cs;
using Testbed;
using Testbed.Visual;
using Testbed.Nudge;
using Testbed.Bepu;
using Testbed.Jolt;

var scenes = new (string name, SceneSetup setup)[]
{
	("1: Box Stack",      VisualScenes.BoxStack50),
	("2: Pyramid",        VisualScenes.Pyramid15),
	("3: Sphere Drop",    VisualScenes.SphereDrop),
	("4: Dominos",        VisualScenes.Dominos),
	("5: Box Wall",       VisualScenes.BoxWall),
	("6: Friction Ramp",  VisualScenes.FrictionRamp),
	("7: Pendulum Chain", VisualScenes.PendulumChain),
	("8: Avalanche",      VisualScenes.ShapeAvalanche),
	("9: Bowling",        VisualScenes.Bowling),
	("A: Bouncy Pit",     VisualScenes.BouncyPit),
	("B: Funnel",         VisualScenes.Funnel),
};

const float ENGINE_SPACING = 30.0f;
const float DT = 1.0f / 60.0f;

Raylib.InitWindow(1600, 900, "Physics Testbed -- Nudge vs Bepu vs Jolt");
Raylib.SetTargetFPS(60);

// Maya-style orbit camera (matches nudge app): left-drag orbit, middle-drag pan, scroll zoom.
float camYaw = 0.0f, camPitch = -0.45f, camDist = 60.0f;
var camFocus = new Vector3(0, 10, 0);
var camera = new Camera3D
{
	Position = new Vector3(0, 30, 60),
	Target = new Vector3(0, 5, 0),
	Up = Vector3.UnitY,
	FovY = 45,
	Projection = CameraProjection.Perspective,
};

Color nudgeColor = new(80, 120, 220, 255);
Color bepuColor = new(80, 200, 100, 255);
Color joltColor = new(220, 140, 50, 255);
Color staticColor = new(160, 160, 160, 255);
Color wireColor = new(40, 40, 40, 255);

EngineSlot[] slots = [];
int currentScene = -1;
bool paused = false;

// Mouse picking state
int[] dragHandles = new int[3] { -1, -1, -1 };
float dragDist = 0;
Vector3 dragTarget = Vector3.Zero;

void LoadScene(int index)
{
	foreach (var s in slots) s.Dispose();

	slots =
	[
		new(new NudgeAdapter(), "Nudge", -ENGINE_SPACING, nudgeColor),
		new(new BepuAdapter(), "Bepu", 0, bepuColor),
		new(new JoltAdapter(), "Jolt", ENGINE_SPACING, joltColor),
	];

	foreach (var slot in slots)
	{
		slot.Adapter.CreateWorld(0, -9.81f, 0);
		scenes[index].setup(slot);
	}

	currentScene = index;
	dragHandles = new int[3] { -1, -1, -1 };
}

LoadScene(0);

while (!Raylib.WindowShouldClose())
{
	// Scene switching: 1-9 for first 9, A/B for 10-11, arrows to cycle
	for (int i = 0; i < Math.Min(scenes.Length, 9); i++)
	{
		if (Raylib.IsKeyPressed(KeyboardKey.One + i))
			LoadScene(i);
	}
	if (scenes.Length > 9 && Raylib.IsKeyPressed(KeyboardKey.A))
		LoadScene(9);
	if (scenes.Length > 10 && Raylib.IsKeyPressed(KeyboardKey.B))
		LoadScene(10);
	if (Raylib.IsKeyPressed(KeyboardKey.Right))
		LoadScene((currentScene + 1) % scenes.Length);
	if (Raylib.IsKeyPressed(KeyboardKey.Left))
		LoadScene((currentScene - 1 + scenes.Length) % scenes.Length);
	if (Raylib.IsKeyPressed(KeyboardKey.Space))
		paused = !paused;

	// Step physics
	if (!paused)
	{
		foreach (var slot in slots)
			slot.Adapter.Step(DT);
	}

	// Render
	Raylib.BeginDrawing();
	Raylib.ClearBackground(new Color(30, 30, 35, 255));

	// --- Mouse picking: right-click drag to interact with bodies ---
	if (Raylib.IsMouseButtonPressed(MouseButton.Right))
	{
		var mousePos = Raylib.GetMousePosition();
		var ray = Raylib.GetScreenToWorldRay(mousePos, camera);

		for (int si = 0; si < slots.Length; si++)
		{
			var slot = slots[si];
			float bestT = float.MaxValue;
			int bestBody = -1;

			foreach (var body in slot.Bodies)
			{
				if (body.IsStatic) continue;
				var (px, py, pz) = slot.Adapter.GetPosition(body.Index);
				var bodyPos = new Vector3(px + slot.OffsetX, py, pz);

				float radius = body.Shape switch
				{
					ShapeType.Box => MathF.Sqrt(body.HX * body.HX + body.HY * body.HY + body.HZ * body.HZ),
					ShapeType.Sphere => body.Radius,
					ShapeType.Capsule => body.HalfHeight + body.Radius,
					_ => 1.0f,
				};

				var oc = ray.Position - bodyPos;
				float b = Vector3.Dot(oc, ray.Direction);
				float c = Vector3.Dot(oc, oc) - radius * radius;
				float disc = b * b - c;
				if (disc < 0) continue;
				float t = -b - MathF.Sqrt(disc);
				if (t < 0) t = -b + MathF.Sqrt(disc);
				if (t >= 0 && t < bestT)
				{
					bestT = t;
					bestBody = body.Index;
				}
			}

			if (bestBody >= 0)
			{
				var hitPoint = ray.Position + ray.Direction * bestT;
				float hitX = hitPoint.X - slot.OffsetX;
				float hitY = hitPoint.Y;
				float hitZ = hitPoint.Z;
				dragHandles[si] = slot.Adapter.BeginDrag(bestBody, hitX, hitY, hitZ);
				dragDist = bestT;
			}
		}
	}

	if (Raylib.IsMouseButtonDown(MouseButton.Right))
	{
		var mousePos = Raylib.GetMousePosition();
		var ray = Raylib.GetScreenToWorldRay(mousePos, camera);
		dragTarget = ray.Position + ray.Direction * dragDist;

		for (int si = 0; si < slots.Length; si++)
		{
			if (dragHandles[si] < 0) continue;
			float tx = dragTarget.X - slots[si].OffsetX;
			float ty = dragTarget.Y;
			float tz = dragTarget.Z;
			slots[si].Adapter.UpdateDrag(dragHandles[si], tx, ty, tz);
		}
	}

	if (Raylib.IsMouseButtonReleased(MouseButton.Right))
	{
		for (int si = 0; si < slots.Length; si++)
		{
			if (dragHandles[si] < 0) continue;
			slots[si].Adapter.EndDrag(dragHandles[si]);
			dragHandles[si] = -1;
		}
	}

	// Maya-style orbit camera
	if (Raylib.IsMouseButtonDown(MouseButton.Left) && !Raylib.IsKeyDown(KeyboardKey.LeftControl))
	{
		var delta = Raylib.GetMouseDelta();
		camYaw += -delta.X * 0.005f;
		camPitch += -delta.Y * 0.005f;
		camPitch = Math.Clamp(camPitch, -1.5f, 1.5f);
	}
	if (Raylib.IsMouseButtonDown(MouseButton.Middle))
	{
		var delta = Raylib.GetMouseDelta();
		float sens = 0.005f * camDist;
		float sy = MathF.Sin(camYaw), cy = MathF.Cos(camYaw);
		camFocus += new Vector3(-cy * delta.X * sens, delta.Y * sens, sy * delta.X * sens);
	}
	float wheel = Raylib.GetMouseWheelMove();
	if (wheel != 0) { camDist *= 1.0f - wheel * 0.1f; camDist = Math.Clamp(camDist, 0.5f, 200.0f); }

	// Rebuild camera from yaw/pitch/dist
	float cpf = MathF.Cos(camPitch), spf = MathF.Sin(camPitch);
	float cyf = MathF.Cos(camYaw), syf = MathF.Sin(camYaw);
	var offset = new Vector3(syf * cpf, -spf, cyf * cpf) * camDist;
	camera.Position = camFocus + offset;
	camera.Target = camFocus;

	Raylib.BeginMode3D(camera);

	// Ground grid per engine
	foreach (var slot in slots)
	{
		for (int i = -10; i <= 10; i++)
		{
			Raylib.DrawLine3D(new Vector3(slot.OffsetX - 10, 0, i), new Vector3(slot.OffsetX + 10, 0, i), new Color(60, 60, 60, 255));
			Raylib.DrawLine3D(new Vector3(slot.OffsetX + i, 0, -10), new Vector3(slot.OffsetX + i, 0, 10), new Color(60, 60, 60, 255));
		}
	}

	// Draw bodies
	foreach (var slot in slots)
	{
		foreach (var body in slot.Bodies)
		{
			var (px, py, pz) = slot.Adapter.GetPosition(body.Index);
			var (qx, qy, qz, qw) = slot.Adapter.GetRotation(body.Index);
			var pos = new Vector3(px + slot.OffsetX, py, pz);
			bool active = body.IsStatic ? false : slot.Adapter.IsBodyActive(body.Index);
			Color color;
			if (body.IsStatic) color = staticColor;
			else if (!active) color = new Color((byte)(slot.DynColor.R / 3), (byte)(slot.DynColor.G / 3), (byte)(slot.DynColor.B / 3), (byte)255); // dim = sleeping
			else color = slot.DynColor;

			switch (body.Shape)
			{
				case ShapeType.Box:
					DrawRotatedCube(pos, qx, qy, qz, qw, body.HX * 2, body.HY * 2, body.HZ * 2, color, wireColor);
					break;
				case ShapeType.Sphere:
					Raylib.DrawSphere(pos, body.Radius, color);
					break;
				case ShapeType.Capsule:
					// Rotate the capsule endpoints by the body quaternion
					var localUp = RotateVec(new Vector3(0, body.HalfHeight, 0), qx, qy, qz, qw);
					var top = pos + localUp;
					var bot = pos - localUp;
					Raylib.DrawCapsule(bot, top, body.Radius, 8, 4, color);
					break;
			}
		}
	}

	// Draw drag indicator
	if (Raylib.IsMouseButtonDown(MouseButton.Right))
	{
		bool anyDrag = false;
		for (int si = 0; si < slots.Length; si++)
			if (dragHandles[si] >= 0) { anyDrag = true; break; }
		if (anyDrag)
			Raylib.DrawSphere(dragTarget, 0.15f, Color.Yellow);
	}

	Raylib.EndMode3D();

	// HUD
	int y = 10;
	Raylib.DrawText($"Scene: {scenes[currentScene].name}", 10, y, 20, Color.White);
	y += 25;
	Raylib.DrawText("Keys: 1-9/A/B switch scene, LEFT/RIGHT cycle, SPACE pause | LMB orbit, RMB drag body", 10, y, 16, Color.Gray);
	y += 25;

	for (int i = 0; i < slots.Length; i++)
	{
		var slot = slots[i];
		double ms = slot.Adapter.GetLastStepTimeMs();
		int awake = slot.Adapter.GetActiveBodyCount();
		int total = slot.Bodies.Count;
		Raylib.DrawText($"{slot.Name}: {ms:F2} ms  ({awake}/{total} awake)", 10, y, 18, slot.DynColor);
		y += 22;
		string breakdown = slot.Adapter.GetPerfBreakdown();
		if (breakdown.Length > 0)
		{
			Raylib.DrawText($"  {breakdown}", 10, y, 14, Color.Gray);
			y += 18;
		}
	}

	if (paused)
		Raylib.DrawText("PAUSED", Raylib.GetScreenWidth() / 2 - 50, Raylib.GetScreenHeight() / 2, 30, Color.Red);

	Raylib.DrawFPS(Raylib.GetScreenWidth() - 90, 10);
	Raylib.EndDrawing();
}

foreach (var s in slots) s.Dispose();
Raylib.CloseWindow();

// ---------------------------------------------------------------------------
// Helpers.

static Vector3 RotateVec(Vector3 v, float qx, float qy, float qz, float qw)
{
	// q * v * q^-1 (quaternion rotation)
	float tx = 2 * (qy * v.Z - qz * v.Y);
	float ty = 2 * (qz * v.X - qx * v.Z);
	float tz = 2 * (qx * v.Y - qy * v.X);
	return new Vector3(
		v.X + qw * tx + (qy * tz - qz * ty),
		v.Y + qw * ty + (qz * tx - qx * tz),
		v.Z + qw * tz + (qx * ty - qy * tx));
}

static unsafe void DrawRotatedCube(Vector3 pos, float qx, float qy, float qz, float qw, float sx, float sy, float sz, Color fill, Color wire)
{
	// Build 4x4 column-major matrix from quaternion + position + scale
	float xx = qx * qx, yy = qy * qy, zz = qz * qz;
	float xy = qx * qy, xz = qx * qz, yz = qy * qz;
	float wx = qw * qx, wy = qw * qy, wz = qw * qz;

	float* m = stackalloc float[16];
	m[0]  = (1 - 2*(yy+zz)) * sx;  m[1]  = (2*(xy+wz)) * sx;      m[2]  = (2*(xz-wy)) * sx;      m[3]  = 0;
	m[4]  = (2*(xy-wz)) * sy;      m[5]  = (1 - 2*(xx+zz)) * sy;  m[6]  = (2*(yz+wx)) * sy;      m[7]  = 0;
	m[8]  = (2*(xz+wy)) * sz;      m[9]  = (2*(yz-wx)) * sz;      m[10] = (1 - 2*(xx+yy)) * sz;  m[11] = 0;
	m[12] = pos.X;                  m[13] = pos.Y;                  m[14] = pos.Z;                  m[15] = 1;

	Rlgl.PushMatrix();
	Rlgl.MultMatrixf(m);
	Raylib.DrawCubeV(Vector3.Zero, Vector3.One, fill);
	Raylib.DrawCubeWiresV(Vector3.Zero, Vector3.One, wire);
	Rlgl.PopMatrix();
}
