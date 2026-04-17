using System.Numerics;
using System.Runtime.InteropServices;
using Raylib_cs;
using Testbed;
using Testbed.Visual;
using Testbed.Nudge;
using Testbed.Bepu;
using Testbed.Jolt;

// ---------------------------------------------------------------------------
// Shadow mapping shaders (GLSL 330).

// Depth-only shader for shadow pass (no shadow map sampling -- avoids read-while-write).
const string DEPTH_VS = @"#version 330
in vec3 vertexPosition;
uniform mat4 mvp;
void main() { gl_Position = mvp * vec4(vertexPosition, 1.0); }
";
const string DEPTH_FS = @"#version 330
out vec4 finalColor;
void main() { finalColor = vec4(0.0); }
";

const string SHADOW_VS = @"#version 330
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;
in vec4 vertexColor;
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;
out vec3 fragPosition;
out vec2 fragTexCoord;
out vec3 fragNormal;
void main()
{
    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
    fragTexCoord = vertexTexCoord;
    fragNormal = normalize(vec3(matNormal * vec4(vertexNormal, 1.0)));
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
";

const string SHADOW_FS = @"#version 330
in vec3 fragPosition;
in vec2 fragTexCoord;
in vec3 fragNormal;

uniform sampler2D texture0;
uniform vec4 colDiffuse;
uniform vec3 lightDir;
uniform vec4 lightColor;
uniform vec4 ambient;
uniform vec3 viewPos;
uniform mat4 lightVP;
uniform sampler2D shadowMap;
uniform int shadowMapResolution;

out vec4 finalColor;

void main()
{
    vec4 texelColor = texture(texture0, fragTexCoord);
    vec3 normal = normalize(fragNormal);
    vec3 viewD = normalize(viewPos - fragPosition);
    vec3 l = -lightDir;

    float NdotL = max(dot(normal, l), 0.0);
    vec3 lightDot = lightColor.rgb * NdotL;

    float specCo = 0.0;
    if (NdotL > 0.0) specCo = pow(max(0.0, dot(viewD, reflect(-l, normal))), 16.0);

    finalColor = texelColor * ((colDiffuse + vec4(specCo)) * vec4(lightDot, 1.0));

    // Shadow calculation
    vec4 fragPosLS = lightVP * vec4(fragPosition, 1.0);
    fragPosLS.xyz /= fragPosLS.w;
    fragPosLS.xyz = (fragPosLS.xyz + 1.0) / 2.0;
    vec2 sampleCoords = fragPosLS.xy;
    float curDepth = fragPosLS.z;

    float shadow = 0.0;
    if (sampleCoords.x >= 0.0 && sampleCoords.x <= 1.0 &&
        sampleCoords.y >= 0.0 && sampleCoords.y <= 1.0 && curDepth <= 1.0)
    {
        float bias = max(0.0005 * (1.0 - dot(normal, l)), 0.00005) + 0.00002;
        int shadowCounter = 0;
        vec2 texelSize = vec2(1.0 / float(shadowMapResolution));
        for (int x = -1; x <= 1; x++)
        {
            for (int y = -1; y <= 1; y++)
            {
                float sampleDepth = texture(shadowMap, sampleCoords + texelSize * vec2(x, y)).r;
                if (curDepth - bias > sampleDepth) shadowCounter++;
            }
        }
        shadow = float(shadowCounter) / 9.0;
    }

    finalColor = mix(finalColor, vec4(0, 0, 0, 1), shadow);

    // Ambient (always applied)
    finalColor += texelColor * ambient * colDiffuse;

    // Gamma correction
    finalColor = pow(finalColor, vec4(1.0 / 2.2));
}
";

// ---------------------------------------------------------------------------
// Constants.

const int SHADOW_MAP_SIZE = 2048;
const float ENGINE_SPACING = 30.0f;
const float DT = 1.0f / 60.0f;

var scenes = new (string name, SceneSetup setup)[]
{
	("1: Box Stack",         VisualScenes.BoxStack50),
	("2: Pyramid",           VisualScenes.Pyramid15),
	("3: Dominos",           VisualScenes.Dominos),
	("4: Box Wall",          VisualScenes.BoxWall),
	("5: Pendulum Chain",    VisualScenes.PendulumChain),
	("6: Newton's Cradle",   VisualScenes.NewtonsCradle),
	("7: Suspension Bridge", VisualScenes.SuspensionBridge),
	("8: Tetherball",        VisualScenes.Tetherball),
};

Raylib.SetConfigFlags(ConfigFlags.Msaa4xHint);
Raylib.InitWindow(1600, 900, "Physics Testbed -- Nudge vs Bepu vs Jolt");
Raylib.SetTargetFPS(60);
Rlgl.EnableBackfaceCulling();

// ---------------------------------------------------------------------------
// Maya-style orbit camera.
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

// ---------------------------------------------------------------------------
// Shadow mapping setup.

var shadowShader = Raylib.LoadShaderFromMemory(SHADOW_VS, SHADOW_FS);
int lightDirLoc = Raylib.GetShaderLocation(shadowShader, "lightDir");
int lightColLoc = Raylib.GetShaderLocation(shadowShader, "lightColor");
int ambientLoc = Raylib.GetShaderLocation(shadowShader, "ambient");
int lightVPLoc = Raylib.GetShaderLocation(shadowShader, "lightVP");
int shadowMapLoc = Raylib.GetShaderLocation(shadowShader, "shadowMap");
int shadowMapResLoc = Raylib.GetShaderLocation(shadowShader, "shadowMapResolution");
int viewPosLoc = Raylib.GetShaderLocation(shadowShader, "viewPos");

// Static uniforms.
var lightDir = Vector3.Normalize(new Vector3(0.6f, -0.75f, -0.4f));
Raylib.SetShaderValue(shadowShader, lightDirLoc, lightDir, ShaderUniformDataType.Vec3);
var lightColorV = new Vector4(1.0f, 1.0f, 1.0f, 1.0f);
Raylib.SetShaderValue(shadowShader, lightColLoc, lightColorV, ShaderUniformDataType.Vec4);
var ambientV = new Vector4(0.15f, 0.15f, 0.15f, 1.0f);
Raylib.SetShaderValue(shadowShader, ambientLoc, ambientV, ShaderUniformDataType.Vec4);
Raylib.SetShaderValue(shadowShader, shadowMapResLoc, SHADOW_MAP_SIZE, ShaderUniformDataType.Int);
Raylib.SetShaderValue(shadowShader, shadowMapLoc, 10, ShaderUniformDataType.Int);

// Depth-only shader for shadow pass.
var depthShader = Raylib.LoadShaderFromMemory(DEPTH_VS, DEPTH_FS);

// Depth-only shadow map FBO.
var shadowMap = LoadShadowmapRenderTexture(SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);

// Models used for shadowed rendering (boxes and spheres).
var cubeModel = Raylib.LoadModelFromMesh(Raylib.GenMeshCube(1.0f, 1.0f, 1.0f));
var sphereModel = Raylib.LoadModelFromMesh(Raylib.GenMeshSphere(1.0f, 16, 16));

unsafe
{
	cubeModel.Materials[0].Shader = shadowShader;
	sphereModel.Materials[0].Shader = shadowShader;
}

// Light camera (orthographic, directional light).
var lightCamera = new Camera3D
{
	Position = Vector3.Normalize(new Vector3(-0.6f, 0.75f, 0.4f)) * 50.0f + new Vector3(0, 5, 0),
	Target = new Vector3(0, 5, 0),
	Up = Vector3.UnitY,
	FovY = 120.0f,
	Projection = CameraProjection.Orthographic,
};

// ---------------------------------------------------------------------------
// State.

EngineSlot[] slots = [];
int currentScene = -1;
bool paused = false;
var capsuleModels = new Dictionary<(float, float), Model>();

// Mouse picking state.
int[] dragHandles = new int[3] { -1, -1, -1 };
float dragDist = 0;
Vector3 dragTarget = Vector3.Zero;

void LoadScene(int index)
{
	// End any active drags before destroying engines
	for (int si = 0; si < slots.Length && si < dragHandles.Length; si++)
	{
		if (dragHandles[si] >= 0) { slots[si].Adapter.EndDrag(dragHandles[si]); dragHandles[si] = -1; }
	}
	foreach (var s in slots) s.Dispose();

	// Free old capsule models.
	foreach (var m in capsuleModels.Values) Raylib.UnloadModel(m);
	capsuleModels.Clear();

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

	// Pre-build capsule models for unique (halfHeight, radius) pairs in this scene.
	foreach (var slot in slots)
	{
		foreach (var body in slot.Bodies)
		{
			if (body.Shape != ShapeType.Capsule) continue;
			var key = (body.HalfHeight, body.Radius);
			if (!capsuleModels.ContainsKey(key))
			{
				var mesh = GenCapsuleMesh(body.Radius, body.HalfHeight, 16, 8);
				var model = Raylib.LoadModelFromMesh(mesh);
				unsafe { model.Materials[0].Shader = shadowShader; }
				capsuleModels[key] = model;
			}
		}
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

	// =======================================================================
	// PASS 1: Render shadow map from light's perspective.
	// Use depth-only shader to avoid sampling the shadow map while writing it.
	// =======================================================================
	unsafe
	{
		cubeModel.Materials[0].Shader = depthShader;
		sphereModel.Materials[0].Shader = depthShader;
		foreach (var m in capsuleModels.Values) { var tmp = m; tmp.Materials[0].Shader = depthShader; }
	}
	Raylib.BeginTextureMode(shadowMap);
	Raylib.ClearBackground(Color.White);
	Raylib.BeginMode3D(lightCamera);
	var lightView = Rlgl.GetMatrixModelview();
	var lightProj = Rlgl.GetMatrixProjection();
	DrawBodiesModel(slots, cubeModel, sphereModel, capsuleModels, staticColor);
	Raylib.EndMode3D();
	Raylib.EndTextureMode();
	var lightViewProj = lightProj * lightView;
	unsafe
	{
		cubeModel.Materials[0].Shader = shadowShader;
		sphereModel.Materials[0].Shader = shadowShader;
		foreach (var m in capsuleModels.Values) { var tmp = m; tmp.Materials[0].Shader = shadowShader; }
	}

	// =======================================================================
	// PASS 2: Render scene with shadows.
	// =======================================================================
	Raylib.BeginDrawing();
	Raylib.ClearBackground(new Color(30, 30, 35, 255));

	// Set per-frame shadow uniforms.
	Raylib.SetShaderValueMatrix(shadowShader, lightVPLoc, lightViewProj);
	Raylib.SetShaderValue(shadowShader, viewPosLoc, camera.Position, ShaderUniformDataType.Vec3);

	// Bind shadow map depth texture to slot 10.
	Rlgl.EnableShader(shadowShader.Id);
	Rlgl.ActiveTextureSlot(10);
	Rlgl.EnableTexture(shadowMap.Depth.Id);

	Raylib.BeginMode3D(camera);

	// Ground grid per engine (unlit lines, drawn before shadowed geometry).
	foreach (var slot in slots)
	{
		for (int i = -10; i <= 10; i++)
		{
			Raylib.DrawLine3D(new Vector3(slot.OffsetX - 10, 0, i), new Vector3(slot.OffsetX + 10, 0, i), new Color(60, 60, 60, 255));
			Raylib.DrawLine3D(new Vector3(slot.OffsetX + i, 0, -10), new Vector3(slot.OffsetX + i, 0, 10), new Color(60, 60, 60, 255));
		}
	}

	// Shadowed bodies (all shapes via Models).
	DrawBodiesModel(slots, cubeModel, sphereModel, capsuleModels, staticColor);

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

// ---------------------------------------------------------------------------
// Cleanup.
foreach (var s in slots) s.Dispose();
foreach (var m in capsuleModels.Values) Raylib.UnloadModel(m);
Raylib.UnloadShader(depthShader);
Raylib.UnloadShader(shadowShader);
Raylib.UnloadModel(cubeModel);
Raylib.UnloadModel(sphereModel);
UnloadShadowmapRenderTexture(shadowMap);
Raylib.CloseWindow();

// ---------------------------------------------------------------------------
// Helpers.

static void DrawBodiesModel(EngineSlot[] slots, Model cubeModel, Model sphereModel, Dictionary<(float, float), Model> capsuleModels, Color staticColor)
{
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
			else if (!active) color = new Color((byte)(slot.DynColor.R / 3), (byte)(slot.DynColor.G / 3), (byte)(slot.DynColor.B / 3), (byte)255);
			else color = slot.DynColor;

			var (axis, angleDeg) = QuatToAxisAngle(qx, qy, qz, qw);

			switch (body.Shape)
			{
				case ShapeType.Box:
					Raylib.DrawModelEx(cubeModel, pos, axis, angleDeg, new Vector3(body.HX * 2, body.HY * 2, body.HZ * 2), color);
					break;
				case ShapeType.Sphere:
					Raylib.DrawModelEx(sphereModel, pos, axis, angleDeg, new Vector3(body.Radius, body.Radius, body.Radius), color);
					break;
				case ShapeType.Capsule:
					if (capsuleModels.TryGetValue((body.HalfHeight, body.Radius), out var capsModel))
						Raylib.DrawModelEx(capsModel, pos, axis, angleDeg, Vector3.One, color);
					break;
			}
		}
	}
}

static (Vector3 axis, float angleDeg) QuatToAxisAngle(float qx, float qy, float qz, float qw)
{
	if (qw < 0) { qx = -qx; qy = -qy; qz = -qz; qw = -qw; }
	float sinHalf = MathF.Sqrt(qx * qx + qy * qy + qz * qz);
	if (sinHalf < 1e-6f) return (Vector3.UnitY, 0.0f);
	float angle = 2.0f * MathF.Atan2(sinHalf, qw);
	return (new Vector3(qx / sinHalf, qy / sinHalf, qz / sinHalf), angle * (180.0f / MathF.PI));
}

static unsafe Mesh GenCapsuleMesh(float radius, float halfHeight, int slices, int rings)
{
	// Capsule = top hemisphere + cylinder body + bottom hemisphere, along Y axis.
	// rings = latitude divisions per hemisphere (not counting pole).
	int totalRings = 2 * rings;
	int vertsPerRing = slices + 1;
	int numVerts = 2 + totalRings * vertsPerRing;
	int numTris = 2 * slices + (totalRings - 1) * slices * 2;

	float* verts = (float*)NativeMemory.AllocZeroed((nuint)(numVerts * 3 * sizeof(float)));
	float* norms = (float*)NativeMemory.AllocZeroed((nuint)(numVerts * 3 * sizeof(float)));
	float* texcs = (float*)NativeMemory.AllocZeroed((nuint)(numVerts * 2 * sizeof(float)));
	ushort* inds = (ushort*)NativeMemory.AllocZeroed((nuint)(numTris * 3 * sizeof(ushort)));

	int vi = 0, ni = 0, ti = 0, ii = 0;

	// Top pole.
	verts[vi++] = 0; verts[vi++] = halfHeight + radius; verts[vi++] = 0;
	norms[ni++] = 0; norms[ni++] = 1; norms[ni++] = 0;
	texcs[ti++] = 0.5f; texcs[ti++] = 0;

	// Latitude rings.
	for (int r = 0; r < totalRings; r++)
	{
		float theta, yOff;
		if (r < rings)
		{
			theta = (r + 1) * MathF.PI / (2 * rings);
			yOff = halfHeight;
		}
		else
		{
			theta = MathF.PI / 2 + (r - rings) * MathF.PI / (2 * rings);
			yOff = -halfHeight;
		}

		float sinT = MathF.Sin(theta), cosT = MathF.Cos(theta);
		float y = yOff + radius * cosT;
		float rr = radius * sinT;

		for (int s = 0; s <= slices; s++)
		{
			float phi = s * 2 * MathF.PI / slices;
			float sinP = MathF.Sin(phi), cosP = MathF.Cos(phi);

			verts[vi++] = rr * cosP;
			verts[vi++] = y;
			verts[vi++] = rr * sinP;

			norms[ni++] = sinT * cosP;
			norms[ni++] = cosT;
			norms[ni++] = sinT * sinP;

			texcs[ti++] = (float)s / slices;
			texcs[ti++] = (float)(r + 1) / (totalRings + 1);
		}
	}

	// Bottom pole.
	verts[vi++] = 0; verts[vi++] = -(halfHeight + radius); verts[vi++] = 0;
	norms[ni++] = 0; norms[ni++] = -1; norms[ni++] = 0;
	texcs[ti++] = 0.5f; texcs[ti++] = 1;

	// Indices: top pole fan.
	for (int s = 0; s < slices; s++)
	{
		inds[ii++] = 0;
		inds[ii++] = (ushort)(1 + s);
		inds[ii++] = (ushort)(1 + s + 1);
	}

	// Indices: quad strips between rings.
	for (int r = 0; r < totalRings - 1; r++)
	{
		int cur = 1 + r * vertsPerRing;
		int next = 1 + (r + 1) * vertsPerRing;
		for (int s = 0; s < slices; s++)
		{
			inds[ii++] = (ushort)(cur + s);
			inds[ii++] = (ushort)(cur + s + 1);
			inds[ii++] = (ushort)(next + s);

			inds[ii++] = (ushort)(cur + s + 1);
			inds[ii++] = (ushort)(next + s + 1);
			inds[ii++] = (ushort)(next + s);
		}
	}

	// Indices: bottom pole fan.
	int lastRing = 1 + (totalRings - 1) * vertsPerRing;
	int botPole = numVerts - 1;
	for (int s = 0; s < slices; s++)
	{
		inds[ii++] = (ushort)botPole;
		inds[ii++] = (ushort)(lastRing + s + 1);
		inds[ii++] = (ushort)(lastRing + s);
	}

	var mesh = new Mesh { VertexCount = numVerts, TriangleCount = numTris, Vertices = verts, Normals = norms, TexCoords = texcs, Indices = inds };
	Raylib.UploadMesh(ref mesh, false);
	return mesh;
}

static RenderTexture2D LoadShadowmapRenderTexture(int width, int height)
{
	var target = new RenderTexture2D();
	target.Id = Rlgl.LoadFramebuffer();
	target.Texture.Width = width;
	target.Texture.Height = height;
	if (target.Id > 0)
	{
		Rlgl.EnableFramebuffer(target.Id);
		target.Depth.Id = Rlgl.LoadTextureDepth(width, height, false);
		target.Depth.Width = width;
		target.Depth.Height = height;
		target.Depth.Format = (PixelFormat)19;
		target.Depth.Mipmaps = 1;
		Rlgl.FramebufferAttach(target.Id, target.Depth.Id, FramebufferAttachType.Depth, FramebufferAttachTextureType.Texture2D, 0);
		if (Rlgl.FramebufferComplete(target.Id))
			Raylib.TraceLog(TraceLogLevel.Info, "SHADOW: Shadow map framebuffer created successfully");
		else
			Raylib.TraceLog(TraceLogLevel.Warning, "SHADOW: Shadow map framebuffer creation failed");
		Rlgl.DisableFramebuffer();
	}
	return target;
}

static void UnloadShadowmapRenderTexture(RenderTexture2D target)
{
	if (target.Id > 0)
		Rlgl.UnloadFramebuffer(target.Id);
}
