// render.c -- GL 3.3 core instanced mesh renderer

// -----------------------------------------------------------------------------
// GL extension loading. Windows gl.h is GL 1.1 only; load everything else.

typedef char GLchar;
typedef ptrdiff_t GLsizeiptr;
typedef ptrdiff_t GLintptr;

#define GL_ARRAY_BUFFER           0x8892
#define GL_ELEMENT_ARRAY_BUFFER   0x8893
#define GL_STATIC_DRAW            0x88E4
#define GL_DYNAMIC_DRAW           0x88E8
#define GL_FRAGMENT_SHADER        0x8B30
#define GL_VERTEX_SHADER          0x8B31
#define GL_COMPILE_STATUS         0x8B81
#define GL_LINK_STATUS            0x8B82
#define GL_INFO_LOG_LENGTH        0x8B84

#define GL_FUNCS \
	X(GLuint,  CreateShader,           GLenum type) \
	X(void,    ShaderSource,           GLuint shader, GLsizei count, const GLchar** string, const GLint* length) \
	X(void,    CompileShader,          GLuint shader) \
	X(void,    GetShaderiv,            GLuint shader, GLenum pname, GLint* params) \
	X(void,    GetShaderInfoLog,       GLuint shader, GLsizei bufSize, GLsizei* length, GLchar* infoLog) \
	X(void,    DeleteShader,           GLuint shader) \
	X(GLuint,  CreateProgram,          ) \
	X(void,    AttachShader,           GLuint program, GLuint shader) \
	X(void,    LinkProgram,            GLuint program) \
	X(void,    GetProgramiv,           GLuint program, GLenum pname, GLint* params) \
	X(void,    UseProgram,             GLuint program) \
	X(GLint,   GetUniformLocation,     GLuint program, const GLchar* name) \
	X(void,    UniformMatrix4fv,       GLint location, GLsizei count, GLboolean transpose, const GLfloat* value) \
	X(void,    Uniform3f,              GLint location, GLfloat v0, GLfloat v1, GLfloat v2) \
	X(void,    BlendFunc,              GLenum sfactor, GLenum dfactor) \
	X(void,    GenVertexArrays,        GLsizei n, GLuint* arrays) \
	X(void,    BindVertexArray,        GLuint array) \
	X(void,    GenBuffers,             GLsizei n, GLuint* buffers) \
	X(void,    BindBuffer,             GLenum target, GLuint buffer) \
	X(void,    BufferData,             GLenum target, GLsizeiptr size, const void* data, GLenum usage) \
	X(void,    BufferSubData,          GLenum target, GLintptr offset, GLsizeiptr size, const void* data) \
	X(void,    DeleteBuffers,          GLsizei n, const GLuint* buffers) \
	X(void,    VertexAttribPointer,    GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const void* pointer) \
	X(void,    EnableVertexAttribArray, GLuint index) \
	X(void,    VertexAttribDivisor,    GLuint index, GLuint divisor) \
	X(void,    DrawElementsInstanced,  GLenum mode, GLsizei count, GLenum type, const void* indices, GLsizei instancecount)

// Declare function pointers: gl_CreateShader, gl_ShaderSource, ...
#define X(ret, name, ...) typedef ret (APIENTRY *PFN_gl##name)(__VA_ARGS__); static PFN_gl##name gl_##name;
GL_FUNCS
#undef X

static void gl_load()
{
#define X(ret, name, ...) gl_##name = (PFN_gl##name)SDL_GL_GetProcAddress("gl" #name);
	GL_FUNCS
#undef X
}

// -----------------------------------------------------------------------------
// Shader compilation.

static const char* s_vert_src =
	"#version 330 core\n"
	"layout(location=0) in vec3 a_pos;\n"
	"layout(location=1) in vec3 a_normal;\n"
	"layout(location=2) in vec4 a_model_c0;\n"
	"layout(location=3) in vec4 a_model_c1;\n"
	"layout(location=4) in vec4 a_model_c2;\n"
	"layout(location=5) in vec4 a_model_c3;\n"
	"layout(location=6) in vec4 a_color;\n"
	"uniform mat4 u_vp;\n"
	"out vec3 v_normal;\n"
	"out vec4 v_color;\n"
	"void main() {\n"
	"    mat4 model = mat4(a_model_c0, a_model_c1, a_model_c2, a_model_c3);\n"
	"    gl_Position = u_vp * model * vec4(a_pos, 1.0);\n"
	"    v_normal = mat3(model) * a_normal;\n"
	"    v_color = a_color;\n"
	"}\n";

static const char* s_frag_src =
	"#version 330 core\n"
	"in vec3 v_normal;\n"
	"in vec4 v_color;\n"
	"out vec4 frag_color;\n"
	"uniform vec3 u_light_dir;\n"
	"uniform vec3 u_ambient;\n"
	"void main() {\n"
	"    vec3 n = normalize(v_normal);\n"
	"    float ndl = max(dot(n, u_light_dir), 0.0);\n"
	"    vec3 lit = v_color.rgb * (u_ambient + (1.0 - u_ambient) * ndl);\n"
	"    frag_color = vec4(lit * v_color.a, v_color.a);\n"
	"}\n";

static GLuint compile_shader(GLenum type, const char* src)
{
	GLuint s = gl_CreateShader(type);
	gl_ShaderSource(s, 1, &src, NULL);
	gl_CompileShader(s);
	GLint ok;
	gl_GetShaderiv(s, GL_COMPILE_STATUS, &ok);
	if (!ok) {
		char buf[512];
		gl_GetShaderInfoLog(s, sizeof(buf), NULL, buf);
		SDL_Log("Shader compile error: %s", buf);
	}
	return s;
}

static GLuint create_program(const char* vs_src, const char* fs_src)
{
	GLuint vs = compile_shader(GL_VERTEX_SHADER, vs_src);
	GLuint fs = compile_shader(GL_FRAGMENT_SHADER, fs_src);
	GLuint prog = gl_CreateProgram();
	gl_AttachShader(prog, vs);
	gl_AttachShader(prog, fs);
	gl_LinkProgram(prog);
	gl_DeleteShader(vs);
	gl_DeleteShader(fs);
	return prog;
}

// -----------------------------------------------------------------------------
// Mesh types.

typedef struct MeshVertex
{
	v3 pos;
	v3 normal;
} MeshVertex;

typedef struct Mesh
{
	GLuint vao;
	GLuint vbo;          // vertex data
	GLuint ebo;          // index data
	GLuint instance_vbo; // per-instance data
	int index_count;
} Mesh;

typedef enum MeshType
{
	MESH_BOX,
	MESH_SPHERE,
	MESH_BUILTIN_COUNT,
} MeshType;

#define MAX_MESH_TYPES 32

// Per-instance data: model matrix + color + opacity.
typedef struct RenderInstance
{
	mat4 model;
	float color[4]; // rgb + alpha (premultiplied)
} RenderInstance;

#define MAX_INSTANCES 4096

// -----------------------------------------------------------------------------
// Mesh VAO setup (shared logic).

static void mesh_upload(Mesh* mesh, MeshVertex* verts, int vert_count, uint16_t* indices, int idx_count)
{
	mesh->index_count = idx_count;

	gl_GenVertexArrays(1, &mesh->vao);
	gl_BindVertexArray(mesh->vao);

	// Vertex buffer
	gl_GenBuffers(1, &mesh->vbo);
	gl_BindBuffer(GL_ARRAY_BUFFER, mesh->vbo);
	gl_BufferData(GL_ARRAY_BUFFER, vert_count * sizeof(MeshVertex), verts, GL_STATIC_DRAW);

	// position
	gl_EnableVertexAttribArray(0);
	gl_VertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)0);
	// normal
	gl_EnableVertexAttribArray(1);
	gl_VertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)sizeof(v3));

	// Index buffer
	gl_GenBuffers(1, &mesh->ebo);
	gl_BindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->ebo);
	gl_BufferData(GL_ELEMENT_ARRAY_BUFFER, idx_count * sizeof(uint16_t), indices, GL_STATIC_DRAW);

	// Instance buffer (dynamic, filled each frame)
	gl_GenBuffers(1, &mesh->instance_vbo);
	gl_BindBuffer(GL_ARRAY_BUFFER, mesh->instance_vbo);
	gl_BufferData(GL_ARRAY_BUFFER, MAX_INSTANCES * sizeof(RenderInstance), NULL, GL_DYNAMIC_DRAW);

	// model matrix: 4 x vec4 at locations 2-5
	for (int i = 0; i < 4; i++) {
		gl_EnableVertexAttribArray(2 + i);
		gl_VertexAttribPointer(2 + i, 4, GL_FLOAT, GL_FALSE, sizeof(RenderInstance), (void*)(uintptr_t)(i * 16));
		gl_VertexAttribDivisor(2 + i, 1);
	}
	// color+alpha at location 6
	gl_EnableVertexAttribArray(6);
	gl_VertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(RenderInstance), (void*)offsetof(RenderInstance, color));
	gl_VertexAttribDivisor(6, 1);

	gl_BindVertexArray(0);
}

// -----------------------------------------------------------------------------
// Box mesh: 24 verts (4 per face for correct normals), 36 indices.

static void mesh_generate_box(Mesh* mesh)
{
	v3 p[8] = {
		{-1,-1,-1}, { 1,-1,-1}, { 1, 1,-1}, {-1, 1,-1},
		{-1,-1, 1}, { 1,-1, 1}, { 1, 1, 1}, {-1, 1, 1},
	};
	v3 n[6] = {
		{ 0, 0,-1}, { 0, 0, 1}, {-1, 0, 0}, { 1, 0, 0}, { 0,-1, 0}, { 0, 1, 0},
	};
	int faces[6][4] = {
		{0,3,2,1}, {4,5,6,7}, {0,4,7,3}, {1,2,6,5}, {0,1,5,4}, {3,7,6,2},
	};

	MeshVertex verts[24];
	uint16_t indices[36];
	for (int f = 0; f < 6; f++) {
		for (int v = 0; v < 4; v++) {
			verts[f*4+v].pos = p[faces[f][v]];
			verts[f*4+v].normal = n[f];
		}
		int b = f * 4;
		int i = f * 6;
		indices[i+0] = b+0; indices[i+1] = b+1; indices[i+2] = b+2;
		indices[i+3] = b+0; indices[i+4] = b+2; indices[i+5] = b+3;
	}
	mesh_upload(mesh, verts, 24, indices, 36);
}

// -----------------------------------------------------------------------------
// Geodesic sphere: recursive subdivision of octahedron.

static int sphere_midpoint(MeshVertex** verts, uint32_t** cache_keys, int** cache_vals, int i0, int i1)
{
	uint32_t lo = i0 < i1 ? i0 : i1;
	uint32_t hi = i0 < i1 ? i1 : i0;
	uint32_t key = (lo << 16) | hi;

	for (int i = 0; i < asize(*cache_keys); i++) {
		if ((*cache_keys)[i] == key) return (*cache_vals)[i];
	}

	v3 mid = norm(scale(add((*verts)[i0].pos, (*verts)[i1].pos), 0.5f));
	int idx = asize(*verts);
	MeshVertex mv = { mid, mid }; // on unit sphere, normal == position
	apush(*verts, mv);
	apush(*cache_keys, key);
	apush(*cache_vals, idx);
	return idx;
}

static void mesh_generate_sphere(Mesh* mesh)
{
	// Octahedron base
	MeshVertex* verts = NULL;
	MeshVertex oct_verts[6] = {
		{{ 0, 1, 0}, { 0, 1, 0}}, {{ 0,-1, 0}, { 0,-1, 0}},
		{{-1, 0, 0}, {-1, 0, 0}}, {{ 1, 0, 0}, { 1, 0, 0}},
		{{ 0, 0, 1}, { 0, 0, 1}}, {{ 0, 0,-1}, { 0, 0,-1}},
	};
	for (int i = 0; i < 6; i++) apush(verts, oct_verts[i]);

	uint16_t* tris = NULL;
	uint16_t oct_tris[] = {
		0,4,3, 0,3,5, 0,5,2, 0,2,4,
		1,3,4, 1,5,3, 1,2,5, 1,4,2,
	};
	for (int i = 0; i < 24; i++) apush(tris, oct_tris[i]);

	// Subdivide 3 times -> 512 triangles, ~258 verts
	for (int level = 0; level < 3; level++) {
		uint16_t* new_tris = NULL;
		uint32_t* cache_keys = NULL;
		int* cache_vals = NULL;

		for (int i = 0; i < asize(tris); i += 3) {
			int a = tris[i], b = tris[i+1], c = tris[i+2];
			int ab = sphere_midpoint(&verts, &cache_keys, &cache_vals, a, b);
			int bc = sphere_midpoint(&verts, &cache_keys, &cache_vals, b, c);
			int ca = sphere_midpoint(&verts, &cache_keys, &cache_vals, c, a);

			apush(new_tris, (uint16_t)a);  apush(new_tris, (uint16_t)ab); apush(new_tris, (uint16_t)ca);
			apush(new_tris, (uint16_t)b);  apush(new_tris, (uint16_t)bc); apush(new_tris, (uint16_t)ab);
			apush(new_tris, (uint16_t)c);  apush(new_tris, (uint16_t)ca); apush(new_tris, (uint16_t)bc);
			apush(new_tris, (uint16_t)ab); apush(new_tris, (uint16_t)bc); apush(new_tris, (uint16_t)ca);
		}

		afree(tris);
		afree(cache_keys);
		afree(cache_vals);
		tris = new_tris;
	}

	mesh_upload(mesh, verts, asize(verts), tris, asize(tris));
	afree(verts);
	afree(tris);
}

// -----------------------------------------------------------------------------
// Capsule mesh: two hemispheres + cylinder, baked radius and half_height.

static void mesh_generate_capsule(Mesh* mesh, float radius, float half_height)
{
	const int slices = 16;
	const int rings = 6;
	const float PI = 3.14159265f;

	MeshVertex* verts = NULL;
	uint16_t* tris = NULL;

	// South pole
	apush(verts, ((MeshVertex){ V3(0, -half_height - radius, 0), V3(0, -1, 0) }));

	// Bottom hemisphere rings (south pole toward equator)
	for (int r = 0; r < rings; r++) {
		float lat = -PI/2 + PI/2 * (float)(r + 1) / (rings + 1);
		float y = -half_height + sinf(lat) * radius;
		float rr = cosf(lat) * radius;
		for (int s = 0; s < slices; s++) {
			float lon = 2 * PI * s / slices;
			float cx = cosf(lon), sz = sinf(lon);
			v3 pos = V3(rr * cx, y, rr * sz);
			v3 nrm = V3(cosf(lat)*cx, sinf(lat), cosf(lat)*sz);
			apush(verts, ((MeshVertex){ pos, nrm }));
		}
	}

	// Cylinder bottom ring
	for (int s = 0; s < slices; s++) {
		float lon = 2 * PI * s / slices;
		float cx = cosf(lon), sz = sinf(lon);
		apush(verts, ((MeshVertex){ V3(radius*cx, -half_height, radius*sz), V3(cx, 0, sz) }));
	}

	// Cylinder top ring
	for (int s = 0; s < slices; s++) {
		float lon = 2 * PI * s / slices;
		float cx = cosf(lon), sz = sinf(lon);
		apush(verts, ((MeshVertex){ V3(radius*cx, half_height, radius*sz), V3(cx, 0, sz) }));
	}

	// Top hemisphere rings (equator toward north pole)
	for (int r = 0; r < rings; r++) {
		float lat = PI/2 * (float)(r + 1) / (rings + 1);
		float y = half_height + sinf(lat) * radius;
		float rr = cosf(lat) * radius;
		for (int s = 0; s < slices; s++) {
			float lon = 2 * PI * s / slices;
			float cx = cosf(lon), sz = sinf(lon);
			v3 pos = V3(rr * cx, y, rr * sz);
			v3 nrm = V3(cosf(lat)*cx, sinf(lat), cosf(lat)*sz);
			apush(verts, ((MeshVertex){ pos, nrm }));
		}
	}

	// North pole
	apush(verts, ((MeshVertex){ V3(0, half_height + radius, 0), V3(0, 1, 0) }));

	// Index layout
	int south = 0;
	int bh = 1;                          // bottom hemisphere start
	int cb = bh + rings * slices;        // cylinder bottom ring
	int ct = cb + slices;                // cylinder top ring
	int th = ct + slices;                // top hemisphere start
	int north = th + rings * slices;     // north pole

	// Helper macro for ring-strip triangles between two rings
	#define RING_STRIP(r0, r1) \
		for (int s = 0; s < slices; s++) { \
			int sn = (s + 1) % slices; \
			apush(tris, (uint16_t)(r0 + s)); \
			apush(tris, (uint16_t)(r0 + sn)); \
			apush(tris, (uint16_t)(r1 + sn)); \
			apush(tris, (uint16_t)(r0 + s)); \
			apush(tris, (uint16_t)(r1 + sn)); \
			apush(tris, (uint16_t)(r1 + s)); \
		}

	// South pole fan
	for (int s = 0; s < slices; s++) {
		int sn = (s + 1) % slices;
		apush(tris, (uint16_t)south);
		apush(tris, (uint16_t)(bh + sn));
		apush(tris, (uint16_t)(bh + s));
	}

	// Bottom hemisphere ring strips
	for (int r = 0; r < rings - 1; r++)
		RING_STRIP(bh + r * slices, bh + (r + 1) * slices)

	// Bottom hemisphere to cylinder bottom
	RING_STRIP(bh + (rings - 1) * slices, cb)

	// Cylinder
	RING_STRIP(cb, ct)

	// Cylinder top to top hemisphere
	RING_STRIP(ct, th)

	// Top hemisphere ring strips
	for (int r = 0; r < rings - 1; r++)
		RING_STRIP(th + r * slices, th + (r + 1) * slices)

	// North pole fan
	for (int s = 0; s < slices; s++) {
		int sn = (s + 1) % slices;
		int lr = th + (rings - 1) * slices;
		apush(tris, (uint16_t)(lr + s));
		apush(tris, (uint16_t)(lr + sn));
		apush(tris, (uint16_t)north);
	}

	#undef RING_STRIP

	mesh_upload(mesh, verts, asize(verts), tris, asize(tris));
	afree(verts);
	afree(tris);
}

// -----------------------------------------------------------------------------
// Hull mesh: triangulate faces as fans, bake scale into vertices.

static void mesh_generate_hull(Mesh* mesh, const Hull* hull, v3 sc)
{
	MeshVertex* verts = NULL;
	uint16_t* tris = NULL;

	for (int f = 0; f < hull->face_count; f++) {
		v3 n = hull->planes[f].normal;
		v3 sn = v3_norm(V3(n.x / sc.x, n.y / sc.y, n.z / sc.z));

		int fan_start = asize(verts);
		int fan_count = 0;
		int start = hull->faces[f].edge;
		int e = start;
		do {
			v3 v = hull->verts[hull->edges[e].origin];
			apush(verts, ((MeshVertex){ V3(v.x*sc.x, v.y*sc.y, v.z*sc.z), sn }));
			fan_count++;
			e = hull->edges[e].next;
		} while (e != start);

		for (int i = 1; i < fan_count - 1; i++) {
			apush(tris, (uint16_t)fan_start);
			apush(tris, (uint16_t)(fan_start + i));
			apush(tris, (uint16_t)(fan_start + i + 1));
		}
	}

	mesh_upload(mesh, verts, asize(verts), tris, asize(tris));
	afree(verts);
	afree(tris);
}

// -----------------------------------------------------------------------------
// Render state.

static GLuint r_program;
static GLint r_loc_vp;
static GLint r_loc_light_dir;
static GLint r_loc_ambient;
static Mesh r_meshes[MAX_MESH_TYPES];
static int r_mesh_ready[MAX_MESH_TYPES];
static RenderInstance* r_instances[MAX_MESH_TYPES]; // dynamic arrays, one per mesh type
static int r_mesh_count = MESH_BUILTIN_COUNT;
static mat4 r_vp;
static v3 r_light_dir;
static v3 r_ambient;

static Mesh* get_mesh(MeshType type)
{
	if (!r_mesh_ready[type]) {
		switch (type) {
		case MESH_BOX:    mesh_generate_box(&r_meshes[type]); break;
		case MESH_SPHERE: mesh_generate_sphere(&r_meshes[type]); break;
		default: break;
		}
		r_mesh_ready[type] = 1;
	}
	return &r_meshes[type];
}

// -----------------------------------------------------------------------------
// Custom mesh registration.

int render_create_capsule_mesh(float radius, float half_height)
{
	assert(r_mesh_count < MAX_MESH_TYPES);
	int idx = r_mesh_count++;
	mesh_generate_capsule(&r_meshes[idx], radius, half_height);
	r_mesh_ready[idx] = 1;
	return idx;
}

int render_create_hull_mesh(const Hull* hull, v3 sc)
{
	assert(r_mesh_count < MAX_MESH_TYPES);
	int idx = r_mesh_count++;
	mesh_generate_hull(&r_meshes[idx], hull, sc);
	r_mesh_ready[idx] = 1;
	return idx;
}

void render_init()
{
	gl_load();
	r_program = create_program(s_vert_src, s_frag_src);
	r_loc_vp = gl_GetUniformLocation(r_program, "u_vp");
	r_loc_light_dir = gl_GetUniformLocation(r_program, "u_light_dir");
	r_loc_ambient = gl_GetUniformLocation(r_program, "u_ambient");
	r_light_dir = norm(V3(0.3f, 1.0f, 0.5f));
	r_ambient = V3(0.15f, 0.15f, 0.18f);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA); // premultiplied alpha
}

void render_set_light_dir(v3 dir)
{
	r_light_dir = norm(dir);
}

void render_set_ambient(v3 color)
{
	r_ambient = color;
}

void render_begin(mat4 vp)
{
	r_vp = vp;
	for (int i = 0; i < r_mesh_count; i++) aclear(r_instances[i]);
}

void render_push(int type, mat4 model, v3 color, float opacity)
{
	RenderInstance inst;
	inst.model = model;
	inst.color[0] = color.x * opacity; // premultiply
	inst.color[1] = color.y * opacity;
	inst.color[2] = color.z * opacity;
	inst.color[3] = opacity;
	apush(r_instances[type], inst);
}

void render_end()
{
	gl_UseProgram(r_program);
	gl_UniformMatrix4fv(r_loc_vp, 1, GL_FALSE, r_vp.m);
	gl_Uniform3f(r_loc_light_dir, r_light_dir.x, r_light_dir.y, r_light_dir.z);
	gl_Uniform3f(r_loc_ambient, r_ambient.x, r_ambient.y, r_ambient.z);

	for (int i = 0; i < r_mesh_count; i++) {
		int count = asize(r_instances[i]);
		if (count == 0) continue;

		Mesh* m = get_mesh(i);
		gl_BindVertexArray(m->vao);

		// Upload instance data
		gl_BindBuffer(GL_ARRAY_BUFFER, m->instance_vbo);
		gl_BufferSubData(GL_ARRAY_BUFFER, 0, count * sizeof(RenderInstance), r_instances[i]);

		gl_DrawElementsInstanced(GL_TRIANGLES, m->index_count, GL_UNSIGNED_SHORT, NULL, count);
	}

	gl_BindVertexArray(0);
}
