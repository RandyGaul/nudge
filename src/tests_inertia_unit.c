// tests_inertia_unit.c -- unit tests for inv_inertia_mul and dinv_inertia_mul.
// These compute R * diag(inv_I) * R^T * v: the world-space inverse inertia
// tensor applied to a vector. Used everywhere in the solver.

#define INERTIA_EPS 1e-6

// Helper: compute R * diag(inv_I) * R^T * v using explicit 3x3 matrix math.
// This is the reference implementation -- no quaternion tricks, just build
// the rotation matrix, multiply it out, compare against the function under test.
static v3 ref_inv_inertia_mul(quat rot, v3 inv_i, v3 v)
{
	// Build rotation matrix columns from quaternion
	float x = rot.x, y = rot.y, z = rot.z, w = rot.w;
	// Column 0 of R
	float r00 = 1 - 2*(y*y + z*z), r10 = 2*(x*y + w*z), r20 = 2*(x*z - w*y);
	// Column 1 of R
	float r01 = 2*(x*y - w*z), r11 = 1 - 2*(x*x + z*z), r21 = 2*(y*z + w*x);
	// Column 2 of R
	float r02 = 2*(x*z + w*y), r12 = 2*(y*z - w*x), r22 = 1 - 2*(x*x + y*y);

	// R^T * v (rotate into body frame)
	float lx = r00*v.x + r10*v.y + r20*v.z;
	float ly = r01*v.x + r11*v.y + r21*v.z;
	float lz = r02*v.x + r12*v.y + r22*v.z;

	// Scale by diagonal inverse inertia
	lx *= inv_i.x; ly *= inv_i.y; lz *= inv_i.z;

	// R * scaled (rotate back to world frame)
	return V3(
		r00*lx + r01*ly + r02*lz,
		r10*lx + r11*ly + r12*lz,
		r20*lx + r21*ly + r22*lz
	);
}

// ============================================================================
// inv_inertia_mul (float) tests

static void test_inv_inertia_identity_rotation()
{
	TEST_BEGIN("inv_inertia_identity_rotation");
	// R = I: result is just diag(inv_I) * v.
	quat q = quat_identity();
	v3 inv_i = V3(2, 3, 5);
	v3 v = V3(1, 1, 1);
	v3 r = inv_inertia_mul(q, inv_i, v);
	TEST_ASSERT_FLOAT(r.x, 2.0f, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.y, 3.0f, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.z, 5.0f, INERTIA_EPS);
}

static void test_inv_inertia_uniform()
{
	TEST_BEGIN("inv_inertia_uniform");
	// Uniform inertia: R * k*I * R^T = k*I for any rotation.
	// Result should be k * v regardless of rotation.
	// Float32 rotate->scale->rotate accumulates ~1e-5 error on large values.
	quat q = quat_axis_angle(norm(V3(1, 1, 1)), 1.23f);
	v3 inv_i = V3(5, 5, 5);
	v3 v = V3(3, -2, 7);
	v3 r = inv_inertia_mul(q, inv_i, v);
	TEST_ASSERT_FLOAT(r.x, 15.0f, 1e-4f);
	TEST_ASSERT_FLOAT(r.y, -10.0f, 1e-4f);
	TEST_ASSERT_FLOAT(r.z, 35.0f, 1e-4f);
}

static void test_inv_inertia_90_about_z()
{
	TEST_BEGIN("inv_inertia_90_about_z");
	// 90 degrees about Z: R swaps X<->Y with sign flip.
	// R = [0 -1 0; 1 0 0; 0 0 1]
	// R * diag(a,b,c) * R^T: XX gets b, YY gets a, ZZ gets c, XY coupling.
	quat q = quat_axis_angle(V3(0, 0, 1), 3.14159265f / 2.0f);
	v3 inv_i = V3(2, 8, 5);
	v3 v = V3(1, 0, 0);
	v3 r = inv_inertia_mul(q, inv_i, v);
	v3 ref = ref_inv_inertia_mul(q, inv_i, v);
	TEST_ASSERT_FLOAT(r.x, ref.x, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.y, ref.y, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.z, ref.z, INERTIA_EPS);
	// After 90 about Z, applying along X should give inv_i.y = 8 in X direction.
	TEST_ASSERT_FLOAT(r.x, 8.0f, 0.01f);
}

static void test_inv_inertia_90_about_x()
{
	TEST_BEGIN("inv_inertia_90_about_x");
	// 90 degrees about X: R swaps Y<->Z.
	quat q = quat_axis_angle(V3(1, 0, 0), 3.14159265f / 2.0f);
	v3 inv_i = V3(2, 8, 5);
	v3 v = V3(0, 1, 0);
	v3 r = inv_inertia_mul(q, inv_i, v);
	v3 ref = ref_inv_inertia_mul(q, inv_i, v);
	TEST_ASSERT_FLOAT(r.x, ref.x, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.y, ref.y, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.z, ref.z, INERTIA_EPS);
	// After 90 about X, applying along Y should give inv_i.z = 5 in Y direction.
	TEST_ASSERT_FLOAT(r.y, 5.0f, 0.01f);
}

static void test_inv_inertia_arbitrary_rotation()
{
	TEST_BEGIN("inv_inertia_arbitrary_rotation");
	// Arbitrary rotation and inertia: compare against reference implementation.
	// Tolerance is looser than INERTIA_EPS because the quaternion-based inv_inertia_mul
	// and the explicit 3x3 matrix ref_inv_inertia_mul accumulate float32 error differently
	// (~1e-5 drift on values of magnitude ~10 after rotate->scale->rotate).
	quat q = quat_axis_angle(norm(V3(1, 2, 3)), 0.7f);
	v3 inv_i = V3(1.5f, 0.3f, 7.0f);
	v3 v = V3(-2, 5, 1);
	v3 r = inv_inertia_mul(q, inv_i, v);
	v3 ref = ref_inv_inertia_mul(q, inv_i, v);
	TEST_ASSERT_FLOAT(r.x, ref.x, 1e-5f);
	TEST_ASSERT_FLOAT(r.y, ref.y, 1e-5f);
	TEST_ASSERT_FLOAT(r.z, ref.z, 1e-5f);
}

static void test_inv_inertia_zero_vector()
{
	TEST_BEGIN("inv_inertia_zero_vector");
	// Zero input vector: result must be zero.
	quat q = quat_axis_angle(V3(0, 1, 0), 1.0f);
	v3 inv_i = V3(2, 3, 5);
	v3 v = V3(0, 0, 0);
	v3 r = inv_inertia_mul(q, inv_i, v);
	TEST_ASSERT_FLOAT(r.x, 0.0f, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.y, 0.0f, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.z, 0.0f, INERTIA_EPS);
}

static void test_inv_inertia_static_body()
{
	TEST_BEGIN("inv_inertia_static_body");
	// Static body: inv_inertia = (0,0,0). Result must be zero.
	quat q = quat_axis_angle(norm(V3(1, 1, 1)), 0.5f);
	v3 inv_i = V3(0, 0, 0);
	v3 v = V3(10, 20, 30);
	v3 r = inv_inertia_mul(q, inv_i, v);
	TEST_ASSERT_FLOAT(r.x, 0.0f, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.y, 0.0f, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.z, 0.0f, INERTIA_EPS);
}

static void test_inv_inertia_one_axis_locked()
{
	TEST_BEGIN("inv_inertia_one_axis_locked");
	// One inertia axis is zero (infinite mass in that direction).
	// E.g. a body constrained to only rotate about X in local frame.
	quat q = quat_identity();
	v3 inv_i = V3(5, 0, 0);
	v3 v = V3(1, 1, 1);
	v3 r = inv_inertia_mul(q, inv_i, v);
	TEST_ASSERT_FLOAT(r.x, 5.0f, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.y, 0.0f, INERTIA_EPS);
	TEST_ASSERT_FLOAT(r.z, 0.0f, INERTIA_EPS);
}

// ============================================================================
// dinv_inertia_mul (double) tests -- must agree with float version

static void test_dinv_inertia_matches_float()
{
	TEST_BEGIN("dinv_inertia_matches_float");
	// Arbitrary rotation: double version must agree with float version.
	// Both use float32 quaternion components, so intermediate precision differs
	// slightly. 1e-4 tolerance accounts for float32 accumulation.
	quat q = quat_axis_angle(norm(V3(3, -1, 2)), 1.1f);
	v3 inv_i = V3(1.5f, 0.3f, 7.0f);
	dv3 dv = DV3(-2, 5, 1);
	v3 fv = V3(-2, 5, 1);
	dv3 dr = dinv_inertia_mul(q, inv_i, dv);
	v3 fr = inv_inertia_mul(q, inv_i, fv);
	TEST_ASSERT_FLOAT((float)dr.x, fr.x, 1e-4f);
	TEST_ASSERT_FLOAT((float)dr.y, fr.y, 1e-4f);
	TEST_ASSERT_FLOAT((float)dr.z, fr.z, 1e-4f);
}

static void test_dinv_inertia_identity_rotation()
{
	TEST_BEGIN("dinv_inertia_identity_rotation");
	quat q = quat_identity();
	v3 inv_i = V3(2, 3, 5);
	dv3 v = DV3(1, 1, 1);
	dv3 r = dinv_inertia_mul(q, inv_i, v);
	TEST_ASSERT(fabs(r.x - 2.0) < 1e-12);
	TEST_ASSERT(fabs(r.y - 3.0) < 1e-12);
	TEST_ASSERT(fabs(r.z - 5.0) < 1e-12);
}

static void test_dinv_inertia_uniform()
{
	TEST_BEGIN("dinv_inertia_uniform");
	quat q = quat_axis_angle(norm(V3(1, 1, 1)), 1.23f);
	v3 inv_i = V3(5, 5, 5);
	dv3 v = DV3(3, -2, 7);
	dv3 r = dinv_inertia_mul(q, inv_i, v);
	TEST_ASSERT(fabs(r.x - 15.0) < 1e-4);
	TEST_ASSERT(fabs(r.y + 10.0) < 1e-4);
	TEST_ASSERT(fabs(r.z - 35.0) < 1e-4);
}

static void test_dinv_inertia_reference()
{
	TEST_BEGIN("dinv_inertia_reference");
	// Compare double version against explicit matrix reference.
	quat q = quat_axis_angle(norm(V3(1, 2, 3)), 0.7f);
	v3 inv_i = V3(1.5f, 0.3f, 7.0f);
	dv3 dv = DV3(-2, 5, 1);
	dv3 dr = dinv_inertia_mul(q, inv_i, dv);
	v3 ref = ref_inv_inertia_mul(q, inv_i, V3(-2, 5, 1));
	TEST_ASSERT(fabs(dr.x - ref.x) < 1e-4);
	TEST_ASSERT(fabs(dr.y - ref.y) < 1e-4);
	TEST_ASSERT(fabs(dr.z - ref.z) < 1e-4);
}

static void test_dinv_inertia_static_body()
{
	TEST_BEGIN("dinv_inertia_static_body");
	quat q = quat_axis_angle(norm(V3(1, 1, 1)), 0.5f);
	v3 inv_i = V3(0, 0, 0);
	dv3 v = DV3(10, 20, 30);
	dv3 r = dinv_inertia_mul(q, inv_i, v);
	TEST_ASSERT(fabs(r.x) < 1e-12);
	TEST_ASSERT(fabs(r.y) < 1e-12);
	TEST_ASSERT(fabs(r.z) < 1e-12);
}

// ============================================================================
// Runner

static void run_inertia_unit_tests()
{
	printf("--- inertia unit tests ---\n");

	// Float version
	test_inv_inertia_identity_rotation();
	test_inv_inertia_uniform();
	test_inv_inertia_90_about_z();
	test_inv_inertia_90_about_x();
	test_inv_inertia_arbitrary_rotation();
	test_inv_inertia_zero_vector();
	test_inv_inertia_static_body();
	test_inv_inertia_one_axis_locked();

	// Double version
	test_dinv_inertia_matches_float();
	test_dinv_inertia_identity_rotation();
	test_dinv_inertia_uniform();
	test_dinv_inertia_reference();
	test_dinv_inertia_static_body();
}
