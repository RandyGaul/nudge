// tests_misc_ldl_unit.c -- unit tests for apply_rotation_delta, ldl_validate_lambda,
// and ldl_condition_check.

// ============================================================================
// apply_rotation_delta tests
//
// Integrates q += 0.5 * (w, 0) * q, then normalizes.
// This is the first-order quaternion integration used for position correction.

static float quat_length(quat q)
{
	return sqrtf(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
}

static void test_rotdelta_zero()
{
	TEST_BEGIN("rotdelta_zero");
	// Zero angular velocity: quaternion unchanged.
	quat q = quat_identity();
	apply_rotation_delta(&q, DV3(0, 0, 0));
	TEST_ASSERT_FLOAT(q.x, 0.0f, 1e-7f);
	TEST_ASSERT_FLOAT(q.y, 0.0f, 1e-7f);
	TEST_ASSERT_FLOAT(q.z, 0.0f, 1e-7f);
	TEST_ASSERT_FLOAT(q.w, 1.0f, 1e-7f);
}

static void test_rotdelta_small_about_z()
{
	TEST_BEGIN("rotdelta_small_about_z");
	// Small rotation about Z. Starting from identity.
	// For small angle theta, the update is approximately a rotation of theta/2
	// (since the integration is q += 0.5 * spin * q).
	// Actually: spin = (0, 0, w, 0), dq = spin * q.
	// For identity q: dq = (0, 0, w, 0) * (0, 0, 0, 1) = (0, 0, w, 0).
	// q_new = q + 0.5*dq = (0, 0, 0.5*w, 1), normalized.
	// This approximates a rotation of angle w about Z.
	double w = 0.1;
	quat q = quat_identity();
	apply_rotation_delta(&q, DV3(0, 0, w));

	// Should be approximately quat_axis_angle(Z, w)
	quat expected = quat_axis_angle(V3(0, 0, 1), (float)w);
	TEST_ASSERT_FLOAT(q.x, expected.x, 0.01f);
	TEST_ASSERT_FLOAT(q.y, expected.y, 0.01f);
	TEST_ASSERT_FLOAT(q.z, expected.z, 0.01f);
	TEST_ASSERT_FLOAT(q.w, expected.w, 0.01f);
	TEST_ASSERT_FLOAT(quat_length(q), 1.0f, 1e-6f);
}

static void test_rotdelta_normalized()
{
	TEST_BEGIN("rotdelta_normalized");
	// After any rotation delta, result must be unit length.
	quat q = quat_axis_angle(norm(V3(1, 2, 3)), 0.5f);
	apply_rotation_delta(&q, DV3(1.0, -0.5, 0.3));
	TEST_ASSERT_FLOAT(quat_length(q), 1.0f, 1e-6f);
}

static void test_rotdelta_large_angle()
{
	TEST_BEGIN("rotdelta_large_angle");
	// Large angular delta (w = 10). First-order integration is inaccurate
	// for large angles, but the result must still be a valid unit quaternion.
	quat q = quat_identity();
	apply_rotation_delta(&q, DV3(10, 0, 0));
	TEST_ASSERT_FLOAT(quat_length(q), 1.0f, 1e-6f);
	// Should have rotated significantly about X
	TEST_ASSERT(fabsf(q.x) > 0.1f);
}

static void test_rotdelta_preserves_non_identity()
{
	TEST_BEGIN("rotdelta_preserves_non_identity");
	// Start from a non-identity rotation, apply a delta, verify still unit length
	// and that the rotation actually changed.
	quat q = quat_axis_angle(V3(0, 1, 0), 1.0f);
	quat q_before = q;
	apply_rotation_delta(&q, DV3(0, 0, 0.5));
	TEST_ASSERT_FLOAT(quat_length(q), 1.0f, 1e-6f);
	// Should have changed
	float diff = fabsf(q.x - q_before.x) + fabsf(q.y - q_before.y) + fabsf(q.z - q_before.z) + fabsf(q.w - q_before.w);
	TEST_ASSERT(diff > 0.01f);
}

static void test_rotdelta_opposite_directions()
{
	TEST_BEGIN("rotdelta_opposite_directions");
	// Apply +w then -w. Should approximately return to original (first-order,
	// so not exact, but close for small angles).
	quat q = quat_identity();
	apply_rotation_delta(&q, DV3(0.05, 0, 0));
	apply_rotation_delta(&q, DV3(-0.05, 0, 0));
	TEST_ASSERT_FLOAT(q.x, 0.0f, 0.01f);
	TEST_ASSERT_FLOAT(q.y, 0.0f, 0.01f);
	TEST_ASSERT_FLOAT(q.z, 0.0f, 0.01f);
	TEST_ASSERT_FLOAT(q.w, 1.0f, 0.01f);
}

static void test_rotdelta_all_axes()
{
	TEST_BEGIN("rotdelta_all_axes");
	// Simultaneous rotation about all three axes. Verify unit length.
	quat q = quat_identity();
	apply_rotation_delta(&q, DV3(0.3, -0.2, 0.4));
	TEST_ASSERT_FLOAT(quat_length(q), 1.0f, 1e-6f);
	// All components should be nonzero
	TEST_ASSERT(fabsf(q.x) > 0.01f);
	TEST_ASSERT(fabsf(q.y) > 0.01f);
	TEST_ASSERT(fabsf(q.z) > 0.01f);
}

// ============================================================================
// ldl_validate_lambda tests

static void test_validate_lambda_nan()
{
	TEST_BEGIN("validate_lambda_nan");
	double zero = 0.0;
	double lam[3] = { 1.0, zero / zero, 2.0 };
	TEST_ASSERT(ldl_validate_lambda(lam, 3) == 0);
}

static void test_validate_lambda_inf()
{
	TEST_BEGIN("validate_lambda_inf");
	double zero = 0.0;
	double lam[3] = { 1.0, 1.0 / zero, 2.0 };
	TEST_ASSERT(ldl_validate_lambda(lam, 3) == 0);
}

// ============================================================================
// ldl_condition_check tests

// ============================================================================
// Runner

static void run_misc_ldl_unit_tests()
{
	printf("--- misc LDL unit tests ---\n");

	// Rotation delta
	test_rotdelta_zero();
	test_rotdelta_small_about_z();
	test_rotdelta_normalized();
	test_rotdelta_large_angle();
	test_rotdelta_preserves_non_identity();
	test_rotdelta_opposite_directions();
	test_rotdelta_all_axes();

	// Validate lambda
	test_validate_lambda_nan();
	test_validate_lambda_inf();

	// Condition check
}
