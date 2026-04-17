// Debug inspection of EPA on canonical inputs.
static void debug_epa()
{
	const Hull* box = hull_unit_box();
	quat id = quat_identity();

	printf("\n=== EPA DEBUG ===\n");
	{
		printf("-- sphere(r=1) at origin vs unit box at origin --\n");
		Manifold m = {0};
		ConvexHull hb = { box, V3(0,0,0), id, V3(1,1,1) };
		int hit = epa_sphere_hull((Sphere){ V3(0,0,0), 1.0f }, hb, &m);
		printf("hit=%d count=%d\n", hit, m.count);
		if (hit && m.count > 0) {
			printf("  pt=(%.3f,%.3f,%.3f) n=(%.3f,%.3f,%.3f) depth=%.4f fid=%u\n",
				m.contacts[0].point.x, m.contacts[0].point.y, m.contacts[0].point.z,
				m.contacts[0].normal.x, m.contacts[0].normal.y, m.contacts[0].normal.z,
				m.contacts[0].penetration, m.contacts[0].feature_id);
		}
	}
	{
		printf("-- box-box overlap +X (at 1.0) --\n");
		Manifold m = {0};
		ConvexHull ha = { box, V3(0,0,0), id, V3(1,1,1) };
		ConvexHull hb = { box, V3(1.0f,0,0), id, V3(1,1,1) };
		int hit = epa_hull_hull(ha, hb, &m);
		printf("hit=%d count=%d\n", hit, m.count);
		if (hit && m.count > 0) {
			printf("  pt=(%.3f,%.3f,%.3f) n=(%.3f,%.3f,%.3f) depth=%.4f fid=%u\n",
				m.contacts[0].point.x, m.contacts[0].point.y, m.contacts[0].point.z,
				m.contacts[0].normal.x, m.contacts[0].normal.y, m.contacts[0].normal.z,
				m.contacts[0].penetration, m.contacts[0].feature_id);
		}
	}
	{
		printf("-- box-box overlap +Y (at 1.0) --\n");
		Manifold m = {0};
		ConvexHull ha = { box, V3(0,0,0), id, V3(1,1,1) };
		ConvexHull hb = { box, V3(0,1.0f,0), id, V3(1,1,1) };
		int hit = epa_hull_hull(ha, hb, &m);
		printf("hit=%d count=%d\n", hit, m.count);
		if (hit && m.count > 0) {
			printf("  pt=(%.3f,%.3f,%.3f) n=(%.3f,%.3f,%.3f) depth=%.4f fid=%u\n",
				m.contacts[0].point.x, m.contacts[0].point.y, m.contacts[0].point.z,
				m.contacts[0].normal.x, m.contacts[0].normal.y, m.contacts[0].normal.z,
				m.contacts[0].penetration, m.contacts[0].feature_id);
		}
	}
	{
		printf("-- sphere-box shallow +X (1.5, 0, 0) r=1 --\n");
		Manifold m = {0};
		ConvexHull hb = { box, V3(0,0,0), id, V3(1,1,1) };
		int hit = epa_sphere_hull((Sphere){ V3(1.5f,0,0), 1.0f }, hb, &m);
		printf("hit=%d count=%d\n", hit, m.count);
		if (hit && m.count > 0) {
			printf("  pt=(%.3f,%.3f,%.3f) n=(%.3f,%.3f,%.3f) depth=%.4f fid=%u\n",
				m.contacts[0].point.x, m.contacts[0].point.y, m.contacts[0].point.z,
				m.contacts[0].normal.x, m.contacts[0].normal.y, m.contacts[0].normal.z,
				m.contacts[0].penetration, m.contacts[0].feature_id);
		}
	}
	printf("=== END EPA DEBUG ===\n\n");
}
