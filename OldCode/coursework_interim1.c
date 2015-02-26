/**
 * The function to optimise as part of the coursework.
 *
 * l0, l1, l2 and l3 record the amount of time spent in each loop
 * and should not be optimised out. :)
 */
void computec() {
	int min(int x, int y) {
		return (x < y) ? x : y;
	}

	double t0, t1;
	float factor = dmp * dt;
	int b = 4;
	int V = (N/b)*b;
	int i, j, k, l;

	// Loop 0.
	// Loop 0.
	t0 = wtime();
	for (i = 0; i < V; i += 4) {
		ax[i] = 0.0f;
		ax[i+1] = 0.0f;
		ax[i+2] = 0.0f;
		ax[i+3] = 0.0f;
	}
	for (i = 0; i < V; i += 4) {
		ay[i] = 0.0f;
		ay[i+1] = 0.0f;
		ay[i+2] = 0.0f;
		ay[i+3] = 0.0f;
	}
	for (i = 0; i < V; i += 4) {
		az[i] = 0.0f;
		az[i+1] = 0.0f;
		az[i+2] = 0.0f;
		az[i+3] = 0.0f;
	}
	for (; i < N; i++) {
		ax[i] = 0.0f;
		ay[i] = 0.0f;
		az[i] = 0.0f;
	}
	t1 = wtime();
	l0 += (t1 - t0);

    // Loop 1.
	t0 = wtime();
	{
		for (i = 0; i < N; i += b) {
			//Replaced by array a
			//float tx = 0.0f;
			//float ty = 0.0f;
			//float tz = 0.0f;
			for (j = 0; j < N; j += b) {
				for (k = i; k < min(i+b, N); k++) {
					for (l = j; l < min(j+b, N); l++) {
						float rx = x[l] - x[k];
						float ry = y[l] - y[k];
						float rz = z[l] - z[k];
						float r2 = rx*rx + ry*ry + rz*rz + eps;
						//float r2inv = 1.0f / sqrt(r2);
						//float r6inv = r2inv * r2inv * r2inv;
						float r6inv = 1.0f / (r2 * sqrt(r2));
						float s = m[l] * r6inv;
						ax[i] += s * rx;
						ay[i] += s * ry;
						az[i] += s * rz;
					}
				}
			}
		}
	}
	t1 = wtime();

	l1 += (t1 - t0);

	// Loop 2.
	t0 = wtime();
	for (i = 0; i < V; i += 4) {
		vx[i] += factor * ax[i];
		vx[i+1] = factor * ax[i+1];
		vx[i+2] = factor * ax[i+2];
		vx[i+3] = factor * ax[i+3];
	}
	for (i = 0; i < V; i += 4) {
		vy[i] += factor * ay[i];
		vy[i+1] = factor * ay[i+1];
		vy[i+2] = factor * ay[i+2];
		vy[i+3] = factor * ay[i+3];
	}
	for (i = 0; i < V; i += 4) {
		vz[i] += factor * az[i];
		vz[i+1] = factor * az[i+1];
		vz[i+2] = factor * az[i+2];
		vz[i+3] = factor * az[i+3];
	}
	for (; i < N; i++) {
		vx[i] += factor * az[i];
		vy[i] += factor * ay[i];
		vz[i] += factor * az[i];
	}
	t1 = wtime();
	l2 += (t1 - t0);

	// Loop 3.
	t0 = wtime();
	for (int i = 0; i < N; i++) {
        x[i] += dt * vx[i];
        y[i] += dt * vy[i];
        z[i] += dt * vz[i];
		if (x[i] >= 1.0f || x[i] <= -1.0f) vx[i] *= -1.0f;
		if (y[i] >= 1.0f || y[i] <= -1.0f) vy[i] *= -1.0f;
		if (z[i] >= 1.0f || z[i] <= -1.0f) vz[i] *= -1.0f;

	}
	t1 = wtime();
	l3 += (t1 - t0);

}
