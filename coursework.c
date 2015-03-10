/**
 * The function to optimise as part of the coursework.
 *
 * l0, l1, l2 and l3 record the amount of time spent in each loop
 * and should not be optimised out. :)
 */
void compute() {

    double t0, t1;

    // Loop 0.
    t0 = wtime();
    #pragma unroll (4)
    for (int i = 0; i < N; i++) {
        ax[i] = 0.0f;
    }
    #pragma unroll (4)
    for (int i = 0; i < N; i++) {
        ay[i] = 0.0f;
    }
    #pragma unroll (4)
    for (int i = 0; i < N; i++) {
        az[i] = 0.0f;
    }
    t1 = wtime();
    l0 += (t1 - t0);

    // Loop 1.
    t0 = wtime();
    for (int i = 0; i < N; i++) {
        float xi = x[i];
        float yi = y[i];
        float zi = z[i];
        float sx = 0.0f;
        float sy = 0.0f;
        float sz = 0.0f;
        for (int j = 0; j < N; j++) {
            float rx = x[j] - xi;
            float ry = y[j] - yi;
            float rz = z[j] - zi;
            float r2 = rx*rx + ry*ry + rz*rz + eps;
            float r2inv = 1.0f / sqrt(r2);
            float r6inv = r2inv * r2inv * r2inv;
            float s = m[j] * r6inv;
            sx += s * rx;
            sy += s * ry;
            sz += s * rz;
        }
        ax[i] = sx;
        ay[i] = sy;
        az[i] = sz;
    }
    t1 = wtime();
    l1 += (t1 - t0);

    // Loop 2.
    t0 = wtime();
    #pragma unroll (4)
    for (int i = 0; i < N; i++) {
        vx[i] += dmp * (dt * ax[i]);
    }
    #pragma unroll (4)
    for (int i = 0; i < N; i++) {
        vy[i] += dmp * (dt * ay[i]);
    }
    #pragma unroll (4)
    for (int i = 0; i < N; i++) {
        vz[i] += dmp * (dt * az[i]);
    }
    t1 = wtime();
    l2 += (t1 - t0);

    // Loop 3.
    t0 = wtime();
    #pragma unroll (4)
    for (int i = 0; i < N; i++) {
        x[i] += dt * vx[i];
        if (x[i] >= 1.0f || x[i] <= -1.0f) vx[i] *= -1.0f;
    }
    #pragma unroll (4)
    for (int i = 0; i < N; i++) {
        y[i] += dt * vy[i];
        if (y[i] >= 1.0f || y[i] <= -1.0f) vy[i] *= -1.0f;
    }
    #pragma unroll (4)
    for (int i = 0; i < N; i++) {
        z[i] += dt * vz[i];
        if (z[i] >= 1.0f || z[i] <= -1.0f) vz[i] *= -1.0f;
    }
    t1 = wtime();
    l3 += (t1 - t0);

}
