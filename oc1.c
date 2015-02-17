/**
 * The function to optimise.
 *
 * l0, l1, l2 and l3 record the amount of time spent in each loop
 * and should not be optimised out. :)
 */
#include <immintrin.h>
#include <omp.h>
#include <memory.h>


void compute() {

    omp_set_num_threads(omp_get_num_threads());
    int n = omp_get_num_threads();
    //Fast inverse square, thank you based Carmack
    float InvSqrt(float x){
        float xhalf = 0.5f * x;
        int i = *(int*)&x;
        i = 0x5f3759d5 - (i >> 1);
        x = *(float*)&i;
        x = x*(1.5f - xhalf*x*x);
        return x;
    }

    double t0, t1;

    // Loop 0.
    t0 = wtime();
    memset(ax,0, sizeof(float)*N);
    memset(ay,0, sizeof(float)*N);
    memset(az,0, sizeof(float)*N);
    t1 = wtime();
    l0 += (t1 - t0);

    // Loop 1.
    t0 = wtime();
    int i, j;
    int i0, j0;
    float rx, ry, rz, r2, r2inv, r6inv, s;


    #pragma omp parallel private (rx, ry, rz, r2, r2inv, r6inv, s)

    {
    for ( j = 0;j < N-2; j+=3) {

    #pragma omp for

        for (i = 0; i < N; i++) {

            rx = x[j] - x[i];
            ry = y[j] - y[i];
            rz = z[j] - z[i];
            r2 = rx*rx + ry*ry + rz*rz + eps;
            r2inv = InvSqrt(r2);
            r6inv = r2inv * r2inv *r2inv;
            s = m[j] * r6inv;
            ax[i] += s * rx;
            ay[i] += s * ry;
            az[i] += s * rz;

            rx = x[j+1] - x[i+1];
            ry = y[j+1] - y[i+1];
            rz = z[j+1] - z[i+1];
            r2 = rx*rx + ry*ry + rz*rz + eps;
            r2inv = InvSqrt(r2);
            r6inv = r2inv * r2inv *r2inv;
            s = m[j+1] * r6inv;
            ax[i+1] += s * rx;
            ay[i+1] += s * ry;
            az[i+1] += s * rz;

            rx = x[j+2] - x[i+2];
            ry = y[j+2] - y[i+2];
            rz = z[j+2] - z[i+2];
            r2 = rx*rx + ry*ry + rz*rz + eps;
            r2inv = InvSqrt(r2);
            r6inv = r2inv * r2inv *r2inv;
            s = m[j+2] * r6inv;
            ax[i+2] += s * rx;
            ay[i+2] += s * ry;
            az[i+2] += s * rz;

         }


}
}
    t1 = wtime();
    l1 += (t1 - t0);

    // Loop 2.
    t0 = wtime();
    float dd = dmp * dt;
    for (int i = 0; i < N; i++) {
        vx[i] += dd * ax[i];
        vy[i] += dd * ay[i];
        vz[i] += dd * az[i];
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
