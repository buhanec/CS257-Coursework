/**
 * The function to optimise as part of the coursework.
 *
 * l0, l1, l2 and l3 record the amount of time spent in each loop
 * and should not be optimised out. :)
 * strip mining?
 */

#define BLOCK_SIZE 4

void compute() {
    void call(float v) {
    }
    /*
    // External loop timers.
    double l0, l1, l2, l3;

    // Simulation parameters.
    int N;
    int steps;
    float eps = 0.00125f;
    float dmp = 0.995f;
    float dt = 0.001f;
    int vis = 0;

    // Solution arrays.
    float *x, *y, *z;
    float *ax, *ay, *az;
    float *vx, *vy, *vz;
    float *m;
    float *c;
    //*/

    // Preponderation
    float factor = dmp * dt;
    int V = (N/BLOCK_SIZE)*BLOCK_SIZE;
    int B = N/BLOCK_SIZE;
    int i, j, k, l;

    // Loop timers
    double t0, t1;


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
    for (i = V; i < N; ++i) {
        ax[i] = 0.0f;
        ay[i] = 0.0f;
        az[i] = 0.0f;
    }
    t1 = wtime();
    l0 += (t1 - t0);


    // Loop 1. http://www.hackersdelight.org/hdcodetxt/rsqrt.c.txt
    t0 = wtime();
    if (1) {
    for (k = 0; k < V; k += 4) {
        for (l = 0; l < V; l += 4) {
            for (i = k; i < k+4; i++) {
                for (j = l; j < l+4; j++) {
                    float rx = x[j] - x[i];
                    float ry = y[j] - y[i];
                    float rz = z[j] - z[i];
                    float r2 = rx*rx + ry*ry + rz*rz + eps;
                    //call(r2);
                    //float r2inv = 1.0f / sqrt(r2);
                    //float r6inv = r2inv * r2inv * r2inv;
                    //float s = m[j] * r6inv;
                    float s = 1.0f * r2;
                    ax[i] += s * rx;
                    ay[i] += s * ry;
                    az[i] += s * rz;
                }
            }
        }
    }
    } else {
    for (int i = 0; i < N; i++) {
        //float sum_x = 0;
        //float sum_y = 0;
        //float sum_z = 0;
        for (int j = 0; j < N; j++) {
            float rx = x[j] - x[i];
            float ry = y[j] - y[i];
            float rz = z[j] - z[i];
            float r2 = rx*rx + ry*ry + rz*rz + eps;
            //float r2inv = 1.0f / sqrt(r2);
            //float r6inv = r2inv * r2inv * r2inv;
            //float s = m[j] * r6inv;
            float s = 1.0f * r2;
            ax[i] += s * rx;
            //vx[i] += factor * s * rx;
            ay[i] += s * ry;
            //vy[i] += factor * s * rx;
            az[i] += s * rz;
            //vz[i] += factor * s * rz;
        }
        //ax[i] += sum_x;
        //ay[i] += sum_y;
        //az[i] += sum_z;
    }
    }
    t1 = wtime();
    l1 += (t1 - t0);

    // Loop 2.
    t0 = wtime();
    for (int i = 0; i < N; i++) {
        vx[i] += factor * ax[i];
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
