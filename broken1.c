#include <immintrin.h>

__m128 _mm_hsum_ps(__m128 a) {
    __m128 b = _mm_add_ps(a, _mm_movehl_ps(a, a));
    return _mm_add_ps(b, _mm_shuffle_ps(b, b, 1));
}

void compute() {
    // Preponderation
    float factor = dmp * dt;
    int V = (N/4)*4;
    int B = N/4;
    int i, j, k, l;

    // packaged floats
    __m128 dmp_ = _mm_load1_ps(&dmp);
    __m128 dt_ = _mm_load1_ps(&dt);
    __m128 eps_ = _mm_load1_ps(&eps);
    __m128 factor_ = _mm_mul_ps(dmp_, dt_);
    __m128 zero_ = _mm_setzero_ps();
    __m128 negzero_ = _mm_set1_ps(-0.0f);
    __m128 one_ = _mm_set1_ps(1.0f);
    __m128 negone_ = _mm_set1_ps(-1.0f);
    __m128 half_ = _mm_set1_ps(0.5f);
    __m128 three_ = _mm_set1_ps(3.0f);
    __m128 fsign_ = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));


    double t0, t1;

    // Loop 0.
    t0 = wtime();
    for (i = 0; i < V; i += 4) {
        _mm_store_ps(ax+i, zero_);
    }
    for (i = 0; i < V; i += 4) {
        _mm_store_ps(ay+i, zero_);
    }
    for (i = 0; i < V; i += 4) {
        _mm_store_ps(az+i, zero_);
    }
    for (; i < N; i++) {
        _mm_store_ss(ax+i, zero_);
        _mm_store_ss(ay+i, zero_);
        _mm_store_ss(az+i, zero_);
    }
    t1 = wtime();
    l0 += (t1 - t0);

    // Loop 1.
    t0 = wtime();
    for (i = 0; i < N; i++) {
        __m128 x1_ = _mm_load1_ps(x+i);
        __m128 y1_ = _mm_load1_ps(y+i);
        __m128 z1_ = _mm_load1_ps(z+i);
        for (j = 0; j < V; j += 4) {
            float rx = x[j] - x[i];
            float rx1 = x[j+1] - x[i];
            float rx2 = x[j+2] - x[i];
            float rx3 = x[j+3] - x[i];
            //__m128 x_ = _mm_load_ps(x+j);
            //__m128 rx_ = _mm_sub_ps(x_, x1_);
            float ry = y[j] - y[i];
            float ry1 = y[j+1] - y[i];
            float ry2 = y[j+2] - y[i];
            float ry3 = y[j+3] - y[i];
            //__m128 y_ = _mm_load_ps(y+j);
            //__m128 ry_ = _mm_sub_ps(y_, y1_);
            float rz = z[j] - z[i];
            float rz1 = z[j+1] - z[i];
            float rz2 = z[j+2] - z[i];
            float rz3 = z[j+3] - z[i];
            //__m128 z_ = _mm_load_ps(z+j);
            //__m128 rz_ = _mm_sub_ps(z_, z1_);
            float r2 = rx*rx + ry*ry + rz*rz + eps;
            float r21 = rx1*rx1 + ry1*ry1 + rz1*rz1 + eps;
            float r22 = rx2*rx2 + ry2*ry2 + rz2*rz2 + eps;
            float r23 = rx3*rx3 + ry3*ry3 + rz3*rz3 + eps;
            //__m128 x2_ = _mm_mul_ps(rx_, rx_);
            //__m128 y2_ = _mm_mul_ps(ry_, ry_);
            //__m128 z2_ = _mm_mul_ps(rz_, rz_);
            //__m128 r2_ = _mm_add_ps(_mm_add_ps(_mm_add_ps(x2_, y2_), z2_), eps_);
            float r2inv = 1.0f / sqrt(r2);
            float r2inv1 = 1.0f / sqrt(r21);
            float r2inv2 = 1.0f / sqrt(r22);
            float r2inv3 = 1.0f / sqrt(r23);
            //__m128 r2inv_ = _mm_div_ps(one_, _mm_sqrt_ps(r2_));
            float r6inv = r2inv * r2inv * r2inv;
            float r6inv1 = r2inv1 * r2inv1 * r2inv1;
            float r6inv2 = r2inv2 * r2inv2 * r2inv2;
            float r6inv3 = r2inv3 * r2inv3 * r2inv3;
            //__m128 r6inv_ = _mm_mul_ps(_mm_mul_ps(r2inv_, r2inv_), r2inv_);
            float s = m[j] * r6inv;
            float s1 = m[j+1] * r6inv1;
            float s2 = m[j+2] * r6inv2;
            float s3 = m[j+3] * r6inv3;
            //__m128 m_ = _mm_load_ps(m+j);
            //__m128 s_ = _mm_mul_ps(m_, r6inv_);
            ax[i] += s * rx;
            ax[i] += s1 * rx1;
            ax[i] += s2 * rx2;
            ax[i] += s3 * rx3;
            //__m128 mx_ = _mm_mul_ps(s_, rx_);
            //__m128 sx1_ = _mm_add_ps(mx_, _mm_movehl_ps(mx_, mx_));
            //__m128 sx2_ = _mm_add_ps(sx1_, _mm_shuffle_ps(sx1_, sx1_, 1));
            //__m128 ax_ = _mm_load1_ps(ax+i);
            //_mm_store_ss(ax+i, _mm_add_ss(ax_, sx2_));
            ay[i] += s * ry;
            ay[i] += s1 * ry1;
            ay[i] += s2 * ry2;
            ay[i] += s3 * ry3;
            //__m128 my_ = _mm_mul_ps(s_, ry_);
            //__m128 sy1_ = _mm_add_ps(my_, _mm_movehl_ps(my_, my_));
            //__m128 sy2_ = _mm_add_ps(sy1_, _mm_shuffle_ps(sy1_, sy1_, 1));
            //__m128 ay_ = _mm_load1_ps(ay+i);
            //_mm_store_ss(ay+i, _mm_add_ss(ay_, sy2_));
            az[i] += s * rz;
            az[i] += s1 * rz1;
            az[i] += s2 * rz2;
            az[i] += s3 * rz3;
            //__m128 mz_ = _mm_mul_ps(s_, rz_);
            //__m128 sz1_ = _mm_add_ps(mz_, _mm_movehl_ps(mz_, mz_));
            //__m128 sz2_ = _mm_add_ps(sz1_, _mm_shuffle_ps(sz1_, sz1_, 1));
            //__m128 az_ = _mm_load1_ps(az+i);
            //_mm_store_ss(az+i, _mm_add_ss(az_, sz2_));
        }
        for (; j < N; j++) {
            float rx = x[j] - x[i];
            float ry = y[j] - y[i];
            float rz = z[j] - z[i];
            float r2 = rx*rx + ry*ry + rz*rz + eps;
            float r2inv = 1.0f / sqrt(r2);
            float r6inv = r2inv * r2inv * r2inv;
            float s = m[j] * r6inv;
            ax[i] += s * rx;
            ay[i] += s * ry;
            az[i] += s * rz;
        }
    }
    t1 = wtime();
    l1 += (t1 - t0);

    // Loop 2.
    t0 = wtime();
    for (i = 0; i < V; i += 4) {
        __m128 ax_ = _mm_load_ps(ax+i);
        __m128 vx_ = _mm_load_ps(vx+i);
        __m128 mult_ = _mm_mul_ps(dmp_, _mm_mul_ps(dt_, ax_));
        __m128 add_ = _mm_add_ps(vx_, mult_);
        _mm_store_ps(vx+i, add_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 ax_ = _mm_load_ps(ay+i);
        __m128 vx_ = _mm_load_ps(vy+i);
        __m128 mult_ = _mm_mul_ps(dmp_, _mm_mul_ps(dt_, ax_));
        __m128 add_ = _mm_add_ps(vx_, mult_);
        _mm_store_ps(vy+i, add_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 ax_ = _mm_load_ps(az+i);
        __m128 vx_ = _mm_load_ps(vz+i);
        __m128 mult_ = _mm_mul_ps(dmp_, _mm_mul_ps(dt_, ax_));
        __m128 add_ = _mm_add_ps(vx_, mult_);
        _mm_store_ps(vz+i, add_);
    }
    for (; i < N; i++) {
        vx[i] += dmp * (dt * ax[i]);
        vy[i] += dmp * (dt * ay[i]);
        vz[i] += dmp * (dt * az[i]);
    }

    t1 = wtime();
    l2 += (t1 - t0);

    // Loop 3.
    t0 = wtime();
    for (i = 0; i < V; i += 4) {
        __m128 x_ = _mm_load_ps(x+i);
        __m128 vx_ = _mm_load_ps(vx+i);

        __m128 mult_ = _mm_mul_ps(dt_, vx_);
        x_ = _mm_add_ps(x_, mult_);
        _mm_store_ps(x+i, x_);

        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 le_ = _mm_cmple_ps(x_, negone_);
        __m128 cond_ = _mm_or_ps(ge_, le_);
        __m128  and_ = _mm_and_ps(cond_, negzero_);
        __m128 invx_ = _mm_xor_ps(vx_, and_);
        _mm_store_ps(vx+i, invx_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 x_ = _mm_load_ps(y+i);
        __m128 vx_ = _mm_load_ps(vy+i);

        __m128 mult_ = _mm_mul_ps(dt_, vx_);
        x_ = _mm_add_ps(x_, mult_);
        _mm_store_ps(y+i, x_);

        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 le_ = _mm_cmple_ps(x_, negone_);
        __m128 cond_ = _mm_or_ps(ge_, le_);
        __m128  and_ = _mm_and_ps(cond_, negzero_);
        __m128 invx_ = _mm_xor_ps(vx_, and_);
        _mm_store_ps(vy+i, invx_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 x_ = _mm_load_ps(z+i);
        __m128 vx_ = _mm_load_ps(vz+i);

        __m128 mult_ = _mm_mul_ps(dt_, vx_);
        x_ = _mm_add_ps(x_, mult_);
        _mm_store_ps(z+i, x_);

        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 le_ = _mm_cmple_ps(x_, negone_);
        __m128 cond_ = _mm_or_ps(ge_, le_);
        __m128  and_ = _mm_and_ps(cond_, negzero_);
        __m128 invx_ = _mm_xor_ps(vx_, and_);
        _mm_store_ps(vz+i, invx_);
    }
    for (; i < N; i++) {
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
