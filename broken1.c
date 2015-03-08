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
    float *t;
    t = (float*) _mm_malloc(4 * sizeof(float), 16);
    for (int i = 0; i < N; i++) {
        __m128 x1_ = _mm_load1_ps(x+i);
        __m128 y1_ = _mm_load1_ps(y+i);
        __m128 z1_ = _mm_load1_ps(z+i);
        for (j = 0; j < V; j += 4) {
            float rx = x[j] - x[i];
            float ry = y[j] - y[i];
            float rz = z[j] - z[i];
            float r2 = rx*rx + ry*ry + rz*rz + eps;
            __m128 x_ = _mm_load_ps(x+j);
            __m128 rx_ = _mm_sub_ps(x_, x1_);
            __m128 y_ = _mm_load_ps(y+j);
            __m128 ry_ = _mm_sub_ps(y_, y1_);
            __m128 z_ = _mm_load_ps(z+j);
            __m128 rz_ = _mm_sub_ps(z_, z1_);
            __m128 x2_ = _mm_mul_ps(rx_, rx_);
            __m128 y2_ = _mm_mul_ps(ry_, ry_);
            __m128 z2_ = _mm_mul_ps(rz_, rz_);
            __m128 r2_ = _mm_add_ps(_mm_add_ps(_mm_add_ps(x2_, y2_), z2_), eps_);
            // Newton-Phelps + fast inverse
            __m128 inv_ = _mm_rsqrt_ps(r2_);
            __m128 top_ = _mm_mul_ps(_mm_mul_ps(r2_, inv_), inv_);
            __m128 r2inv_ = _mm_mul_ps(_mm_mul_ps(half_, inv_), _mm_sub_ps(three_, top_));
            // Accurate inverse square root
            //float r2inv = 1.0f / sqrt(r2);
            __m128 r2inv_ =_mm_div_ps(one_, _mm_sqrt_ps(r2_));
            //__m128 r2inv_ = _mm_rsqrt_ps(r2_);
            // Magic bullshit
            //__m128 r2inv_ = _mm_rsqrt_ps(r2_);
            __m128 r6inv_ = _mm_mul_ps(_mm_mul_ps(r2inv_, r2inv_), r2inv_);
            __m128 m_ = _mm_load_ps(m+j);
            __m128 s_ = _mm_mul_ps(m_, r6inv_);
            __m128 mx_ = _mm_mul_ps(s_, rx_);
            __m128 sx_ = _mm_hsum_ps(mx_);
            __m128 ax_ = _mm_load1_ps(ax+i);
            _mm_store_ss(ax+i, _mm_add_ss(ax_, sx_));
            __m128 my_ = _mm_mul_ps(s_, ry_);
            __m128 sy_ = _mm_hsum_ps(my_);
            __m128 ay_ = _mm_load1_ps(ay+i);
            _mm_store_ss(ay+i, _mm_add_ss(ay_, sy_));
            __m128 mz_ = _mm_mul_ps(s_, rz_);
            __m128 sz_ = _mm_hsum_ps(mz_);
            __m128 az_ = _mm_load1_ps(az+i);
            _mm_store_ss(az+i, _mm_add_ss(az_, sz_));
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
    _mm_free(t);
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
