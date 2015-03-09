#include <immintrin.h>

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

__m128 _mm_hsum_ps(__m128 a) {
    __m128 b = _mm_add_ps(a, _mm_movehl_ps(a, a));
    return _mm_add_ps(b, _mm_shuffle_ps(b, b, 1));
}

void compute() {
    // Preponderationing
    float factor = dmp * dt;
    int V = (N/4)*4;
    int i, j;

    // Packed preponderationing
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

    // Timers
    double t0, t1;

    // Acceleration and velocity loop
    t0 = wtime();
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        float *t1;
        float *t;
        t = (float*) _mm_malloc(4 * sizeof(float), 16);
        // prep x
        __m128 xi_ = _mm_load1_ps(x+i);
        __m128 sx_ = zero_;
        // prep y
        __m128 yi_ = _mm_load1_ps(y+i);
        __m128 sy_ = zero_;
        // prep z
        __m128 zi_ = _mm_load1_ps(z+i);
        __m128 sz_ = zero_;
        // the nasty part
        for (j = 0; j < V; j += 4) {
            // rx/x2
            __m128 rx_ = _mm_load_ps(x+j);
                   rx_ = _mm_sub_ps(rx_, xi_);
            __m128 r2_ = _mm_mul_ps(rx_, rx_);
            // r2 start
                   r2_ = _mm_add_ps(r2_, eps_);
            // ry/y2
            __m128 ry_ = _mm_load_ps(y+j);
                   ry_ = _mm_sub_ps(ry_, yi_);
            __m128 y2_ = _mm_mul_ps(ry_, ry_);
            // r2 con
                   r2_ = _mm_add_ps(y2_, r2_);
            //rz/z2
            __m128 rz_ = _mm_load_ps(z+j);
                   rz_ = _mm_sub_ps(rz_, zi_);
            __m128 z2_ = _mm_mul_ps(rz_, rz_);
            // r2 fin
                   r2_ = _mm_add_ps(z2_, r2_);
            // s start
            __m128 s_ = _mm_load_ps(m+j);
            // Fast inverse - source of error
            //__m128 r2inv_ = _mm_rsqrt_ps(r2_);
            // Newton-Raphson step - source of error
            //__m128 top_ = _mm_mul_ps(_mm_mul_ps(r2_, r2inv_), r2inv_);
            //__m128 r2inv_ = _mm_mul_ps(_mm_mul_ps(half_, r2inv_), _mm_sub_ps(three_, top_));
            // Accurate inverse square root
            __m128 r2inv_ = _mm_div_ps(one_, _mm_sqrt_ps(r2_));
            // r6inv
            __m128 r6inv_ = _mm_mul_ps(_mm_mul_ps(r2inv_, r2inv_), r2inv_);
            // s fin
                   s_ = _mm_mul_ps(s_, r6inv_);
            // Directly calculate velocity - source of error
            //       s_  = _mm_mul_ps(factor_, s_);
            // Calculate results
            __m128 mx_ = _mm_mul_ps(s_, rx_);
                   sx_ = _mm_add_ps(mx_, sx_);
            __m128 my_ = _mm_mul_ps(s_, ry_);
                   sy_ = _mm_add_ps(my_, sy_);
            __m128 mz_ = _mm_mul_ps(s_, rz_);
                   sz_ = _mm_add_ps(mz_, sz_);

        }
        // Horizontal sum - source of error
        __m128 vx_  = _mm_load1_ps(vx+i);
               sx_  = _mm_add_ps(sx_, vx_);
        __m128 sx1_ = _mm_add_ps(sx_, _mm_movehl_ps(sx_, sx_));
        __m128 sx2_ = _mm_add_ps(sx1_, _mm_shuffle_ps(sx1_, sx1_, 1));
        _mm_store_ss(ax+i, sx2_);
        __m128 vy_  = _mm_load1_ps(vy+i);
               sy_  = _mm_add_ps(sy_, vy_);
        __m128 sy1_ = _mm_add_ps(sy_, _mm_movehl_ps(sy_, sy_));
        __m128 sy2_ = _mm_add_ps(sy1_, _mm_shuffle_ps(sy1_, sy1_, 1));
        _mm_store_ss(ay+i, sy2_);
        __m128 vz_  = _mm_load1_ps(vz+i);
               sz_  = _mm_add_ps(sz_, vz_);
        __m128 sz1_ = _mm_add_ps(sz_, _mm_movehl_ps(sz_, sz_));
        __m128 sz2_ = _mm_add_ps(sz1_, _mm_shuffle_ps(sz1_, sz1_, 1));
        _mm_store_ss(az+i, sz2_);
    }

    t1 = wtime();
    l0 += (t1 - t0);


    // Position Loops
    t0 = wtime();
    for (i = 0; i < V; i += 4) {
        __m128 vx_ = _mm_load_ps(vx+i);
        __m128 x_  = _mm_load_ps(x+i);
        __m128 m_  = _mm_mul_ps(dt_, vx_);
               x_ = _mm_add_ps(vx_, m_);
        _mm_store_ps(x+i, x_);

        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 le_ = _mm_cmple_ps(x_, negone_);
        __m128 cond_ = _mm_or_ps(ge_, le_);
        __m128  and_ = _mm_and_ps(cond_, negzero_);
        __m128 invx_ = _mm_xor_ps(vx_, and_);
        _mm_store_ps(vx+i, invx_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 vx_ = _mm_load_ps(vy+i);
        __m128 x_  = _mm_load_ps(y+i);
        __m128 m_  = _mm_mul_ps(dt_, vx_);
               x_ = _mm_add_ps(vx_, m_);
        _mm_store_ps(y+i, x_);

        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 le_ = _mm_cmple_ps(x_, negone_);
        __m128 cond_ = _mm_or_ps(ge_, le_);
        __m128  and_ = _mm_and_ps(cond_, negzero_);
        __m128 invx_ = _mm_xor_ps(vx_, and_);
        _mm_store_ps(vy+i, invx_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 vx_ = _mm_load_ps(vz+i);
        __m128 x_  = _mm_load_ps(z+i);
        __m128 m_  = _mm_mul_ps(dt_, vx_);
               x_ = _mm_add_ps(vx_, m_);
        _mm_store_ps(z+i, x_);

        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 le_ = _mm_cmple_ps(x_, negone_);
        __m128 cond_ = _mm_or_ps(ge_, le_);
        __m128  and_ = _mm_and_ps(cond_, negzero_);
        __m128 invx_ = _mm_xor_ps(vx_, and_);
        _mm_store_ps(vz+i, invx_);
    }
    for (; i < N; i++) {
        vx[i] += dmp * (dt * ax[i]);
        vy[i] += dmp * (dt * ay[i]);
        vz[i] += dmp * (dt * az[i]);
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
