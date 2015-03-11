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
    const __m128 dmp_     = _mm_load1_ps(&dmp);
    const __m128 dt_      = _mm_load1_ps(&dt);
    const __m128 eps_     = _mm_load1_ps(&eps);
    const __m128 factor_  = _mm_mul_ps(dmp_, dt_);
    const __m128 zero_    = _mm_setzero_ps();
    const __m128 negzero_ = _mm_set1_ps(-0.0f);
    const __m128 one_     = _mm_set1_ps(1.0f);
    const __m128 negone_  = _mm_set1_ps(-1.0f);
    const __m128 half_    = _mm_set1_ps(0.5f);
    const __m128 three_   = _mm_set1_ps(3.0f);

    // Timers
    double t0, t1;

    // Loop 1 (and 0)
    t0 = wtime();
    #pragma omp parallel for
    for (i = 0; i < N; i++) {
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
            __m128 s_  = _mm_load_ps(m+j);
            // Fast inverse - source of error
            __m128 r2inv_ = _mm_rsqrt_ps(r2_);
            // Newton-Raphson step - source of error
            //__m128 top_ = _mm_mul_ps(_mm_mul_ps(r2_, r2inv_), r2inv_);
            //__m128 r2inv_ = _mm_mul_ps(_mm_mul_ps(half_, r2inv_), _mm_sub_ps(three_, top_));
            // Accurate inverse square root
            //__m128 r2inv_ = _mm_div_ps(one_, _mm_sqrt_ps(r2_));
            // r6inv
            __m128 r6inv_ = _mm_mul_ps(_mm_mul_ps(r2inv_, r2inv_), r2inv_);
            // s fin
                   s_  = _mm_mul_ps(s_, r6inv_);
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
        // cleaup loop
        for (; j < N; j ++) {
            // rx/x2
            __m128 rx_ = _mm_load_ss(x+j);
                   rx_ = _mm_sub_ss(rx_, xi_);
            __m128 r2_ = _mm_mul_ss(rx_, rx_);
            // r2 start
                   r2_ = _mm_add_ss(r2_, eps_);
            // ry/y2
            __m128 ry_ = _mm_load_ss(y+j);
                   ry_ = _mm_sub_ss(ry_, yi_);
            __m128 y2_ = _mm_mul_ss(ry_, ry_);
            // r2 con
                   r2_ = _mm_add_ss(y2_, r2_);
            //rz/z2
            __m128 rz_ = _mm_load_ss(z+j);
                   rz_ = _mm_sub_ss(rz_, zi_);
            __m128 z2_ = _mm_mul_ss(rz_, rz_);
            // r2 fin
                   r2_ = _mm_add_ss(z2_, r2_);
            // s start
            __m128 s_  = _mm_load_ss(m+j);
            // Fast inverse - source of error
            __m128 r2inv_ = _mm_rsqrt_ss(r2_);
            // Newton-Raphson step - source of error
            //__m128 top_ = _mm_mul_ss(_mm_mul_ss(r2_, r2inv_), r2inv_);
            //__m128 r2inv_ = _mm_mul_ss(_mm_mul_ss(half_, r2inv_), _mm_sub_ss(three_, top_));
            // Accurate inverse square root
            //__m128 r2inv_ = _mm_div_ss(one_, _mm_sqrt_ss(r2_));
            // r6inv
            __m128 r6inv_ = _mm_mul_ss(_mm_mul_ss(r2inv_, r2inv_), r2inv_);
            // s fin
                   s_  = _mm_mul_ss(s_, r6inv_);
            // Calculate results
            __m128 mx_ = _mm_mul_ss(s_, rx_);
                   sx_ = _mm_add_ss(mx_, sx_);
            __m128 my_ = _mm_mul_ss(s_, ry_);
                   sy_ = _mm_add_ss(my_, sy_);
            __m128 mz_ = _mm_mul_ss(s_, rz_);
                   sz_ = _mm_add_ss(mz_, sz_);
        }
        // Horizontal sum - source of error
        __m128 sx1_ = _mm_add_ps(sx_, _mm_movehl_ps(sx_, sx_));
        __m128 sx2_ = _mm_add_ps(sx1_, _mm_shuffle_ps(sx1_, sx1_, 1));
        _mm_store_ss(ax+i, sx2_);
        __m128 sy1_ = _mm_add_ps(sy_, _mm_movehl_ps(sy_, sy_));
        __m128 sy2_ = _mm_add_ps(sy1_, _mm_shuffle_ps(sy1_, sy1_, 1));
        _mm_store_ss(ay+i, sy2_);
        __m128 sz1_ = _mm_add_ps(sz_, _mm_movehl_ps(sz_, sz_));
        __m128 sz2_ = _mm_add_ps(sz1_, _mm_shuffle_ps(sz1_, sz1_, 1));
        _mm_store_ss(az+i, sz2_);
    }
    t1 = wtime();
    l1 += (t1 - t0);


    // Loop 3 and 4
    t0 = wtime();
    for (i = 0; i < V; i += 4) {
        __m128 ax_ = _mm_load_ps(ax+i);
        __m128 vx_ = _mm_load_ps(vx+i);
        __m128 m_  = _mm_mul_ps(factor_, ax_);
               vx_ = _mm_add_ps(vx_, m_);

        __m128 x_  = _mm_load_ps(x+i);
               m_  = _mm_mul_ps(dt_, vx_);
               x_  = _mm_add_ps(x_, m_);
        _mm_store_ps(x+i, x_);

               x_  = _mm_andnot_ps(negzero_, x_);
        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 mx_ = _mm_and_ps(ge_, negzero_);
               vx_ = _mm_xor_ps(vx_, mx_);
        _mm_store_ps(vx+i, vx_);
    }
    for (; i < N; i++) {
        __m128 ax_ = _mm_load_ss(ax+i);
        __m128 vx_ = _mm_load_ss(vx+i);
        __m128 m_  = _mm_mul_ss(factor_, ax_);
               vx_ = _mm_add_ss(vx_, m_);

        __m128 x_  = _mm_load_ss(x+i);
               m_  = _mm_mul_ss(dt_, vx_);
               x_  = _mm_add_ss(x_, m_);
        _mm_store_ss(x+i, x_);

               x_  = _mm_andnot_ps(negzero_, x_);
        __m128 ge_ = _mm_cmpge_ss(x_, one_);
        __m128 mx_ = _mm_and_ps(ge_, negzero_);
               vx_ = _mm_xor_ps(vx_, mx_);
        _mm_store_ss(vx+i, vx_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 ay_ = _mm_load_ps(ay+i);
        __m128 vy_ = _mm_load_ps(vy+i);
        __m128 m_  = _mm_mul_ps(factor_, ay_);
               vy_ = _mm_add_ps(vy_, m_);

        __m128 y_  = _mm_load_ps(y+i);
               m_  = _mm_mul_ps(dt_, vy_);
               y_  = _mm_add_ps(y_, m_);
        _mm_store_ps(y+i, y_);

               y_  = _mm_andnot_ps(negzero_, y_);
        __m128 ge_ = _mm_cmpge_ps(y_, one_);
        __m128 my_ = _mm_and_ps(ge_, negzero_);
               vy_ = _mm_xor_ps(vy_, my_);
        _mm_store_ps(vy+i, vy_);
    }
    for (; i < N; i++) {
        __m128 ay_ = _mm_load_ss(ay+i);
        __m128 vy_ = _mm_load_ss(vy+i);
        __m128 m_  = _mm_mul_ss(factor_, ay_);
               vy_ = _mm_add_ss(vy_, m_);

        __m128 y_  = _mm_load_ss(y+i);
               m_  = _mm_mul_ss(dt_, vy_);
               y_  = _mm_add_ss(y_, m_);
        _mm_store_ss(y+i, y_);

               y_  = _mm_andnot_ps(negzero_, y_);
        __m128 ge_ = _mm_cmpge_ss(y_, one_);
        __m128 my_ = _mm_and_ps(ge_, negzero_);
               vy_ = _mm_xor_ps(vy_, my_);
        _mm_store_ss(vy+i, vy_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 az_ = _mm_load_ps(az+i);
        __m128 vz_ = _mm_load_ps(vz+i);
        __m128 m_  = _mm_mul_ps(factor_, az_);
               vz_ = _mm_add_ps(vz_, m_);

        __m128 z_  = _mm_load_ps(z+i);
               m_  = _mm_mul_ps(dt_, vz_);
               z_  = _mm_add_ps(z_, m_);
        _mm_store_ps(z+i, z_);

               z_  = _mm_andnot_ps(negzero_, z_);
        __m128 ge_ = _mm_cmpge_ps(z_, one_);
        __m128 mz_ = _mm_and_ps(ge_, negzero_);
               vz_ = _mm_xor_ps(vz_, mz_);
        _mm_store_ps(vz+i, vz_);
    }
    for (; i < N; i++) {
        __m128 az_ = _mm_load_ss(az+i);
        __m128 vz_ = _mm_load_ss(vz+i);
        __m128 m_  = _mm_mul_ss(factor_, az_);
               vz_ = _mm_add_ss(vz_, m_);

        __m128 z_  = _mm_load_ss(z+i);
               m_  = _mm_mul_ss(dt_, vz_);
               z_  = _mm_add_ss(z_, m_);
        _mm_store_ss(z+i, z_);

               z_  = _mm_andnot_ps(negzero_, z_);
        __m128 ge_ = _mm_cmpge_ss(z_, one_);
        __m128 mz_ = _mm_and_ps(ge_, negzero_);
               vz_ = _mm_xor_ps(vz_, mz_);
        _mm_store_ss(vz+i, vz_);
    }
    t1 = wtime();
    l3 += (t1 - t0);
}
