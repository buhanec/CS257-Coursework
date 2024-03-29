#include <immintrin.h>

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

__m128 _mm_hsum_ps(__m128 a) {
    __m128 b = _mm_add_ps(a, _mm_movehl_ps(a, a));
    return _mm_add_ps(b, _mm_shuffle_ps(b, b, 1));
}

// - prior to merging loops l0 had the highest LLC miss ratio accounting for
//   40% +- 4% of the cache misses
// - most likely due to extremely memory-bound operation
// - after merging loops 1+2 and 3+4, new test performed
// - l1+2 now relatively stall-free and cache misses are not that common
// - l3+4 accounts for 80% of the cache misses

// http://www.agner.org/optimize/instruction_tables.pdf
// Intel Intrinsics guide
// Some reference for square root precisions
// Some reference for order of operations precision
// Some other reference for something else

// after removing l0 and further optimising l1, l3+4 had 80% of the cache
// misses, mostly stemming from comparisons
// additionally the 3+4

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
        for (int j = 0; j < V; j += 4) {
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
        // Horizontal sum - source of error
        //__m128 vx_  = _mm_load1_ps(vx+i);
        //       sx_  = _mm_add_ps(sx_, vx_);
        __m128 sx1_ = _mm_add_ps(sx_, _mm_movehl_ps(sx_, sx_));
        __m128 sx2_ = _mm_add_ps(sx1_, _mm_shuffle_ps(sx1_, sx1_, 1));
        _mm_store_ss(ax+i, sx2_);
        //__m128 vy_  = _mm_load1_ps(vy+i);
        //       sy_  = _mm_add_ps(sy_, vy_);
        __m128 sy1_ = _mm_add_ps(sy_, _mm_movehl_ps(sy_, sy_));
        __m128 sy2_ = _mm_add_ps(sy1_, _mm_shuffle_ps(sy1_, sy1_, 1));
        _mm_store_ss(ay+i, sy2_);
        //__m128 vz_  = _mm_load1_ps(vz+i);
        //       sz_  = _mm_add_ps(sz_, vz_);
        __m128 sz1_ = _mm_add_ps(sz_, _mm_movehl_ps(sz_, sz_));
        __m128 sz2_ = _mm_add_ps(sz1_, _mm_shuffle_ps(sz1_, sz1_, 1));
        _mm_store_ss(az+i, sz2_);
    }
    t1 = wtime();
    l1 += (t1 - t0);


    // Loop 3.
    t0 = wtime();
    for (i = 0; i < V; i += 4) {
        __m128 ax_ = _mm_load_ps(ax+i);
        __m128 vx_ = _mm_load_ps(vx+i);
        __m128 m1_ = _mm_mul_ps(dmp_, _mm_mul_ps(dt_, ax_));
               vx_ = _mm_add_ps(vx_, m1_);

        __m128 x_  = _mm_load_ps(x+i);
        __m128 m_  = _mm_mul_ps(dt_, vx_);
               x_  = _mm_add_ps(x_, m_);
        _mm_store_ps(x+i, x_);

               x_  = _mm_andnot_ps(negzero_, x_);
        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 mx_ = _mm_and_ps(ge_, negzero_);
               vx_ = _mm_xor_ps(vx_, mx_);
        _mm_store_ps(vx+i, vx_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 ax_ = _mm_load_ps(ay+i);
        __m128 vx_ = _mm_load_ps(vy+i);
        __m128 m1_ = _mm_mul_ps(dmp_, _mm_mul_ps(dt_, ax_));
               vx_ = _mm_add_ps(vx_, m1_);

        __m128 x_  = _mm_load_ps(y+i);
        __m128 m_  = _mm_mul_ps(dt_, vx_);
               x_  = _mm_add_ps(x_, m_);
        _mm_store_ps(y+i, x_);

               x_  = _mm_andnot_ps(negzero_, x_);
        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 mx_ = _mm_and_ps(ge_, negzero_);
               vx_ = _mm_xor_ps(vx_, mx_);
        _mm_store_ps(vy+i, vx_);
    }
    for (i = 0; i < V; i += 4) {
        __m128 ax_ = _mm_load_ps(az+i);
        __m128 vx_ = _mm_load_ps(vz+i);
        __m128 m1_ = _mm_mul_ps(dmp_, _mm_mul_ps(dt_, ax_));
               vx_ = _mm_add_ps(vx_, m1_);

        __m128 x_  = _mm_load_ps(z+i);
        __m128 m_  = _mm_mul_ps(dt_, vx_);
               x_  = _mm_add_ps(x_, m_);
        _mm_store_ps(z+i, x_);

               x_  = _mm_andnot_ps(negzero_, x_);
        __m128 ge_ = _mm_cmpge_ps(x_, one_);
        __m128 mx_ = _mm_and_ps(ge_, negzero_);
               vx_ = _mm_xor_ps(vx_, mx_);
        _mm_store_ps(vz+i, vx_);
    }
    t1 = wtime();
    l3 += (t1 - t0);
}
