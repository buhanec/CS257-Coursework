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
