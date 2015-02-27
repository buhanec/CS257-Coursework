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
