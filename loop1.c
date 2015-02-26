    float *tx, *ty, *s;
    tx = (float*) _mm_malloc(N * sizeof(float), 16);
    ty = (float*) _mm_malloc(N * sizeof(float), 16);
    s = (float*) _mm_malloc(N * sizeof(float), 16);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < V; j += 4) {
            __m128 s_ = _mm_load_ps(s+j);
            __m128 x1_ = _mm_load1_ps(x+i);
            __m128 x_ = _mm_load_ps(x+j);
            __m128 sub_ = _mm_sub_ps(x_, x1_);
            __m128 mul_ = _mm_mul_ps(sub_, sub_);
            _mm_store_ps(tx+j, sub_);
            _mm_store_ps(s+j, mul_);
        }
        for (int j = 0; j < N; j++) {
            ty[j] = y[j] - y[i];
            s[j] += ty[j]*ty[j];
        }
        for (int j = 0; j < N; j++) {
            float tz = z[j] - z[i];
            s[j] += tz*tz + eps;
            s[j] = 1.f/sqrt(s[j]);
            s[j] = m[j]*s[j]*s[j]*s[j];
            az[i] += s[j] * tz;
        }
        for (int j = 0; j < N; j++) {
            ax[i] += s[j] * tx[j];
        }
        for (int j = 0; j < N; j++) {
            ay[i] += s[j] * ty[j];
        }
    }
