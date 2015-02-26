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
