    // Preponderation
    float factor = dmp * dt;
    int V = (N/4)*4;
    int B = N/4;
    int i, j, k, l;

    // packaged floats
    __m128 dmp_ = _mm_load1_ps(&dmp);
    __m128 dt_ = _mm_load1_ps(&dt);
    __m128 factor_ = _mm_mul_ps(dmp_, dt_);
    __m128 zero_ = _mm_setzero_ps();
    __m128 negzero_ = _mm_set1_ps(-0.0f);
    __m128 one_ = _mm_set1_ps(1.0f);
    __m128 negone_ = _mm_set1_ps(-1.0f);
    __m128 fsign_ = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
