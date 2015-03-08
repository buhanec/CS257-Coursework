    split up r2 after loaing
    for (int i = 0; i < N; i++) {
        __m128 x1_ = _mm_load1_ps(x+i);
        __m128 y1_ = _mm_load1_ps(y+i);
        __m128 z1_ = _mm_load1_ps(z+i);
        for (j = 0; j < V; j += 4) {
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
            __m128 r6inv_ = _mm_mul_ps(_mm_mul_ps(r2inv_, r2inv_), r2inv_);
            __m128 m_ = _mm_load_ps(m+j);
            __m128 s_ = _mm_mul_ps(m_, r6inv_);
            __m128 ax_ = _mm_load1_ps(ax+i);
            __m128 mx_ = _mm_mul_ps(s_, rx_);
            __m128 sx1_ = _mm_add_ps(mx_, _mm_movehl_ps(mx_, mx_));
            __m128 sx2_ = _mm_add_ps(sx1_, _mm_shuffle_ps(sx1_, sx1_, 1));
            _mm_store_ss(ax+i, _mm_add_ss(ax_, sx2_));
            __m128 my_ = _mm_mul_ps(s_, ry_);
            __m128 sy1_ = _mm_add_ps(my_, _mm_movehl_ps(my_, my_));
            __m128 sy2_ = _mm_add_ps(sy1_, _mm_shuffle_ps(sy1_, sy1_, 1));
            __m128 ay_ = _mm_load1_ps(ay+i);
            _mm_store_ss(ay+i, _mm_add_ss(ay_, sy2_));
            __m128 mz_ = _mm_mul_ps(s_, rz_);
            __m128 sz1_ = _mm_add_ps(mz_, _mm_movehl_ps(mz_, mz_));
            __m128 sz2_ = _mm_add_ps(sz1_, _mm_shuffle_ps(sz1_, sz1_, 1));
            __m128 az_ = _mm_load1_ps(az+i);
            _mm_store_ss(az+i, _mm_add_ss(az_, sz2_));
        }
        //*/
        /*
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
        */
    }
