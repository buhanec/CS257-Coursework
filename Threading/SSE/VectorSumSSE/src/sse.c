#include "example.h"
#include <stdlib.h>
#include <immintrin.h>

void sse(float* a, float* b, float* c, int n)
{

	int i;

  int unroll = (n/4)*4;

  for(i=0;i<unroll;i+=4)
  {
/*
    __m128 vec_a = _mm_load_ps(a);
    __m128 vec_b = _mm_load_ps(b);
    __m128 vec_c = _mm_add_ps(vec_a,vec_b);
    _mm_store_ps(c,vec_c);
*/
    c[i] = a[i] + b[i];
    c[i+1] = a[i+1] + b[i+1];
    c[i+2] = a[i+2] + b[i+2];
    c[i+3] = a[i+3] + b[i+3];

  }

  for(;i<n;i++)
  {
    c[i] = a[i] + b[i];
  }

}
