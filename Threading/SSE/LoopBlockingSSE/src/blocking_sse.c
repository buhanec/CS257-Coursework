#include "example.h"
#include <stdlib.h>
#include <immintrin.h>

void blocking_sse(float* a, float* b, float* c, int n, int block)
{

	int i,j,k;
	int i0,j0,k0;
  __m128 vec_a, vec_b, vec_c;

	for(i0=0;i0<n;i0+=block)
	{
		for(j0=0;j0<n;j0+=block)
		{
			for(k0=0;k0<n;k0+=block)
			{
				for(i=i0;i< (min(i0+block,n));i++)
				{

          int vecsize = ((min(j0+block,n)) / 4 ) * 4;

					for(j=j0;j<vecsize;j+=4)
					{
            vec_c = _mm_load_ps(&(c[j+(i*n)]));

						for(k=k0;k< (min(k0+block,n));k++)
						{
              vec_a = _mm_set1_ps(a[k + (i*n)]);
              vec_b = _mm_load_ps(&(b[j + (k*n)]));
              vec_c = _mm_add_ps(vec_c,_mm_mul_ps(vec_a,vec_b));
						}

            _mm_store_ps(&(c[j+(i*n)]),vec_c);
					}

					for(;j<(min(j0+block,n));j++)
					{
						for(k=k0;k< (min(k0+block,n));k++)
						{
							c[j+(i*n)] += a[k + (i*n)] * b[j + (k*n)];
						}
					}


				}
			}
		}
	}
}
