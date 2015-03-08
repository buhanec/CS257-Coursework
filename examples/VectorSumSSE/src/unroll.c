#include "example.h"
#include <stdlib.h>
#include <immintrin.h>

void unroll(float* a, float* b, float* c, int n)
{

	int i;

  int unroll = (n/4)*4;

  for(i=0;i<unroll;i+=4)
  {
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
