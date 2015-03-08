
#include "example.h"

void init_arrays(float* a, float* b, float* c, int n)
{
    int i;

		for(i=0;i<n;i++)
		{
				a[i] = ((i % 5)+1) * 0.0573;
				b[i] = ((i % 3)+1) * 0.0346;
				c[i] = 0;
		}
}
