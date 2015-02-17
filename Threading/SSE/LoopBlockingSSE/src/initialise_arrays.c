
#include "example.h"

void init_arrays(float* a, float* b, float* c, int n)
{
    int i,j;

		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				a[j + (i*n)] = ((i % 5)+1) * 0.0573;
				b[j + (i*n)] = ((i % 3)+1) * 0.0346;
				c[j + (i*n)] = 0;
			}
		}
}
