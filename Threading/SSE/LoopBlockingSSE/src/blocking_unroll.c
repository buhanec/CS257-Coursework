#include "example.h"
#include <stdlib.h>

void blocking_unroll(float* a, float* b, float* c, int n, int block)
{

	int i,j,k;
	int i0,j0,k0;

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
						for(k=k0;k< (min(k0+block,n));k++)
						{
							c[j+(i*n)] += a[k + (i*n)] * b[j + (k*n)];
							c[j+1+(i*n)] += a[k + (i*n)] * b[j + 1 + (k*n)];
							c[j+2+(i*n)] += a[k + (i*n)] * b[j + 2 + (k*n)];
							c[j+3+(i*n)] += a[k + (i*n)] * b[j + 3 + (k*n)];
						}
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
