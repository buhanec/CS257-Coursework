#include "example.h"
#include <stdlib.h>

void blocking(float* a, float* b, float* c, int n, int block)
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
					for(j=j0;j<(min(j0+block,n));j++)
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
