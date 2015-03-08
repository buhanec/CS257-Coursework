#include "example.h"
#include <stdlib.h>

void nonblocking(float* a, float* b, float* c, int n, int empty)
{
	int i,j,k;

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			for(k=0;k<n;k++)
			{
				c[j+(i*n)] += a[k + (i*n)] * b[j + (k*n)];		
			}
		}
	}
}
