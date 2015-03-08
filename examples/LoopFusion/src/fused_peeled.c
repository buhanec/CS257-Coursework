#include "example.h"

void fused_peeled(int* a, int* b, int* c, int* d, int n)
{
	int i;

	b[0] = a[0]*2;
	c[0] = b[0] + 4;

	for(i=1;i<n;i++)
	{
		b[i] = a[i]*2;
		c[i] = b[i] + 4;
		d[i] = c[i-1] - 5;
	}
}
