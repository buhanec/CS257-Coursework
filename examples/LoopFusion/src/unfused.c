#include "example.h"

void unfused(int* a, int* b, int* c, int* d, int n)
{
	int i;

	for(i=0;i<n;i++)
	{
		b[i] = a[i]*2;
	}
	
	for(i=0;i<n;i++)
	{
		c[i] = b[i] + 4;
	}

	for(i=1;i<n;i++)
	{
		d[i] = c[i-1] - 5;
	}
}
