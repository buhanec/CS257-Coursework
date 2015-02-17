#include "example.h"

void init_arrays(int* a, int* b, int* c, int* d,int n)
{
    int i;
		for(i=0;i<n;i++)
		{
			a[i] = i;
			b[i] = 0;
			c[i] = 0;
			d[i] = 0;
		}
}
