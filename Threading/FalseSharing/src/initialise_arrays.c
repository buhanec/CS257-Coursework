#include "example.h"

void init_arrays(int* a, int* b, int n, int threadcount, int gap)
{
    int i;

		for(i=0;i<n;i++)
		{
			a[i] = i % 4;
		}

		for(i=0;i<threadcount*gap;i++)
		{
			b[i] = 0;
		}
}
