#include "example.h"

void fused_peeled_nointer(int* a, int* empty, int* empty2,int* d, int n)
{
	int i;

	for(i=1;i<n;i++)
	{
		d[i] = (a[i-1]*2)+4 - 5;
	}
}
