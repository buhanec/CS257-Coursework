#include <stdlib.h>

int rowfirst(int* a, int m, int n)
{
	int i,j;
	int sum = 0;
	for(j=0;j<n;j++)
	{
		for(i=0;i<m;i++)
		{
			sum += a[j + (i*n)];
		}
	}
	return sum;
}
