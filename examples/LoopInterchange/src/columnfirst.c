#include <stdlib.h>

int columnfirst(int* a, int m, int n)
{
	int i,j;
	int sum = 0;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{

			sum += a[j + (i*n)];
		}
	}
	return sum;
}
