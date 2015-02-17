#include "example.h"

void init_arrays(int* a, int n)
{
  int i,j;
  for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			a[j+(i*n)] = i % 4;
		}
	}
}
