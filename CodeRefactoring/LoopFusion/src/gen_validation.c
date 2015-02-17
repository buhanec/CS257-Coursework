#include "example.h"

int genVal(int* d,int n)
{
  int i;
  int sum=0;

	for(i=1;i<n;i++)
	{
		sum += d[i];
  }

  return sum;
}
