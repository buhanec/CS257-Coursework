#include <omp.h>
#include "example.h"

void nofalsesharing(int* a, int* result, int n, int gap)
{
	int i;

#	pragma omp parallel for default(none) shared(a,result,n) private(i) firstprivate(gap)
	for(i=0;i<n;i++)
	{
		int threadid = omp_get_thread_num();
		result[threadid*gap] += a[i];
	}

}
