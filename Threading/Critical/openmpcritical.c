#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char** argv)
{
	int a = 14;
	int b = 12;
	int tcount;

	if(argc>1)
	{
		tcount = atoi(argv[1]);
	}
	else
	{
		printf("Please call using the format: ./openmpprivate <threadcount>\n");
		return 0;
	}

	omp_set_num_threads(tcount);

#	pragma omp parallel default(none)
	{
		int tid = omp_get_thread_num();

#		pragma omp critical
		{
			printf("I am thread %d\n",tid);
			printf("I am still thread %d\n",tid);
			printf("I am thread %d, about to leave\n",tid);
		}
	}

	return 0;
}
