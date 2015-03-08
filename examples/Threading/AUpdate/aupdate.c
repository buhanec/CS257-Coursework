#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char** argv)
{
	int a = 14;
	int b = 12;
	int tcount;

	printf("Usage:./openmpprivate <threadcount>\n");
	printf("Demonstrates race condition without atomic operation\n");
	printf("Demonstrates no race condition with atomic operation\n");
	printf("Note: Race condition easier to see with more threads (can go higher than core count)\n");
	printf("May need to run a few times to see it - with race condition number of threads is incorrect\n");

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
	int count=0;
#	pragma omp parallel
	{
#		pragma omp atomic
		count++;
	}
	printf("Number of Threads, Atomic:%d\n",count);

	count = 0;
#	pragma omp parallel
	{
		count++;
	}
	printf("Number of Threads, Non-Atomic:%d\n",count);


	return 0;
}
