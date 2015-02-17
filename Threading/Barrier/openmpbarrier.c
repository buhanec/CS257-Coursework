#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char** argv)
{
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
		if(tid==0)
		{
			printf("Thread 0 is going to sleep for 10 seconds\n");
			sleep(10);
			printf("Thread 0 is awake again!\n");
		}
		printf("I'm thread %d\n",tid);

# 		pragma omp barrier
		printf("I'm thread %d after the barrier\n",tid);
	}

	return 0;
}
