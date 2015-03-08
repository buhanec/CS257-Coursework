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

	printf("Starting Static\n");

#	pragma omp parallel default(none)
	{
		int tid = omp_get_thread_num();
		int i;
#		pragma omp for schedule(static)
		for(i=0;i<10;i++)
		{
			if(tid==0)
			{
				printf("Thread 0 going to sleep\n");
				sleep(5);
				printf("Thread 0 is awake\n");
				printf("T0:i=%d\n",i);
			}
			else
			{
				printf("T%d:i=%d\n",tid,i);
			}
		}
	}

	printf("Starting Dynamic\n");

#	pragma omp parallel default(none)
	{
		int tid = omp_get_thread_num();
		int i;
#		pragma omp for schedule(dynamic)
		for(i=0;i<10;i++)
		{
			if(tid==0)
			{
				printf("Thread 0 going to sleep\n");
				sleep(5);
				printf("Thread 0 is awake\n");
				printf("T0:i=%d\n",i);
			}
			else
			{
				printf("T%d:i=%d\n",tid,i);
			}
		}
	}

	return 0;
}
