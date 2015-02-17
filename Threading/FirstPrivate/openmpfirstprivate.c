#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char** argv)
{
	int a = 14;
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

	printf("Start:a=%d\n",a);

#	pragma omp parallel default(none) firstprivate(a)
	{
		int tid = omp_get_thread_num();
		if(tid == 0)
		{
			a = a + 4;
			printf("T%d:a=%d\n",tid,a);
		}
		else if(tid == 1)
		{
			a = a + 5;
			printf("T%d:a=%d\n",tid,a);
		}
	}

	printf("Stop:a=%d\n",a);

	return 0;
}
