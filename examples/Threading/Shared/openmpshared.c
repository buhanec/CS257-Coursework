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

	printf("Start:a=%d,b=%d\n",a,b);

#	pragma omp parallel default(none) shared(a,b)
	{
		int tid = omp_get_thread_num();
		if(tid == 0)
		{
			a = a + 4;
			printf("T%d:a=%d\n",tid,a);
		}
		else if(tid==1)
		{
			b = b + 5;
			printf("T%d:b=%d\n",tid,b);
		}
	}

	printf("Stop:a=%d,b=%d\n",a,b);

	return 0;
}
