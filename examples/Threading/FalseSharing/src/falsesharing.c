#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include "example.h"

#define REPEAT 1000
#define MAXSIZE 4194304
#define GAPSIZE 32

void timeKernel(int* a,int* result, int n, int gap, void (*pt2Func)(int*,int*,int,int), double* time);


int main(int argc, char** argv)
{
	// Store Time at start, time at end and difference
	double t_start, t_stop;
	int *a;
	int *b;
	int n,r,i;
	int fs_sum,nfs_sum;
	double fs_t,nfs_t;
	int threadcount;
	int gap = GAPSIZE;

  printf("Usage: ./falsesharing <threadcount>\n");
  printf("Will default to 1 thread if not specified.\n\n");

	// Get number of threads to use.
	// OpenMP would normally pick up the OMP_NUM_THREADS environemnt variable
	// Including this here just as a quick means of fixing the upper limit on the
	// number of threads
	if(argc>1)
	{
		threadcount = atoi(argv[1]);
	}
	else
	{
		threadcount=1;
	}

	// Set max thread count
	omp_set_num_threads(threadcount);

#	pragma omp parallel
	{
		if(omp_get_thread_num() == 0)
		{
			int threads = omp_get_num_threads();
			printf("OpenMP Test: Set to %d, OpenMP Reports %d\n",threadcount,threads);
		}
	}

	printf("N,GapSize,FalseSharing,NoFalseSharing\n");

	for(n=1;n<=MAXSIZE;n=n*2)
	{
		fs_sum = 0;	fs_t = 0.0;
		nfs_sum = 0;	nfs_t = 0.0;

		// Setup Arrays
		a = malloc(sizeof(int) * n);
		b = malloc(sizeof(int) * threadcount * gap);

    init_arrays(a,b,n,threadcount,gap);
    timeKernel(a,b,n,1,&nofalsesharing,&fs_t);

		for(i=0;i<threadcount;i++)
		{
			fs_sum += b[i];
		}

    init_arrays(a,b,n,threadcount,gap);
    timeKernel(a,b,n,gap,&nofalsesharing,&nfs_t);

		for(i=0;i<threadcount*gap;i++)
		{
			nfs_sum += b[i];
		}

		if((fs_sum != nfs_sum))
		{
			printf("For Value N, Kernels do not validate\n");
			printf("%d,%d,%d,\n",n,fs_sum,nfs_sum);
		}

		printf("%d,%d,%E,%E\n",n,gap,fs_t,nfs_t);

		// Free Memory
		free(a);
		free(b);

	}

	return 1;
}


void timeKernel(int* a,int* result, int n, int gap, void (*pt2Func)(int*,int*,int,int), double* time)
{
    double t_start = 0.0;
	  double t_stop = 0.0;
    double t = 0.0;
    int r;

		t_start = timer();
	  for(r=0;r<REPEAT;r++)
	  {
      (*pt2Func)(a,result,n,gap);
    }
    t_stop = timer();
    t = (t_stop - t_start)/REPEAT;

    *time = t;


}


