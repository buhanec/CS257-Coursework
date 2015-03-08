#include <stdlib.h>
#include <stdio.h>
#include "example.h"

#define REPEAT 1000
#define MAXSIZE 32768

void timeKernel(int* a,int m, int n, int (*pt2Func)(int*,int,int), int *sum, double* time);

int main(int argc, char* argv)
{
	// Store Time at start, time at end and difference
	double t_start, t_stop;
	int *a;
	int n,r,j,i;
	int rf_sum,cf_sum;
	double rf_t,cf_t;

	printf("N,RowFirstTraversal,ColumnFirstTraversal\n");

	for(n=8;n<=MAXSIZE;n=n*2)
	{
		rf_sum = 0;	rf_t = 0.0;
		cf_sum = 0;	cf_t = 0.0;

		// Setup Arrays
		a = malloc(sizeof(int) * n * n);

    // Time Different Kernels
    init_arrays(a,n);
    timeKernel(a,n,n,&rowfirst,&rf_sum,&rf_t);

    init_arrays(a,n);
    timeKernel(a,n,n,&columnfirst,&cf_sum,&cf_t);

    // Validate Answers Match
		if((rf_sum != cf_sum))
		{
			printf("For Value N, Kernels do not validate\n");
			printf("%d,%d,%d,\n",n,rf_sum,cf_sum);
		}

    // Print Times
		printf("%d,%E,%E\n",n,rf_t,cf_t);

		// Free Memory
		free(a);
	}
}

void timeKernel(int* a,int m, int n, int (*pt2Func)(int*,int,int), int *sum, double* time)
{
    double t_start = 0.0;
	  double t_stop = 0.0;
    double t = 0.0;
    int r;

		t_start = timer();
	  for(r=0;r<REPEAT;r++)
	  {
      *sum = (*pt2Func)(a,m,n);
    }
    t_stop = timer();
    t = (t_stop - t_start)/REPEAT;

    *time = t;
}


