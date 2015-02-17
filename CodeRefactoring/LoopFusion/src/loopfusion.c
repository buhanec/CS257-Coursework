#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include "example.h"

#define REPEAT 1000
#define MAXSIZE 2097152

void timeKernel(int* a,int* b, int* c, int* d, int n, void (*pt2Func)(int*,int*,int*,int*,int),double* time);


int main(int argc, char* argv)
{
	// Store Time at start, time at end and difference
	double t_start, t_stop;
	int *a, *b, *c, *d;
	int n,r,j,i;
	int uf_sum,f_sum,fp_sum,fpn_sum;
	double uf_t,f_t,fp_t,fpn_t;

	printf("N,Fused,Unfused,FusedPeeled,FusedPeeledNoInter\n");

	for(n=1;n<=MAXSIZE;n=n*2)
	{
		uf_sum = 0;	uf_t = 0.0;
		f_sum = 0;	f_t = 0.0;
		fp_sum = 0;	fp_t = 0.0;
		fpn_sum = 0;	fpn_t = 0.0;

		// Setup Arrays
		a = malloc(sizeof(int) * n);
		b = malloc(sizeof(int) * n);
		c = malloc(sizeof(int) * n);
		d = malloc(sizeof(int) * n);

		// Time Kernel
    init_arrays(a,b,c,d,n);
    timeKernel(a,b,c,d,n,&unfused,&uf_t);
    uf_sum = genVal(d,n);

    // Time Kernel
    init_arrays(a,b,c,d,n);
    timeKernel(a,b,c,d,n,&fused,&f_t);
    f_sum = genVal(d,n);

    // Time Kernel
    init_arrays(a,b,c,d,n);
    timeKernel(a,b,c,d,n,&fused_peeled,&fp_t);
    fp_sum = genVal(d,n);

    // Time Kernel
    init_arrays(a,b,c,d,n);
    timeKernel(a,b,c,d,n,&fused_peeled_nointer,&fpn_t);
    fpn_sum = genVal(d,n);

    // Validate
		if((uf_sum != f_sum) || (uf_sum != fp_sum) || (uf_sum != fpn_sum))
		{
			printf("For Value N, Kernels do not validate:");
			printf("%d,%d,%d,%d\n",uf_sum,f_sum,fp_sum,fpn_sum);
		}

    // Print Times
		printf("%d,%E,%E,%E,%E\n",n,uf_t,f_t,fp_t,fpn_t);

		// Free Memory
		free(a);
		free(b);
		free(c);
		free(d);
	}
}


void timeKernel(int* a,int* b, int* c, int* d, int n, void (*pt2Func)(int*,int*,int*,int*,int),double* time)
{
    double t_start = 0.0;
	  double t_stop = 0.0;
    double t = 0.0;
    int r;

		t_start = timer();
	  for(r=0;r<REPEAT;r++)
	  {
      (*pt2Func)(a,b,c,d,n);
    }
    t_stop = timer();
    t = (t_stop - t_start)/REPEAT;

    *time = t;


}











