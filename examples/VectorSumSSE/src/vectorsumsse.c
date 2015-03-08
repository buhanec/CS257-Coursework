#include <stdlib.h>
#include <stdio.h>
#include "example.h"
#include <malloc.h>

#define REPEAT 1000
#define MAXSIZE 2091392

void timeKernel(float* a,float* b, float* c,int n, void (*pt2Func)(float*,float*,float*,int),double* time);

int main(int argc, char* argv)
{
	// Store Time at start, time at end and difference
	double t_start, t_stop;
	float *a, *b, *c;
	int n,r,j,i;
	double nonsse_sum,sse_sum,unroll_sum;
	double nonsse_t,sse_t,unroll_t;

	printf("N,NonSSE,Unroll4,SSE\n");

	for(n=1;n<=MAXSIZE;n=n*2)
	{
		nonsse_sum = 0.0;	nonsse_t = 0.0;
		sse_sum = 0.0;	sse_t = 0.0;
		unroll_sum = 0.0;	unroll_t = 0.0;

		// Setup Arrays
		a = memalign(64,sizeof(float) * n);
		b = memalign(64,sizeof(float) * n);
		c = memalign(64,sizeof(float) * n);

    // Time NonSSE
    init_arrays(a,b,c,n);
    timeKernel(a,b,c,n,&nonsse,&nonsse_t);

		// Gen Validation
		for(i=0;i<n;i++)
		{
				nonsse_sum += c[i];
		}


    // Time Unroll
    init_arrays(a,b,c,n);
    timeKernel(a,b,c,n,&unroll,&unroll_t);

		// Gen Validation
		for(i=0;i<n;i++)
		{
		  unroll_sum += c[i];
		}


    // Time SSE
    init_arrays(a,b,c,n);
    timeKernel(a,b,c,n,&sse,&sse_t);

		// Gen Validation
		for(i=0;i<n;i++)
		{
		  sse_sum += c[i];
		}

    // Validation

		if((nonsse_sum != sse_sum) || (nonsse_sum != unroll_sum))
		{
			printf("-----\n");
			printf("For Value N, Kernels do not validate\n");
			printf("%d,%E,%E,%E\n",n,nonsse_sum,sse_sum,unroll_sum);
			printf("-----\n");

		}

		printf("%d,%E,%E,%E\n",n,nonsse_t,unroll_t,sse_t);


		// Free Memory
		free(a);
		free(b);
		free(c);
	}
}



void timeKernel(float* a,float* b, float* c,int n, void (*pt2Func)(float*,float*,float*,int),double* time)
{
    double t_start = 0.0;
	  double t_stop = 0.0;
    double t = 0.0;
    int r;

		t_start = timer();
	  for(r=0;r<REPEAT;r++)
	  {
      (*pt2Func)(a,b,c,n);
    }
    t_stop = timer();
    t = (t_stop - t_start)/REPEAT;

    *time = t;
}

