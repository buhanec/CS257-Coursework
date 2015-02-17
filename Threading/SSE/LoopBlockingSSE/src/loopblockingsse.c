#include <stdlib.h>
#include <stdio.h>
#include "example.h"
#include <malloc.h>

#define REPEAT 2
#define MAXSIZE 2048
#define BLOCKSIZE 32

void timeKernel(float* a,float* b, float* c,int n, void (*pt2Func)(float*,float*,float*,int,int),double* time);

int main(int argc, char* argv)
{
	// Store Time at start, time at end and difference
	double t_start, t_stop;
	float *a, *b, *c;
	int n,r,j,i;
	double nb_sum,b_sum,br_sum,brv_sum;
	double nb_t,b_t,br_t,brv_t;

	printf("N,NonBlocking,Blocking,BlockingUnroll4,BlockingSSE\n");

	for(n=8;n<=MAXSIZE;n=n*2)
	{
		nb_sum = 0.0;	nb_t = 0.0;
		b_sum = 0.0;	b_t = 0.0;
		br_sum = 0.0;	br_t = 0.0;
    brv_sum = 0.0; brv_t = 0.0;

		// Setup Arrays
		a = memalign(64,sizeof(float) * n * n);
		b = memalign(64,sizeof(float) * n * n);
		c = memalign(64,sizeof(float) * n * n);

    // Time NonBlocking
    init_arrays(a,b,c,n);
    timeKernel(a,b,c,n,&nonblocking,&nb_t);

		// Gen Validation
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				nb_sum += c[j + (i*n)];
			}
		}

    // Time Blocking
    init_arrays(a,b,c,n);
    timeKernel(a,b,c,n,&blocking,&b_t);

		// Gen Validation
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				b_sum += c[j + (i*n)];
			}
		}

    // Time Blocking Unroll
    init_arrays(a,b,c,n);
    timeKernel(a,b,c,n,&blocking_unroll,&br_t);

		// Gen Validation
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				br_sum += c[j + (i*n)];
			}
		}


    // Time Blocking SSE
    init_arrays(a,b,c,n);
    timeKernel(a,b,c,n,&blocking_sse,&brv_t);

		// Gen Validation
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				brv_sum += c[j + (i*n)];
			}
		}



		if((nb_sum != b_sum) || (nb_sum != br_sum) || (nb_sum != brv_sum) )
		{
			printf("-----\n");
			printf("For Value N, Kernels do not validate\n");
			printf("%d,%E,%E,%E,%E\n",n,nb_sum,b_sum,br_sum,brv_sum);
			printf("-----\n");

		}

		printf("%d,%E,%E,%E,%E\n",n,nb_t,b_t,br_t,brv_t);


		// Free Memory
		free(a);
		free(b);
		free(c);
	}
}



void timeKernel(float* a,float* b, float* c,int n, void (*pt2Func)(float*,float*,float*,int,int),double* time)
{
    double t_start = 0.0;
	  double t_stop = 0.0;
    double t = 0.0;
    int r;

		t_start = timer();
	  for(r=0;r<REPEAT;r++)
	  {
      (*pt2Func)(a,b,c,n,BLOCKSIZE);
    }
    t_stop = timer();
    t = (t_stop - t_start)/REPEAT;

    *time = t;
}

