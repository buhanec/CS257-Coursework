#include <stdlib.h>
#include <stdio.h>
#include "example.h"

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
	double nb_sum,b_sum;
	double nb_t,b_t;

	printf("N,NonBlocking,Blocking\n");

	for(n=16;n<=MAXSIZE;n=n*2)
	{
		nb_sum = 0.0;	nb_t = 0.0;
		b_sum = 0.0;	b_t = 0.0;

		// Setup Arrays
		a = malloc(sizeof(float) * n * n);
		b = malloc(sizeof(float) * n * n);
		c = malloc(sizeof(float) * n * n);

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


		if((nb_sum != b_sum))
		{
			printf("-----\n");
			printf("For Value N, Kernels do not validate\n");
			printf("%d,%E,%E\n",n,nb_sum,b_sum);
			printf("-----\n");

		}

		printf("%d,%E,%E\n",n,nb_t,b_t);


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

