#ifndef EXAMPLE_H
#define EXAMPLE_H

#define min(a,b) (a) < (b) ? (a) : (b)
#define max(a,b) (a) > (b) ? (a) : (b)

double timer();
void init_arrays(float* a, float* b, float* c, int n);
void nonblocking(float* a, float* b, float* c, int n, int empty);
void blocking(float* a, float* b, float* c, int n, int block);

#endif
