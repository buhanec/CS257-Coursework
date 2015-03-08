#ifndef EXAMPLE_H
#define EXAMPLE_H

#define min(a,b) (a) < (b) ? (a) : (b)
#define max(a,b) (a) > (b) ? (a) : (b)

double timer();
void init_arrays(float* a, float* b, float* c, int n);
void sse(float* a, float* b, float* c, int n);
void nonsse(float* a, float* b, float* c, int n);
void unroll(float* a, float* b, float* c, int n);

#endif
