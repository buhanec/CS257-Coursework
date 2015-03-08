#ifndef EXAMPLE_H
#define EXAMPLE_H

double timer();
void init_arrays(int* a, int* b, int n, int threadcount, int gap);
void falsesharing(int* a, int* result, int n, int gap);
void nofalsesharing(int* a, int* result, int n, int gap);

#endif
