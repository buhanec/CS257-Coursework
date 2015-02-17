#ifndef EXAMPLE_H
#define EXAMPLE_H

double timer();
void unfused(int* a, int* b, int* c, int* d, int n);
void fused(int* a, int* b, int* c, int* d, int n);
void fused_peeled(int* a, int* b, int* c, int* d, int n);
void fused_peeled_nointer(int* a, int* empty, int* empty2,int* d, int n);
void init_arrays(int* a, int* b, int* c, int* d,int n);
int genVal(int* d,int n);

#endif
