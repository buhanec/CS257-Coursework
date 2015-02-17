#include <stdlib.h>
#include <immintrin.h>
#include <stdio.h>

int main(int arc, char** argv)
{

	float* a = memalign(64,4*sizeof(float));

	a[0]=3.2;
	a[1]=2.1;
	a[2]=12.01;
	a[3]=101.85;

	__m128 a_vec = _mm_load_ps(a);
	__m128 temp = _mm_add_ps(a_vec, _mm_movehl_ps(a_vec,a_vec));
	__m128 sum = _mm_add_ss(temp, _mm_shuffle_ps(temp,temp,1));

	// Final result stored in a[0]
	_mm_store_ps(a,sum);

	printf("a[0]:%f,a[1]:%f,a[2]:%f,a[3]:%f\n",a[0],a[1],a[2],a[3]);

	free(a);

	return 1;
}
