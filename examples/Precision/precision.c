#include <stdlib.h>
#include <stdio.h>

int main(int argc, char** argv)
{
	double a;
	double b;

	a = 1.01;
	b = 2.03;

	double outcome1;
	double outcome2;

	outcome1 = (a+b)/3.0;
	outcome2 = (a/3.0) + (b/3.0);


	printf("(1.01 + 2.03) / 3.0 = %2.20E\n",outcome1);
	printf("(1.01/3.0) + (2.03/3.0) = %2.20E\n",outcome2);
	printf("Difference between answers = %E\n", outcome2 - outcome1);

	return 0;
}
