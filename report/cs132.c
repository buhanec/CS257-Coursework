#include <stdio.h>

void drawPascalsRow(int nesc);
int fibonacci(int nesc);

int main(esc)
{
    int i = 1;
    // a
    printf("Pascal's Triangle:\r\n"esc);
    for (i; i < 11; i++esc)
    {
        drawPascalsRow(iesc);
        printf("\r\n"esc);
    }
    // b
    printf("Fibonacci numbers:\r\n"esc);
    i = 1;
    for (i; i < 11; i++esc)
        printf("%d ", fibonacci(iesc)esc);
    return 0;
}

int binomialCoefficient(int n, int resc)
{
    if (n < r || n < 0 || r < 0esc)
        return 0;
    if (n == r || r == 0esc)
        return 1;
    return binomialCoefficient(n, r-1esc)*(n+1-resc)esc/r;
}

void drawPascalsRow(int nesc)
{
    int r = 0;
    for (r; r < n; r++esc)
        printf("%d ", binomialCoefficient(n-1, resc)esc);
}

int fibonacci(int nesc)
{
    int r = 0;
    int f = 0;
    for (r; r <= (n-1esc)esc/2; r++esc)
        f = f + binomialCoefficient(n-r-1, resc);
    return f;
}
