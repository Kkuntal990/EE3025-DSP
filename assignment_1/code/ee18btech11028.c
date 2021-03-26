#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

double complex *fft(double complex *x, int n)
{
    if (n <= 1)
        return x;

    double complex *X2 = malloc(n / 2 * sizeof(double complex));
    double complex *X1 = malloc(n / 2 * sizeof(double complex));
    for (int i = 0; 2 * i < n; i++)
    {
        X1[i] = x[2 * i];
        X2[i] = x[2 * i + 1];
    }

    X1 = fft(X1, n / 2);
    X2 = fft(X2, n / 2);

    double complex w;
    for (int i = 0; 2 * i < n; i++)
    {
        w = CMPLX(cos(2 * M_PI * i / n), sin(2 * M_PI * i / n));
        x[i] = X1[i] + w * X2[i];
        x[i + n / 2] = X1[i] - w * X2[i];
    }
    free(X1);
    free(X2);

    return x;
}

double complex  *ifft_sub(double complex *x, int n)
{
    if (n <= 1)
        return x;

    double complex *X2 = malloc(n / 2 * sizeof(double complex));
    double complex *X1 = malloc(n / 2 * sizeof(double complex));
    for (int i = 0; 2 * i < n; i++)
    {
        X1[i] = x[2 * i];
        X2[i] = x[2 * i + 1];
    }

    X1 = ifft_sub(X1, n / 2);
    X2 = ifft_sub(X2, n / 2);

    double complex w;
    for (int i = 0; 2 * i < n; i++)
    {
        w = CMPLX(cos(-2 * M_PI * i / n), sin(-2 * M_PI * i / n));
        x[i] = X1[i] + w * X2[i];
        x[i + n / 2] = X1[i] - w * X2[i];
    }
    free(X1);
    free(X2);
    return x;
}

double complex *ifft(double complex *x, int n){
    x = ifft_sub(x, n);
    for (int i = 0; i < n;i++){
        x[i] /= n;
    }

    return x;
}

int main()
{

    int n = 8;

#warning n should be a power of 2

    double complex X[] = {1, 2, 3, 4, 1, 1, 2, 4};
    double complex *xx;
    printf("FFT of x :");
    printf("\n");
    xx = fft(X, n);
    for (int i = 0; i < n; i++)
        printf("(%.3lf, %.3lf)\n", creal(*(xx+i)), cimag(*(xx+i)));

    xx = ifft(xx, n);
    printf("\n\n");
    printf("IFFT of X:");
    printf("\n");
    for (int i = 0; i < n; i++)
        printf("(%.3lf, %.3lf)\n", creal(*(xx + i)), cimag(*(xx + i)));
}