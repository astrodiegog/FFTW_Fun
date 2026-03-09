#include "testfunctions.h"

extern void TestFunctionOne(int Nx, double *x_arr, fftw_complex *y_arr, double a)
{
	int i;
	double x;
	for (i = 0; i < Nx; i++)
	{
		x = x_arr[i];
		if ( fabs(x) > (1. / (2. * a)) )
		{
			y_arr[i] = 0. + 0. * I;
		}
		else
		{
			y_arr[i] = 1. + 0. * I;
		}
	} 

	return;
}


extern void TestFunctionOne_FFT(int Nx, double *kx_arr, fftw_complex *FFT_analytic, double a)
{
    int i;
    double kx;
    for (i = 0; i < Nx; i++)
    {
        kx = kx_arr[i];
		if (kx == 0)
		{
			FFT_analytic[i] = 1. + 0. * I;
		}
		else
		{
			FFT_analytic[i] = sin((M_PI * kx) / a) / ((kx * M_PI) / a)  + 0. * I;
		}
    }

    return;
}


