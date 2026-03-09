#include <stdlib.h>
#include <math.h>

#include <complex.h>

#include <fftw3.h>


/* \fn void TestFunctionOne(int, double, fftw_complex, float) */
/* Routine for testing fftw */
extern void TestFunctionOne(int Nx, double *x_arr, fftw_complex *y_arr, double a);


/* \fn void TestFunctionOne_FFT(int, double, fftw_complex, float) */
/* Routine for analytic fft */
extern void TestFunctionOne_FFT(int Nx, double *kx_arr, fftw_complex *FFT_analytic, double a);


