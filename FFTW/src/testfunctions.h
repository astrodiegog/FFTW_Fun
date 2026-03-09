#include <stdlib.h>
#include <math.h>

#include <complex.h>

#include <fftw3.h>


/* \fn void TestFunctionOne(int, double, fftw_complex, float) */
/* Routine for testing fftw */
extern void TestFunctionOne(int Nx, double *x_arr, fftw_complex *y_arr, double a);

