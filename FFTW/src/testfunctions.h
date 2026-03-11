#include <stdlib.h>
#include <math.h>

#include <complex.h>

#include "add1D_utils.h"


#include "hdf5.h"

#include <fftw3.h>




/* \fn void TestFunctionOne(int, double, fftw_complex, float) */
/* Routine for testing fftw */
extern void TestFunctionOne(int Nx, double *x_arr, fftw_complex *y_arr, double a);


/* \fn void TestFunctionOne_FFT(int, double, fftw_complex, float) */
/* Routine for analytic fft */
extern void TestFunctionOne_FFT(int Nx, double *kx_arr, fftw_complex *FFT_analytic, double a);


/* \fn void RunTestOne(hid_t, double *, double *, double *, int, double) */
/* Routine to complete Test One */
extern void RunTestOne(hid_t grp_test_id, double *x_arr, hid_t dataspace_id_c, double *kx_arr_c, hid_t dataspace_id_r, double *kx_arr_r, int Nx, int Nx_r, double a);

