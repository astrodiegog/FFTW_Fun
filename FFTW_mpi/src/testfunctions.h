#include <stdlib.h>
#include <math.h>

#include <complex.h>

#include "HDF5_utils.h"


#include "hdf5.h"

#include <fftw3-mpi.h>






/* \fn fftw_complex TestFunctionOne(double, double) */
/* Routine for testing fftw */
fftw_complex TestFunctionOne(double x, double a);

/* \fn fftw_complex TestFunctionOne_FFT(double, double) */
/* Routine for testing analytic fftw */
fftw_complex TestFunctionOne_FFT(double x, double a);


/* \fn fftw_complex TestFunctionTwo(double, double) */
/* Routine for testing fftw */
fftw_complex TestFunctionTwo(double x, double a);

/* \fn fftw_complex TestFunctionTwo_FFT(double, double) */
/* Routine for testing analytic fftw */
fftw_complex TestFunctionTwo_FFT(double x, double a);


/* \fn fftw_complex TestFunctionThree(double, double) */
/* Routine for testing fftw */
fftw_complex TestFunctionThree(double x, double a);

/* \fn fftw_complex TestFunctionThree_FFT(double, double) */
/* Routine for testing analytic fftw */
fftw_complex TestFunctionThree_FFT(double x, double a);


/* \fn void TestFunctionFour(double, double) */
/* Routine for testing fftw */
fftw_complex TestFunctionFour(double x, double y);

/* \fn void TestFunctionFour_FFT(double, double) */
/* Routine for analytic fft */
fftw_complex TestFunctionFour_FFT(double kx, double ky);


/* \fn void TestFunctionFive(double, double) */
/* Routine for testing fftw */
fftw_complex TestFunctionFive(double x, double y);

/* \fn void TestFunctionFive_FFT(double, double) */
/* Routine for analytic fft */
fftw_complex TestFunctionFive_FFT(double kx, double ky);

/* \fn void TestFunctionSix(int, int, double *, double *, fftw_complex *, double, double ) */
/* Routine for testing fftw */
extern void TestFunctionSix(int Nx, int Ny, double *x_arr, double *y_arr, fftw_complex *fxy_arr, double a, double b);

/* \fn void TestFunctionSix_FFT(int, int, double *, double *, fftw_complex *, double, double) */
/* Routine for analytic fft */
extern void TestFunctionSix_FFT(int Nx, int Ny, double *kx_arr, double *ky_arr, fftw_complex *FFT_analytic, double a, double b);



/* \fn void RunTestFive(hid_t, double *, double *, hid_t, double *, double *, hid_t, double *, int, int, int, double) */
/* Routine to complete Test Five */
extern void RunTestFive(hid_t grp_test_id, double *x_arr, double *y_arr, hid_t dataspace_id_c, double *kx_arr_c, double *ky_arr_c, hid_t dataspace_id_r, double *ky_arr_r, int Nx, int Ny, int Ny_r);

/* \fn void RunTestSix(hid_t, double *, double *, hid_t, double *, double *, hid_t, double *, int, int, int, double) */
/* Routine to complete Test Six */
extern void RunTestSix(hid_t grp_test_id, double *x_arr, double *y_arr, hid_t dataspace_id_c, double *kx_arr_c, double *ky_arr_c, hid_t dataspace_id_r, double *ky_arr_r, int Nx, int Ny, int Ny_r, double a, double b);



/* \fn void RunOneDimensionalTest(hid_t, int, double, double) */
/* Routine to run all 1D FFT tests */
extern void RunOneDimensionalTests(hid_t grp_1D_id, int Nx, double xmin, double dx);

/* \fn void RunTwoDimensionalTest(hid_t, int, int, double, double, double, double) */
/* Routine to run all 2D FFT tests */
extern void RunTwoDimensionalTests(hid_t grp_2D_id, int Nx, int Ny, double xmin, double ymin, double dx, double dy);

