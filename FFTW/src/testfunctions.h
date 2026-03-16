#include <stdlib.h>
#include <math.h>

#include <complex.h>

#include "HDF5_utils.h"


#include "hdf5.h"

#include <fftw3.h>




/* \fn void TestFunctionOne(int, double, fftw_complex, float) */
/* Routine for testing fftw */
extern void TestFunctionOne(int Nx, double *x_arr, fftw_complex *y_arr, double a);


/* \fn void TestFunctionOne_FFT(int, double, fftw_complex, float) */
/* Routine for analytic fft */
extern void TestFunctionOne_FFT(int Nx, double *kx_arr, fftw_complex *FFT_analytic, double a);


/* \fn void TestFunctionTwo(int, double, fftw_complex, float) */
/* Routine for testing fftw */
extern void TestFunctionTwo(int Nx, double *x_arr, fftw_complex *y_arr, double a);

/* \fn void TestFunctionTwo_FFT(int, double, fftw_complex, float) */
/* Routine for analytic fft */
extern void TestFunctionTwo_FFT(int Nx, double *kx_arr, fftw_complex *FFT_analytic, double a);



/* \fn void TestFunctionThree(int, double, fftw_complex, float) */
/* Routine for testing fftw */
extern void TestFunctionThree(int Nx, double *x_arr, fftw_complex *y_arr, double a);

/* \fn void TestFunctionThree_FFT(int, double, fftw_complex, float) */
/* Routine for analytic fft */
extern void TestFunctionThree_FFT(int Nx, double *kx_arr, fftw_complex *FFT_analytic, double a);


/* \fn void TestFunctionFour(int, int, double *, double *, fftw_complex *) */
/* Routine for testing fftw */
extern void TestFunctionFour(int Nx, int Ny, double *x_arr, double *y_arr, fftw_complex *fxy_arr);

/* \fn void TestFunctionFour_FFT(int, int, double *, double *, fftw_complex *) */
/* Routine for analytic fft */
extern void TestFunctionFour_FFT(int Nx, int Ny, double *kx_arr, double *ky_arr, fftw_complex *FFT_analytic);


/* \fn void TestFunctionFive(int, int, double *, double *, fftw_complex *) */
/* Routine for testing fftw */
extern void TestFunctionFive(int Nx, int Ny, double *x_arr, double *y_arr, fftw_complex *fxy_arr);

/* \fn void TestFunctionFive_FFT(int, int, double *, double *, fftw_complex *) */
/* Routine for analytic fft */
extern void TestFunctionFive_FFT(int Nx, int Ny, double *kx_arr, double *ky_arr, fftw_complex *FFT_analytic);

/* \fn void TestFunctionSix(int, int, double *, double *, fftw_complex *, double, double ) */
/* Routine for testing fftw */
extern void TestFunctionSix(int Nx, int Ny, double *x_arr, double *y_arr, fftw_complex *fxy_arr, double a, double b);

/* \fn void TestFunctionSix_FFT(int, int, double *, double *, fftw_complex *, double, double) */
/* Routine for analytic fft */
extern void TestFunctionSix_FFT(int Nx, int Ny, double *kx_arr, double *ky_arr, fftw_complex *FFT_analytic, double a, double b);


/* \fn void RunTestOne(hid_t, double *, hid_t, double *, hid_t, double *, int, int, double) */
/* Routine to complete Test One */
extern void RunTestOne(hid_t grp_test_id, double *x_arr, hid_t dataspace_id_c, double *kx_arr_c, hid_t dataspace_id_r, double *kx_arr_r, int Nx, int Nx_r, double a);

/* \fn void RunTestTwo(hid_t, double *, hid_t, double *, hid_t, double *, int, int, double) */
/* Routine to complete Test Two */
extern void RunTestTwo(hid_t grp_test_id, double *x_arr, hid_t dataspace_id_c, double *kx_arr_c, hid_t dataspace_id_r, double *kx_arr_r, int Nx, int Nx_r, double a);

/* \fn void RunTestThree(hid_t, double *, hid_t, double *, hid_t, double *, int, int, double) */
/* Routine to complete Test Three */
extern void RunTestThree(hid_t grp_test_id, double *x_arr, hid_t dataspace_id_c, double *kx_arr_c, hid_t dataspace_id_r, double *kx_arr_r, int Nx, int Nx_r, double a);

/* \fn void RunTestFour(hid_t, double *, double *, hid_t, double *, double *, hid_t, double *, int, int, int, double) */
/* Routine to complete Test Four */
extern void RunTestFour(hid_t grp_test_id, double *x_arr, double *y_arr, hid_t dataspace_id_c, double *kx_arr_c, double *ky_arr_c, hid_t dataspace_id_r, double *ky_arr_r, int Nx, int Ny, int Ny_r);

/* \fn void RunTestFive(hid_t, double *, double *, hid_t, double *, double *, hid_t, double *, int, int, int, double) */
/* Routine to complete Test Five */
extern void RunTestFive(hid_t grp_test_id, double *x_arr, double *y_arr, hid_t dataspace_id_c, double *kx_arr_c, double *ky_arr_c, hid_t dataspace_id_r, double *ky_arr_r, int Nx, int Ny, int Ny_r);

/* \fn void RunTestSix(hid_t, double *, double *, hid_t, double *, double *, hid_t, double *, int, int, int, double) */
/* Routine to complete Test Six */
extern void RunTestSix(hid_t grp_test_id, double *x_arr, double *y_arr, hid_t dataspace_id_c, double *kx_arr_c, double *ky_arr_c, hid_t dataspace_id_r, double *ky_arr_r, int Nx, int Ny, int Ny_r, double a, double b);
