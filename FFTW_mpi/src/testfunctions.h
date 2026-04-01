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


/* \fn fftw_complex TestFunctionFour(double, double) */
/* Routine for testing fftw */
fftw_complex TestFunctionFour(double x, double y);

/* \fn fftw_complex TestFunctionFour_FFT(double, double) */
/* Routine for analytic fft */
fftw_complex TestFunctionFour_FFT(double kx, double ky);


/* \fn fftw_complex TestFunctionFive(double, double) */
/* Routine for testing fftw */
fftw_complex TestFunctionFive(double x, double y);

/* \fn fftw_complex TestFunctionFive_FFT(double, double) */
/* Routine for analytic fft */
fftw_complex TestFunctionFive_FFT(double kx, double ky);

/* \fn fftw_complex TestFunctionSix(double, double, double, double) */
/* Routine for testing fftw */
fftw_complex TestFunctionSix(double x, double y, double a, double b);

/* \fn fftw_complex TestFunctionSix_FFT(double, double, double, double) */
/* Routine for analytic fft */
fftw_complex TestFunctionSix_FFT(double kx, double ky, double a, double b);



/* \fn void RunOneDimensionalTest(hid_t, int, double, double) */
/* Routine to run all 1D FFT tests */
extern void RunOneDimensionalTests(hid_t grp_1D_id, int Nx, double xmin, double dx);

/* \fn void RunTwoDimensionalTest(hid_t, int, int, double, double, double, double) */
/* Routine to run all 2D FFT tests */
extern void RunTwoDimensionalTests(hid_t grp_2D_id, int Nx, int Ny, double xmin, double ymin, double dx, double dy);

