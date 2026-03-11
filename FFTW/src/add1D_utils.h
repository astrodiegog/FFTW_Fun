#include <math.h>

#include <complex.h>
#include "hdf5.h"
#include <fftw3.h>



void Write_HDF5_1Dgrouptest(hid_t grp_test_id, char *arr_name, hid_t dataspace_id, double *data_arr);


/* \fn void add_yarr_1Dgrouptest(hid_t, hid_t, fftw_complex *, int) */
/* Routine to add y array */
extern void add_yarr_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *y_arr, int Nx);

/* \fn void add_FFTc2c_1Dgrouptest(hid_t, hid_t, fftw_complex *, int) */
/* Routine to add complex2complex FFT array */
extern void add_FFTc2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *FFT_c2c, int Nx);

/* \fn void add_FFTr2c_1Dgrouptest(hid_t, hid_t, fftw_complex *, int) */
/* Routine to add real2complex FFT array */
extern void add_FFTr2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *FFT_r2c, int Nx);

/* \fn void add_FFTanalyticc2c_1Dgrouptest(hid_t, hid_t, fftw_complex *, int) */
/* Routine to add analytic complex2complex FFT array */
extern void add_FFTanalyticc2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *FFT_analytic_c2c, int Nx);

/* \fn void add_FFTanalyticc2c_1Dgrouptest(hid_t, hid_t, fftw_complex *, int) */
/* Routine to add analytic real2complex FFT array */
extern void add_FFTanalyticr2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *FFT_analytic_r2c, int Nx);

/* \fn void add_iFFTc2c_1Dgrouptest(hid_t, hid_t, fftw_complex *, int) */
/* Routine to add inverse complex2complex FFT array */
extern void add_iFFTc2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *iFFT_c2c, int Nx);

/* \fn void add_iFFTc2r_1Dgrouptest(hid_t, hid_t, fftw_complex *) */
/* Routine to add inverse complex2real FFT array */
extern void add_iFFTc2r_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, double *iFFT_c2r);

