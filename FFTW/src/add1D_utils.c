#include "add1D_utils.h"


void Write_HDF5_1Dgrouptest(hid_t grp_test_id, char *arr_name, hid_t dataspace_id, double *data_arr)
{
	hid_t dataset_id;
    herr_t status;

    dataset_id = H5Dcreate(grp_test_id, arr_name, H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_arr);
    status = H5Dclose(dataset_id);

    return;
}


extern void add_yarr_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *y_arr, int Nx)
{
    int i;
    double y_arr_Real[Nx], y_arr_Imag[Nx];

    for (i = 0; i < Nx; i++)
    {
        y_arr_Real[i] = creal(y_arr[i]);
        y_arr_Imag[i] = cimag(y_arr[i]);
    }

	Write_HDF5_1Dgrouptest(grp_test_id, "y_arr_Real", dataspace_id, &y_arr_Real[0]);
	Write_HDF5_1Dgrouptest(grp_test_id, "y_arr_Imag", dataspace_id, &y_arr_Imag[0]);

    return;
}

extern void add_FFTc2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *FFT_c2c, int Nx)
{
    int i;
    double FFT_c2c_Real[Nx], FFT_c2c_Imag[Nx];

    for (i = 0; i < Nx; i++)
    {
        FFT_c2c_Real[i] = creal(FFT_c2c[i]);
        FFT_c2c_Imag[i] = cimag(FFT_c2c[i]);
    }

	Write_HDF5_1Dgrouptest(grp_test_id, "FFT_c2c_Real", dataspace_id, &FFT_c2c_Real[0]);
    Write_HDF5_1Dgrouptest(grp_test_id, "FFT_c2c_Imag", dataspace_id, &FFT_c2c_Imag[0]);

    return;
}

extern void add_FFTr2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *FFT_r2c, int Nx)
{
    int i;
    double FFT_r2c_Real[Nx], FFT_r2c_Imag[Nx];

    for (i = 0; i < Nx; i++)
    {
        FFT_r2c_Real[i] = creal(FFT_r2c[i]);
        FFT_r2c_Imag[i] = cimag(FFT_r2c[i]);
    }

	Write_HDF5_1Dgrouptest(grp_test_id, "FFT_r2c_Real", dataspace_id, &FFT_r2c_Real[0]);
    Write_HDF5_1Dgrouptest(grp_test_id, "FFT_r2c_Imag", dataspace_id, &FFT_r2c_Imag[0]);

    return;
}

extern void add_FFTanalyticc2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *FFT_analytic_c2c, int Nx)
{
    int i;
    double FFT_analytic_c2c_Real[Nx], FFT_analytic_c2c_Imag[Nx];

    for (i = 0; i < Nx; i++)
    {
        FFT_analytic_c2c_Real[i] = creal(FFT_analytic_c2c[i]);
        FFT_analytic_c2c_Imag[i] = cimag(FFT_analytic_c2c[i]);
    }

    Write_HDF5_1Dgrouptest(grp_test_id, "FFT_analytic_c2c_Real", dataspace_id, &FFT_analytic_c2c_Real[0]);
    Write_HDF5_1Dgrouptest(grp_test_id, "FFT_analytic_c2c_Imag", dataspace_id, &FFT_analytic_c2c_Imag[0]);

    return;
}



extern void add_FFTanalyticr2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *FFT_analytic_r2c, int Nx)
{
    int i;
    double FFT_analytic_r2c_Real[Nx], FFT_analytic_r2c_Imag[Nx];

    for (i = 0; i < Nx; i++)
    {
        FFT_analytic_r2c_Real[i] = creal(FFT_analytic_r2c[i]);
        FFT_analytic_r2c_Imag[i] = cimag(FFT_analytic_r2c[i]);
    }

    Write_HDF5_1Dgrouptest(grp_test_id, "FFT_analytic_r2c_Real", dataspace_id, &FFT_analytic_r2c_Real[0]);
    Write_HDF5_1Dgrouptest(grp_test_id, "FFT_analytic_r2c_Imag", dataspace_id, &FFT_analytic_r2c_Imag[0]);

    return;
}


extern void add_iFFTc2c_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, fftw_complex *iFFT_c2c, int Nx)
{
    int i;
    double iFFT_c2c_Real[Nx], iFFT_c2c_Imag[Nx];

    for (i = 0; i < Nx; i++)
    {
        iFFT_c2c_Real[i] = creal(iFFT_c2c[i]);
        iFFT_c2c_Imag[i] = cimag(iFFT_c2c[i]);
    }

	Write_HDF5_1Dgrouptest(grp_test_id, "iFFT_c2c_Real", dataspace_id, &iFFT_c2c_Real[0]);
    Write_HDF5_1Dgrouptest(grp_test_id, "iFFT_c2c_Imag", dataspace_id, &iFFT_c2c_Imag[0]);

    return;
}


extern void add_iFFTc2r_1Dgrouptest(hid_t grp_test_id, hid_t dataspace_id, double *iFFT_c2r)
{

    Write_HDF5_1Dgrouptest(grp_test_id, "iFFT_c2r", dataspace_id, &iFFT_c2r[0]);

    return;
}




