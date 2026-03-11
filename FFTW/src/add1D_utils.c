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


extern void Write_FFTWarr_1Dgrouptest(hid_t grp_test_id, char *arr_prefix, hid_t dataspace_id, fftw_complex *FFTW_arr, int Nx)
{
	int i;
    double FFTW_arr_Real[Nx], FFTW_arr_Imag[Nx];
	char realarr_name[MAXLEN], imagarr_name[MAXLEN];

    for (i = 0; i < Nx; i++)
    {
        FFTW_arr_Real[i] = creal(FFTW_arr[i]);
        FFTW_arr_Imag[i] = cimag(FFTW_arr[i]);
    }

	strcpy(realarr_name, arr_prefix);
	strcpy(imagarr_name, arr_prefix);

	strcat(realarr_name, "_Real");
	strcat(imagarr_name, "_Imag");

	Write_HDF5_1Dgrouptest(grp_test_id, realarr_name, dataspace_id, &FFTW_arr_Real[0]);
    Write_HDF5_1Dgrouptest(grp_test_id, imagarr_name, dataspace_id, &FFTW_arr_Imag[0]);

    return;
}




