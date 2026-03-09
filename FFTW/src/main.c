#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <complex.h>

#include "hdf5.h"
#include <fftw3.h>

#include "testfunctions.h"



int main(int argc, char **argv)
{
	printf("waddup !\n");

#ifdef HOWDY
    printf(" --- HOWDY ! --- \n");
#endif //HOWDY

	if (argc > 1)
	{
		fprintf(stderr, "too much info! \n");
		return 0;
	}

	char *FileName = "FFTWFun_out.h5";
	hid_t file_id;
	hid_t grp_1D_id;
	hid_t grp_test_id;
	hid_t dataspace_id_c, dataspace_id_r;

	/* Create file */
	file_id = H5Fcreate(FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t dims_c[1], dims_r[1];
	int Nx = 2048;
	int Rank = 1; 

	/* Create x-arr */
	double xmin = -5.;
	double xmax = 5.;
	double x_arr[Nx];
	double dx = (xmax - xmin) / Nx;

	int i;
	for (i=0; i < Nx; i++)
	{
		x_arr[i] = xmin + i * dx;
	}

	/* Create kx-arr */
	double kx_arr_c[Nx];
	int Nx_r = (int) ((Nx / 2.) + 1.); // number of fft bins for real fft
	double kx_arr_r[Nx_r];

	/* Fill in positive frequencies */
	for (i=0; i < Nx_r; i++)
	{
		kx_arr_r[i] = i / (Nx * dx);
		kx_arr_c[i] = i / (Nx * dx);
	}

	/* Flip first negative frequency */
	if (Nx % 2)
	{
		kx_arr_c[Nx_r] *= -1.;
	}
	else
	{
		kx_arr_c[Nx_r - 1] *= -1.;
	}

	/* Fill in negative frequencies */
	for (i=Nx_r; i < Nx; i++)
	{
		kx_arr_c[i] = -1. * (Nx - i) / (Nx * dx);
	}


	/* Create dataspace */
	dims_c[0] = Nx;
	dims_r[0] = Nx_r;
	dataspace_id_c = H5Screate_simple(Rank, dims_c, NULL);
    dataspace_id_r = H5Screate_simple(Rank, dims_r, NULL);


	/* Create group */
	grp_1D_id = H5Gcreate(file_id, "/OneDimensionalTests", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/* Create+Write xarray & kxarrays attributes */
	hid_t attr_id;
	herr_t status;

	attr_id = H5Acreate(grp_1D_id, "x_arr", H5T_IEEE_F64BE, dataspace_id_c, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &x_arr);
	status = H5Aclose(attr_id);

	attr_id = H5Acreate(grp_1D_id, "kx_arr_c", H5T_IEEE_F64BE, dataspace_id_c, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &kx_arr_c);
    status = H5Aclose(attr_id);

	attr_id = H5Acreate(grp_1D_id, "kx_arr_r", H5T_IEEE_F64BE, dataspace_id_r, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &kx_arr_r);
    status = H5Aclose(attr_id);

	/* Test 1 */
	hid_t dataset_id;
	grp_test_id = H5Gcreate(grp_1D_id, "TestOne", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	double a = 2.;
	fftw_complex *y_arr, *FFT_c2c;
	fftw_plan plan;

	y_arr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    plan = fftw_plan_dft_1d(Nx, y_arr, FFT_c2c, FFTW_FORWARD, FFTW_ESTIMATE);
	

	TestFunctionOne(Nx, &x_arr[0], &y_arr[0], a);

	fftw_execute(plan);



	
	double y_arr_Real[Nx], y_arr_Imag[Nx];

	for (i = 0; i < Nx; i++)
	{
		y_arr_Real[i] = creal(y_arr[i]);
		y_arr_Imag[i] = cimag(y_arr[i]);
	}

	dataset_id = H5Dcreate(grp_test_id, "y_arr_Real", H5T_IEEE_F64BE, dataspace_id_c, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y_arr_Real);
	status = H5Dclose(dataset_id);

	dataset_id = H5Dcreate(grp_test_id, "y_arr_Imag", H5T_IEEE_F64BE, dataspace_id_c, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y_arr_Imag);
    status = H5Dclose(dataset_id);


	double FFT_c2c_Real[Nx], FFT_c2c_Imag[Nx];

	for (i = 0; i < Nx; i++)
	{
		FFT_c2c_Real[i] = creal(FFT_c2c[i]);
		FFT_c2c_Imag[i] = cimag(FFT_c2c[i]);
	}

	dataset_id = H5Dcreate(grp_test_id, "FFT_c2c_Real", H5T_IEEE_F64BE, dataspace_id_c, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &FFT_c2c_Real);
    status = H5Dclose(dataset_id);

    dataset_id = H5Dcreate(grp_test_id, "FFT_c2c_Imag", H5T_IEEE_F64BE, dataspace_id_c, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &FFT_c2c_Imag);
    status = H5Dclose(dataset_id);



	// what do I add for each test?
	// f(x)
	// FFT_c2c
	// iFFT_c2c
	// FFT_analytic_c2c
	// FFT_r2c
	// iFFT_r2c
	// FFT_analytic_r2c

	fftw_destroy_plan(plan);
    fftw_free(y_arr); 
	fftw_free(FFT_c2c);

	status = H5Gclose(grp_test_id);
	status = H5Fclose(file_id);

	return 0;
}

