#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <complex.h>

#include "hdf5.h"
#include <fftw3.h>

#include "testfunctions.h"



double *create_xarr(int Nx, double xmin, double xmax)
{
	double *x_arr = (double *) malloc(sizeof(double) * Nx);
	double dx = (xmax - xmin) / Nx;

	int i;
    for (i=0; i < Nx; i++)
    {
        x_arr[i] = xmin + i * dx;
    }

	return x_arr;
}






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
	int i; // iterator

	/* Create x-arr */
	double xmin = -5.;
	double xmax = 5.;
    double dx = (xmax - xmin) / Nx;
	double *x_arr;

	x_arr = create_xarr(Nx, xmin, xmax);

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


	/* Create 1D group */
	grp_1D_id = H5Gcreate(file_id, "/OneDimensionalTests", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/* Create+Write xarray & kxarrays attributes */
	hid_t attr_id;
	herr_t status;

	attr_id = H5Acreate(grp_1D_id, "x_arr", H5T_IEEE_F64BE, dataspace_id_c, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, x_arr);
	status = H5Aclose(attr_id);

	attr_id = H5Acreate(grp_1D_id, "kx_arr_c", H5T_IEEE_F64BE, dataspace_id_c, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &kx_arr_c);
    status = H5Aclose(attr_id);

	attr_id = H5Acreate(grp_1D_id, "kx_arr_r", H5T_IEEE_F64BE, dataspace_id_r, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &kx_arr_r);
    status = H5Aclose(attr_id);

	/* Run Test 1 */
	hid_t dataset_id;
	grp_test_id = H5Gcreate(grp_1D_id, "TestOne", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	double a = 2.;

	RunTestOne(grp_test_id, &x_arr[0], dataspace_id_c, &kx_arr_c[0], dataspace_id_r, &kx_arr_r[0], Nx, Nx_r, a);



	free(x_arr);
	status = H5Gclose(grp_test_id);
	status = H5Gclose(grp_1D_id);
	status = H5Fclose(file_id);

	return 0;
}

