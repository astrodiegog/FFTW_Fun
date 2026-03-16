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




double *create_kxarr_c(int Nx, double dx)
{   
	int i;
    int Nx_r = (int) ((Nx / 2.) + 1.); // number of fft bins for real fft

	double *kx_arr_c = (double *) malloc(sizeof(double) * Nx);

    /* Fill in positive frequencies */
    for (i=0; i < Nx_r; i++)
    {
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
    
    return kx_arr_c;
}


double *create_kxarr_r(int Nx, double dx)
{   
    int i;
    int Nx_r = (int) ((Nx / 2.) + 1.); // number of fft bins for real fft
    
    double *kx_arr_r = (double *) malloc(sizeof(double) * Nx_r);

    /* Fill in positive frequencies */
    for (i=0; i < Nx_r; i++)
    {
        kx_arr_r[i] = i / (Nx * dx);
    }

	return kx_arr_r;

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

	/* Declare info for HDF5 file*/
	char *FileName = "FFTWFun_out.h5";
	hid_t file_id;
	hid_t grp_1D_id, grp_2D_id;
	hid_t grp_test_id;
	hid_t dataspace1D_id_c, dataspace1D_id_r;
	hid_t dataspace1D_x_id_c, dataspace1D_y_id_c;
	hid_t dataspace1D_x_id_r, dataspace1D_y_id_r;
	hid_t dataspace2D_id_c, dataspace2D_id_r;

	/* Declare array of dimensions */
	hsize_t dims1D_c[1], dims1D_r[1];
	hsize_t dims1D_x_c[1], dims1D_y_c[1];
	hsize_t dims1D_x_r[1], dims1D_y_r[1];
	hsize_t dims2D_c[2], dims2D_r[2];
	int Rank = 1;

	/* Declare info for x-arr and kx arrays (for complex and real FFT) */
	int Nx = 2048;
    int Nx_r = (int) ((Nx / 2.) + 1.); // number of fft bins for real fft
	double xmin = -5.;
	double xmax = 5.;
    double dx = (xmax - xmin) / Nx;
	double *x_arr, *kx_arr_c, *kx_arr_r;

	/* Declare same info for y-arrs, place data later on */
	int Ny, Ny_r ;
    double ymin, ymax, dy;
    double *y_arr, *ky_arr_c, *ky_arr_r;


	/* Create file */
    file_id = H5Fcreate(FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Create x-arr and kx arrays (for complex and real FFT) */
	x_arr = create_xarr(Nx, xmin, xmax);
	kx_arr_c = create_kxarr_c(Nx, dx);
    kx_arr_r = create_kxarr_r(Nx, dx);

	/* Create dataspace for real and complex FFT */
	dims1D_c[0] = Nx;
	dims1D_r[0] = Nx_r;
	dataspace1D_id_c = H5Screate_simple(Rank, dims1D_c, NULL);
    dataspace1D_id_r = H5Screate_simple(Rank, dims1D_r, NULL);

	/* Create group for 1D tests */
	grp_1D_id = H5Gcreate(file_id, "/OneDimensionalTests", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/* Create+Write xarray & kxarrays attributes */
	hid_t attr_id;
	herr_t status;

	attr_id = H5Acreate(grp_1D_id, "x_arr", H5T_IEEE_F64BE, dataspace1D_id_c, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, x_arr);
	status = H5Aclose(attr_id);

	attr_id = H5Acreate(grp_1D_id, "kx_arr_c", H5T_IEEE_F64BE, dataspace1D_id_c, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, kx_arr_c);
    status = H5Aclose(attr_id);

	attr_id = H5Acreate(grp_1D_id, "kx_arr_r", H5T_IEEE_F64BE, dataspace1D_id_r, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, kx_arr_r);
    status = H5Aclose(attr_id);

	/* Run Test 1 */
	grp_test_id = H5Gcreate(grp_1D_id, "TestOne", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	double a = 2.;
	RunTestOne(grp_test_id, &x_arr[0], dataspace1D_id_c, &kx_arr_c[0], dataspace1D_id_r, &kx_arr_r[0], Nx, Nx_r, a);


	/* Run Test 2 */
	grp_test_id = H5Gcreate(grp_1D_id, "TestTwo", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    a = 30.;
    RunTestTwo(grp_test_id, &x_arr[0], dataspace1D_id_c, &kx_arr_c[0], dataspace1D_id_r, &kx_arr_r[0], Nx, Nx_r, a);


	/* Run Test 3 */
    grp_test_id = H5Gcreate(grp_1D_id, "TestThree", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    a = 17.;
    RunTestThree(grp_test_id, &x_arr[0], dataspace1D_id_c, &kx_arr_c[0], dataspace1D_id_r, &kx_arr_r[0], Nx, Nx_r, a);
	free(x_arr);
	free(kx_arr_c);
    free(kx_arr_r);

	/* Define info for 2D tests */
	Nx = 256;
    Nx_r = (int) ((Nx / 2.) + 1.); // number of fft bins for real fft
	dx = (xmax - xmin) / Nx;

	Ny = 512;
    Ny_r = (int) ((Ny / 2.) + 1.); // number of fft bins for real fft
	ymin = -10.;
	ymax = 10.;
	dy = (ymax - ymin) / Ny;

	/* Create x and y arrays */
	x_arr = create_xarr(Nx, xmin, xmax);
    kx_arr_c = create_kxarr_c(Nx, dx);
    kx_arr_r = create_kxarr_r(Nx, dx);

	y_arr = create_xarr(Ny, ymin, ymax);
    ky_arr_c = create_kxarr_c(Ny, dy);
    ky_arr_r = create_kxarr_r(Ny, dy);
	
	/* Create dataspace for real and complex FFT */
	dims1D_x_c[0] = Nx;
	dims1D_y_c[0] = Ny;
	dataspace1D_x_id_c = H5Screate_simple(Rank, dims1D_x_c, NULL);
	dataspace1D_y_id_c = H5Screate_simple(Rank, dims1D_y_c, NULL);
	dims1D_x_r[0] = Nx_r;
    dims1D_y_r[0] = Ny_r;
	dataspace1D_x_id_r = H5Screate_simple(Rank, dims1D_x_r, NULL);
    dataspace1D_y_id_r = H5Screate_simple(Rank, dims1D_y_r, NULL);
	dims2D_c[0] = Nx;
    dims2D_c[1] = Ny;
    dims2D_r[0] = Nx;
    dims2D_r[1] = Ny_r;
	Rank = 2;
    dataspace2D_id_c = H5Screate_simple(Rank, dims2D_c, NULL);
    dataspace2D_id_r = H5Screate_simple(Rank, dims2D_r, NULL);

	/* Create group for 2D tests */
    grp_2D_id = H5Gcreate(file_id, "/TwoDimensionalTests", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    /* Create+Write xarray & kxarrays attributes */
	attr_id = H5Acreate(grp_2D_id, "x_arr", H5T_IEEE_F64BE, dataspace1D_x_id_c, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, x_arr);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(grp_2D_id, "kx_arr_c", H5T_IEEE_F64BE, dataspace1D_x_id_c, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, kx_arr_c);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(grp_2D_id, "kx_arr_r", H5T_IEEE_F64BE, dataspace1D_x_id_r, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, kx_arr_r);
    status = H5Aclose(attr_id);

	attr_id = H5Acreate(grp_2D_id, "y_arr", H5T_IEEE_F64BE, dataspace1D_y_id_c, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, y_arr);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(grp_2D_id, "ky_arr_c", H5T_IEEE_F64BE, dataspace1D_y_id_c, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, ky_arr_c);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(grp_2D_id, "ky_arr_r", H5T_IEEE_F64BE, dataspace1D_y_id_r, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, ky_arr_r);
    status = H5Aclose(attr_id);

	/* Run Test 4 */
    grp_test_id = H5Gcreate(grp_2D_id, "TestFour", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    RunTestFour(grp_test_id, &x_arr[0], &y_arr[0], dataspace2D_id_c, &kx_arr_c[0], &ky_arr_c[0], dataspace2D_id_r, &kx_arr_r[0], &ky_arr_r[0], Nx, Ny, Nx_r, Ny_r);



	/* Free malloc-ed arrays*/
	free(x_arr);
	free(kx_arr_c);
	free(kx_arr_r);

	free(y_arr);
    free(ky_arr_c);
    free(ky_arr_r);

	/* Close HDF5 objects */
	status = H5Gclose(grp_test_id);
	status = H5Gclose(grp_1D_id);
	status = H5Gclose(grp_2D_id);
	status = H5Fclose(file_id);

	return 0;
}

