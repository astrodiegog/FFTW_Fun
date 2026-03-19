#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <complex.h>
#include <time.h>

#include <mpi.h>
#include "hdf5.h"

#include "HDF5_utils.h"
#include "mpi_utils.h"


int procID;

fftw_complex TestFunctionOne_FR(double x, double a)
{
	if ( fabs(x) > (1. / (2. * a)) )
    {
		return  0. + 0. * I;
	}
	else
	{
		return 1. + 0. * I;
	}
}

fftw_complex TestFunctionTwo_FR(double x, double a)
{
    double x2;
	x2 = x * x;

    return exp(-1. * a * x2) + 0. * I;
}


fftw_complex TestFunctionOne_FFT_FR(double kx, double a)
{
	if (kx == 0)
	{
		return 1. + 0. * I;
	}
	else
	{
		return sin((M_PI * kx) / a) / ((kx * M_PI) / a)  + 0. * I;
	}
}


fftw_complex TestFunctionTwo_FFT_FR(double kx, double a)
{
    double kx2, exp_arg;
	kx2 = kx * kx;
	exp_arg = -1. * M_PI * M_PI * kx2 / a;
	return sqrt(M_PI / a) * exp(exp_arg);

}

int main(int argc, char **argv)
{
	/* Program info */
	int nprocs;
	MPI_Status status_mpi;
	
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
	char FileName[16] = "FFTWFun_out.h5.";
	
	char FileName_appendix[MAXLEN];
	hid_t file_id;
	hid_t grp_1D_id;
	hid_t grp_test_id;
	hid_t dataspace1D_id_c;
	hid_t attrs1D_id;

	herr_t status;

	/* Declare array of dimensions */
	hsize_t dims1D_c[1];
	hsize_t attrs1D[1], attrs2D[2];
	int Rank = 1;

	/* Declare info for x-arr */
	int Nx = 2048;
	double xmin = -5.;
	double xmax = 5.;
    double dx = (xmax - xmin) / Nx;
	double *x_arr_local, *kx_arr_c, *kx_arr_r;

	/* Declare x-arr info for processor */
	int np_x;
	int Nx_local, Nx_offset;

	/* Declare time info */
	struct timeval t_start, t_end;
	double time_elapsed_us;

	/* Call MPI routines & set nprocs and nprocID*/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	/* Start time tracking */
    gettimeofday(&t_start, NULL);

	printf("--- Rank %d : hello! --- \n", procID);

	/* fftw init */
	fftw_mpi_init();

	/* Assign number of processes to pointer of x */
	Tile_Decomposition1D(nprocs, &np_x);

	/* Use number of processes along x axis to place number of cells for each process*/
	Domain_Decomposition1D(Nx, np_x, &Nx_local);
	Nx_offset = Nx_local * procID;

	/* Create file */
	printf("--- Rank %d : I have %d cells from the global %d cells with %d offsets --- \n", procID, Nx_local, Nx, Nx_offset);
	sprintf(FileName_appendix, "%d", procID);
    strcat(FileName, FileName_appendix);
	printf("File name: %s \n", FileName);
    file_id = H5Fcreate(FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


	/* Create dataspace for real and complex FFT */
    dims1D_c[0] = Nx_local;
	attrs1D[0] = 1;
	int int_data[1];
	double double_data[1];
    dataspace1D_id_c = H5Screate_simple(Rank, dims1D_c, NULL);
	attrs1D_id = H5Screate_simple(Rank, attrs1D, NULL);

	/* Create group for 1D tests */
    grp_1D_id = H5Gcreate(file_id, "/OneDimensionalTests", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Create+Write xarray attributes */
    hid_t attr_id;

	x_arr_local = (double *) malloc(sizeof(double) * Nx_local);
	int i;
    for (i=0; i < Nx_local; i++)
    {
        x_arr_local[i] = xmin + (Nx_offset * dx) + i * dx;
    }

    Write_HDF5_dataset(grp_1D_id, "x_arr_local", dataspace1D_id_c, x_arr_local);

	int_data[0] = Nx;
	Write_HDF5_int_attribute(grp_1D_id, "dims_global", attrs1D_id, &int_data[0]);

	int_data[0] = Nx_local;
	Write_HDF5_int_attribute(grp_1D_id, "dims_local", attrs1D_id, &int_data[0]);	

	int_data[0] = Nx_offset;
    Write_HDF5_int_attribute(grp_1D_id, "dims_offset", attrs1D_id, &int_data[0]); 

	double_data[0] = dx;
	Write_HDF5_double_attribute(grp_1D_id, "dx", attrs1D_id, &double_data[0]);




	/* Declare fftw info*/
    ptrdiff_t N0;
    fftw_plan plan_FFT_c2c, plan_iFFT_c2c;
    ptrdiff_t alloc_local_FFT, local_ni_FFT, local_i_start_FFT, local_no_FFT, local_o_start_FFT;
	ptrdiff_t alloc_local_iFFT, local_ni_iFFT, local_i_start_iFFT, local_no_iFFT, local_o_start_iFFT;


	/* get local data size and allocate */
	double *x_arr_local_FFT, *kx_arr_local_FFT, *x_arr_local_iFFT;
	fftw_complex *fx_arr_local, *FFT_c2c_local, *FFT_analytic_c2c_local, *iFFT_c2c_local, *iFFT_analytic_c2c_local;
	hid_t dataspace1D_id_local_in_c_FFT, dataspace1D_id_local_out_c_FFT;
	hid_t dataspace1D_id_local_in_c_iFFT, dataspace1D_id_local_out_c_iFFT;
	N0 = Nx;

    alloc_local_FFT = fftw_mpi_local_size_1d(N0, MPI_COMM_WORLD,
									FFTW_FORWARD, FFTW_ESTIMATE,
									&local_ni_FFT, &local_i_start_FFT,
                                    &local_no_FFT, &local_o_start_FFT);

	alloc_local_iFFT = fftw_mpi_local_size_1d(N0, MPI_COMM_WORLD,
                                    FFTW_BACKWARD, FFTW_ESTIMATE,
                                    &local_ni_iFFT, &local_i_start_iFFT,
                                    &local_no_iFFT, &local_o_start_iFFT);

	dims1D_c[0] = local_ni_FFT;
    dataspace1D_id_local_in_c_FFT = H5Screate_simple(Rank, dims1D_c, NULL);
    dims1D_c[0] = local_no_FFT;
    dataspace1D_id_local_out_c_FFT = H5Screate_simple(Rank, dims1D_c, NULL);

	dims1D_c[0] = local_ni_iFFT;
    dataspace1D_id_local_in_c_iFFT = H5Screate_simple(Rank, dims1D_c, NULL);
    dims1D_c[0] = local_no_iFFT;
    dataspace1D_id_local_out_c_iFFT = H5Screate_simple(Rank, dims1D_c, NULL);


	x_arr_local_FFT = (double *) fftw_malloc(sizeof(double) * alloc_local_FFT);
	kx_arr_local_FFT = (double *) fftw_malloc(sizeof(double) * alloc_local_FFT);
	x_arr_local_iFFT = (double *) fftw_malloc(sizeof(double) * alloc_local_iFFT);
	fx_arr_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_FFT);
	FFT_c2c_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_FFT);
	FFT_analytic_c2c_local = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * alloc_local_FFT);
	iFFT_c2c_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_iFFT);

	plan_FFT_c2c = fftw_mpi_plan_dft_1d(N0, fx_arr_local, FFT_c2c_local, MPI_COMM_WORLD, 
										FFTW_FORWARD, FFTW_ESTIMATE);
	plan_iFFT_c2c = fftw_mpi_plan_dft_1d(N0, FFT_c2c_local, iFFT_c2c_local, MPI_COMM_WORLD, 
										 FFTW_BACKWARD, FFTW_ESTIMATE);

	/* fill in x and kx info for FFT and iFFT */
    for (i = 0; i < local_ni_FFT; i++)
	{
		x_arr_local_FFT[i] = xmin + (local_i_start_FFT * dx) + i * dx;
	}

	for (i = 0; i < local_no_FFT; i++)
	{
		/* Assigning kmodes assumes even number of local_ni */
        if ( (int) (i + local_i_start_FFT) > (int) ((N0 / 2) - 1) )
        {
            /* Negative frequencies */
            kx_arr_local_FFT[i] = -( N0 - (i + local_o_start_FFT)) / (dx * Nx);
        }
        else
        {
            /* Positive frequencies*/
            kx_arr_local_FFT[i] = (i + local_o_start_FFT) / (dx * Nx);
        }
	}

	for (i = 0; i < local_no_iFFT; i++)
    {
        x_arr_local_iFFT[i] = xmin + (local_o_start_iFFT * dx) + i * dx;
    }



	printf("--- Rank %d: I am allocating %d number of complex numbers \n", procID, alloc_local_FFT);
	printf("--- Rank %d: \t fx_local has %d cells with offset of %d \n", procID, local_ni_FFT, local_i_start_FFT);
    printf("--- Rank %d: \t FFT_local has %d cells with offset of %d \n", procID, local_no_FFT, local_o_start_FFT);
	printf("--- Rank %d: \t iFFT_local has %d cells with offset of %d \n", procID, local_no_iFFT, local_o_start_iFFT);

	Write_HDF5_dataset(grp_1D_id, "x_arr_local_FFT", dataspace1D_id_local_in_c_FFT, &x_arr_local_FFT[0]);
	Write_HDF5_dataset(grp_1D_id, "kx_arr_local_FFT", dataspace1D_id_local_in_c_FFT, &kx_arr_local_FFT[0]);
	Write_HDF5_dataset(grp_1D_id, "x_arr_local_iFFT", dataspace1D_id_local_in_c_iFFT, &x_arr_local_iFFT[0]);

	/* Run Test 1 */
	double a = 2.;
    /* initialize data to some function my_function(x,y) */
    for (i = 0; i < local_ni_FFT; i++)
    {
        fx_arr_local[i] = TestFunctionOne_FR(x_arr_local_FFT[i], a);
    }

	for (i = 0; i < local_no_FFT; i++)
    {
		FFT_analytic_c2c_local[i] = TestFunctionOne_FFT_FR(kx_arr_local_FFT[i], a);
	}


    grp_test_id = H5Gcreate(grp_1D_id, "TestOne", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	Write_FFTWarr_1Dgrouptest(grp_test_id, "fx_arr", dataspace1D_id_local_in_c_FFT, &fx_arr_local[0], local_ni_FFT);
	Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace1D_id_local_out_c_FFT, &FFT_analytic_c2c_local[0], local_no_FFT);
	fftw_execute(plan_FFT_c2c);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_c2c", dataspace1D_id_local_out_c_FFT, &FFT_c2c_local[0], local_no_FFT);
	fftw_execute(plan_iFFT_c2c);
	Write_FFTWarr_1Dgrouptest(grp_test_id, "iFFT_c2c", dataspace1D_id_local_out_c_iFFT, &iFFT_c2c_local[0], local_no_iFFT);


	/* Run Test 2 */
	a = 30.;
	for (i = 0; i < local_ni_FFT; i++)
    {
        fx_arr_local[i] = TestFunctionTwo_FR(x_arr_local_FFT[i], a);
		FFT_analytic_c2c_local[i] = TestFunctionTwo_FFT_FR(kx_arr_local_FFT[i], a);
    }
	for (i = 0; i < local_no_FFT; i++)
    {
        FFT_analytic_c2c_local[i] = TestFunctionTwo_FFT_FR(kx_arr_local_FFT[i], a);
    }

    grp_test_id = H5Gcreate(grp_1D_id, "TestTwo", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	Write_FFTWarr_1Dgrouptest(grp_test_id, "fx_arr", dataspace1D_id_local_in_c_FFT, &fx_arr_local[0], local_ni_FFT);
	Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace1D_id_local_out_c_FFT, &FFT_analytic_c2c_local[0], local_no_FFT);
    fftw_execute(plan_FFT_c2c);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_c2c", dataspace1D_id_local_out_c_FFT, &FFT_c2c_local[0], local_no_FFT);
    fftw_execute(plan_iFFT_c2c);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "iFFT_c2c", dataspace1D_id_local_out_c_iFFT, &iFFT_c2c_local[0], local_no_iFFT);
	

	// RunTestOne(grp_test_id, &x_arr[0], dataspace1D_id_c, &kx_arr_c[0], dataspace1D_id_r, &kx_arr_r[0], Nx, Nx_r, a);


	/* Destroy plan */
	fftw_destroy_plan(plan_FFT_c2c);
	fftw_destroy_plan(plan_iFFT_c2c);

	/* Free memory */
	free(x_arr_local);
	free(x_arr_local_FFT);
	free(x_arr_local_iFFT);
    fftw_free(fx_arr_local);
	fftw_free(FFT_c2c_local);
	fftw_free(FFT_analytic_c2c_local);
	fftw_free(iFFT_c2c_local);
	


	/* Close HDF5 objects */
    status = H5Gclose(grp_test_id);
    status = H5Gclose(grp_1D_id);
	status = H5Fclose(file_id);

	gettimeofday(&t_end, NULL);

    time_elapsed_us = (t_end.tv_sec - t_start.tv_sec) * 1.e6;
    time_elapsed_us += t_end.tv_usec - t_start.tv_usec;
    printf("--- Rank %d: Total elapsed time : %.6f secs \n", procID, time_elapsed_us * 1e-6);

	MPI_Finalize();

	return 0;
}

