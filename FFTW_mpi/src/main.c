#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <complex.h>
#include <time.h>

#include <mpi.h>
#include "hdf5.h"

#include "HDF5_utils.h"
#include "testfunctions.h"
#include "mpi_utils.h"


int procID;


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


	RunOneDimensionalTests(grp_1D_id, Nx, xmin, dx);



	/* Free memory */
	free(x_arr_local);
	
	/* Close HDF5 info */
    status = H5Gclose(grp_1D_id);
	status = H5Fclose(file_id);

	gettimeofday(&t_end, NULL);

    time_elapsed_us = (t_end.tv_sec - t_start.tv_sec) * 1.e6;
    time_elapsed_us += t_end.tv_usec - t_start.tv_usec;
    printf("--- Rank %d: Total elapsed time : %.6f secs \n", procID, time_elapsed_us * 1e-6);

	MPI_Finalize();

	return 0;
}

