#include "testfunctions.h"



extern void TestFunctionFive(int Nx, int Ny, double *x_arr, double *y_arr, fftw_complex *fxy_arr)
{
    int i, j, indx;
    double x, y;
    for (i = 0; i < Nx; i++)
    {   
        x = x_arr[i];
        for (j = 0; j < Ny; j++)
        {   
            y = y_arr[j];
            indx = j + Ny * i; 
            if (x == 0 && y == 0)
            {   
                /* Avoid dividing by zero */
                fxy_arr[indx] = 0. + 0. * I;
            }
            else
            {   
                fxy_arr[indx] = (1. * I )/ (x + y * I);
            }
        }
    }
    return;
}

extern void TestFunctionFive_FFT(int Nx, int Ny, double *kx_arr, double *ky_arr, fftw_complex *FFT_analytic)
{
    int i, j, indx;
    double kx, ky;
    for (i = 0; i < Nx; i++)
    {
        kx = kx_arr[i];
        for (j = 0; j < Ny; j++)
        {
            ky = ky_arr[j];
            indx = j + Ny * i;
            if (kx == 0 && ky == 0)
            {
                /* Avoid dividing by zero */
                FFT_analytic[indx] = 0. + 0. * I;
            }
            else
            {
                FFT_analytic[indx] = 1. / (kx + ky * I) ;
            }
        }
    }

    return;
}


extern void TestFunctionSix(int Nx, int Ny, double *x_arr, double *y_arr, fftw_complex *fxy_arr, double a, double b)
{   
    int i, j, indx;
    double x, x2, y, y2, a2, b2, exp_arg;
	a2 = a * a;
	b2 = b * b;
    for (i = 0; i < Nx; i++)
    {   
        x = x_arr[i];
		x2 = x * x;
        for (j = 0; j < Ny; j++)
        {   
            y = y_arr[j]; 
			y2 = y * y;
            indx = j + Ny * i; 
   			exp_arg = -1. * M_PI * (a2 * x2 + b2 * y2);
            fxy_arr[indx] = exp(exp_arg);
        }
    }
    return;
}


extern void TestFunctionSix_FFT(int Nx, int Ny, double *kx_arr, double *ky_arr, fftw_complex *FFT_analytic, double a, double b)
{
    int i, j, indx;
    double kx, kx2, ky, ky2, a2, b2, exp_arg;
	a2 = a * a;
	b2 = b * b;
    for (i = 0; i < Nx; i++)
    {
        kx = kx_arr[i];
		kx2 = kx * kx;
        for (j = 0; j < Ny; j++)
        {
            ky = ky_arr[j];
			ky2 = ky * ky;
            indx = j + Ny * i;
			exp_arg = -1. * M_PI * (kx2 / a2 + ky2 / b2);
			FFT_analytic[indx] = (1. / fabs( a * b)) * exp(exp_arg);

        }       
    }       
        
    return;
}




extern void RunTestFive(hid_t grp_test_id, double *x_arr, double *y_arr, hid_t dataspace_id_c, double *kx_arr_c, double *ky_arr_c, hid_t dataspace_id_r, double *ky_arr_r, int Nx, int Ny, int Ny_r)
{
    /* Declare variables */
    fftw_complex *fxy_arr;
    fftw_complex *FFT_c2c, *iFFT_c2c, *FFT_analytic_c2c;
    fftw_complex *FFT_r2c, *FFT_analytic_r2c;
    double *iFFT_c2r;

	double *fxy_arr_Real;

    fftw_plan plan_FFT_c2c, plan_FFT_r2c;
    fftw_plan plan_iFFT_c2c, plan_iFFT_c2r;

    int i, j, indx;


    /* Allocate memory for arrays+plan */
    fxy_arr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
	fxy_arr_Real = (double *) fftw_malloc(sizeof(double) * Nx * Ny);
    FFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    FFT_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny_r);
    FFT_analytic_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    FFT_analytic_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny_r);
    iFFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    iFFT_c2r = (double *) fftw_malloc(sizeof(double) * Nx * Ny);

    /* Create FFT c2c plans */
    plan_FFT_c2c = fftw_plan_dft_2d(Nx, Ny, fxy_arr, FFT_c2c, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_iFFT_c2c = fftw_plan_dft_2d(Nx, Ny, FFT_c2c, iFFT_c2c, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Create FFT r2c and c2r plans */
    plan_FFT_r2c = fftw_plan_dft_r2c_2d(Nx, Ny, fxy_arr_Real, FFT_r2c, FFTW_ESTIMATE);
    plan_iFFT_c2r = fftw_plan_dft_c2r_2d(Nx, Ny, FFT_r2c, iFFT_c2r, FFTW_ESTIMATE );

    /* Evaluate f(x) , grab real values for r2c */
    TestFunctionFive(Nx, Ny, &x_arr[0], &y_arr[0], &fxy_arr[0]);
    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            indx = j + Ny * i;
            fxy_arr_Real[indx] = creal(fxy_arr[indx]);
        }
    }

    /* Evalute F(k) */
    TestFunctionFive_FFT(Nx, Ny, &kx_arr_c[0], &ky_arr_c[0], &FFT_analytic_c2c[0]);
    TestFunctionFive_FFT(Nx, Ny_r, &kx_arr_c[0], &ky_arr_r[0], &FFT_analytic_r2c[0]);

    /* Execute FFT */
    fftw_execute(plan_FFT_c2c);
    fftw_execute(plan_FFT_r2c);

    /* Write f(x), FFT_c2c, FFT_analytic_c2c, FFT_r2c, FFT_analytic_r2c as iFFT destroys input */
    Write_FFTWarr_2Dgrouptest(grp_test_id, "fxy_arr", dataspace_id_c, &fxy_arr[0], Nx, Ny);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_c2c", dataspace_id_c, &FFT_c2c[0], Nx, Ny);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_r2c", dataspace_id_r, &FFT_r2c[0], Nx, Ny_r);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace_id_c, &FFT_analytic_c2c[0], Nx, Ny);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_analytic_r2c", dataspace_id_r, &FFT_analytic_r2c[0], Nx, Ny_r);

	/* Execute iFFT */
    fftw_execute(plan_iFFT_c2c);
    fftw_execute(plan_iFFT_c2r);

    /* Write iFFT_c2c, iFFT_r2c */
    Write_FFTWarr_2Dgrouptest(grp_test_id, "iFFT_c2c", dataspace_id_c, &iFFT_c2c[0], Nx, Ny);
    Write_HDF5_dataset(grp_test_id, "iFFT_c2r", dataspace_id_c, &iFFT_c2r[0]);


    /* Destroy FFT plans */
    fftw_destroy_plan(plan_FFT_c2c);
    fftw_destroy_plan(plan_FFT_r2c);
    fftw_destroy_plan(plan_iFFT_c2c);
    fftw_destroy_plan(plan_iFFT_c2r);

    /* Free memory */
    fftw_free(fxy_arr);
    fftw_free(FFT_c2c);
    fftw_free(FFT_analytic_c2c);
    fftw_free(iFFT_c2c);
    fftw_free(FFT_r2c);
    fftw_free(FFT_analytic_r2c);
    fftw_free(iFFT_c2r);
}

extern void RunTestSix(hid_t grp_test_id, double *x_arr, double *y_arr, hid_t dataspace_id_c, double *kx_arr_c, double *ky_arr_c, hid_t dataspace_id_r, double *ky_arr_r, int Nx, int Ny, int Ny_r, double a, double b)
{
    /* Declare variables */
    fftw_complex *fxy_arr;
    fftw_complex *FFT_c2c, *iFFT_c2c, *FFT_analytic_c2c;
    fftw_complex *FFT_r2c, *FFT_analytic_r2c;
    double *iFFT_c2r;

    double *fxy_arr_Real;

    fftw_plan plan_FFT_c2c, plan_FFT_r2c;
    fftw_plan plan_iFFT_c2c, plan_iFFT_c2r;

    int i, j, indx;


    /* Allocate memory for arrays+plan */
    fxy_arr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    fxy_arr_Real = (double *) fftw_malloc(sizeof(double) * Nx * Ny);
    FFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    FFT_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny_r);
    FFT_analytic_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    FFT_analytic_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny_r);
    iFFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    iFFT_c2r = (double *) fftw_malloc(sizeof(double) * Nx * Ny);

    /* Create FFT c2c plans */
    plan_FFT_c2c = fftw_plan_dft_2d(Nx, Ny, fxy_arr, FFT_c2c, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_iFFT_c2c = fftw_plan_dft_2d(Nx, Ny, FFT_c2c, iFFT_c2c, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Create FFT r2c and c2r plans */
    plan_FFT_r2c = fftw_plan_dft_r2c_2d(Nx, Ny, fxy_arr_Real, FFT_r2c, FFTW_ESTIMATE);
    plan_iFFT_c2r = fftw_plan_dft_c2r_2d(Nx, Ny, FFT_r2c, iFFT_c2r, FFTW_ESTIMATE );

    /* Evaluate f(x) , grab real values for r2c */
    TestFunctionSix(Nx, Ny, &x_arr[0], &y_arr[0], &fxy_arr[0], a, b);
    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            indx = j + Ny * i;
            fxy_arr_Real[indx] = creal(fxy_arr[indx]);
        }
    }

    /* Evalute F(k) */
    TestFunctionSix_FFT(Nx, Ny, &kx_arr_c[0], &ky_arr_c[0], &FFT_analytic_c2c[0], a, b);
    TestFunctionSix_FFT(Nx, Ny_r, &kx_arr_c[0], &ky_arr_r[0], &FFT_analytic_r2c[0], a, b);

    /* Execute FFT */
    fftw_execute(plan_FFT_c2c);
    fftw_execute(plan_FFT_r2c);

    /* Write f(x), FFT_c2c, FFT_analytic_c2c, FFT_r2c, FFT_analytic_r2c as iFFT destroys input */
    Write_FFTWarr_2Dgrouptest(grp_test_id, "fxy_arr", dataspace_id_c, &fxy_arr[0], Nx, Ny);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_c2c", dataspace_id_c, &FFT_c2c[0], Nx, Ny);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_r2c", dataspace_id_r, &FFT_r2c[0], Nx, Ny_r);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace_id_c, &FFT_analytic_c2c[0], Nx, Ny);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_analytic_r2c", dataspace_id_r, &FFT_analytic_r2c[0], Nx, Ny_r);

    /* Execute iFFT */
    fftw_execute(plan_iFFT_c2c);
    fftw_execute(plan_iFFT_c2r);

	/* Write iFFT_c2c, iFFT_r2c */
    Write_FFTWarr_2Dgrouptest(grp_test_id, "iFFT_c2c", dataspace_id_c, &iFFT_c2c[0], Nx, Ny);
    Write_HDF5_dataset(grp_test_id, "iFFT_c2r", dataspace_id_c, &iFFT_c2r[0]);


    /* Destroy FFT plans */
    fftw_destroy_plan(plan_FFT_c2c);
    fftw_destroy_plan(plan_FFT_r2c);
    fftw_destroy_plan(plan_iFFT_c2c);
    fftw_destroy_plan(plan_iFFT_c2r);

    /* Free memory */
    fftw_free(fxy_arr);
    fftw_free(FFT_c2c);
    fftw_free(FFT_analytic_c2c);
    fftw_free(iFFT_c2c);
    fftw_free(FFT_r2c);
    fftw_free(FFT_analytic_r2c);
    fftw_free(iFFT_c2r);
}



fftw_complex TestFunctionOne(double x, double a)
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

fftw_complex TestFunctionTwo(double x, double a)
{
    double x2;
    x2 = x * x;

    return exp(-1. * a * x2) + 0. * I;
}


fftw_complex TestFunctionOne_FFT(double kx, double a)
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


fftw_complex TestFunctionTwo_FFT(double kx, double a)
{
    double kx2, exp_arg;
    kx2 = kx * kx;
    exp_arg = -1. * M_PI * M_PI * kx2 / a;
    return sqrt(M_PI / a) * exp(exp_arg);
}

fftw_complex TestFunctionThree(double x, double a)
{   
    return exp(-1. * a * fabs(x));
}


fftw_complex TestFunctionThree_FFT(double kx, double a)
{   
    double kx2, num, den;
	kx2 = kx * kx;
	num = 2. * a; 
	den = (a * a) + (4. * M_PI * M_PI * kx2);
    
    return num / den;
}



extern void RunOneDimensionalTests(hid_t grp_1D_id, int Nx, double xmin, double dx)
{
	/* Declare info for HDF5 file */
	hid_t grp_test_id;
	herr_t status;

	/* Declare an iterator */
	int i;

	/* Declare testing info */
	double a;

	/* Declare array of dimensions for datasets */
    hsize_t dims1D_c[1];
    int Rank = 1;

	/* Declare fftw info*/
    ptrdiff_t N0;
    fftw_plan plan_FFT_c2c, plan_iFFT_c2c;
    ptrdiff_t alloc_local_FFT, local_ni_FFT, local_i_start_FFT, local_no_FFT, local_o_start_FFT;
    ptrdiff_t alloc_local_iFFT, local_ni_iFFT, local_i_start_iFFT, local_no_iFFT, local_o_start_iFFT;


	/* Declare all FFT-related arrays*/
    double *x_arr_local_FFT, *kx_arr_local_FFT, *x_arr_local_iFFT;
    fftw_complex *fx_arr_local, *FFT_c2c_local, *FFT_analytic_c2c_local, *iFFT_c2c_local;
    hid_t dataspace1D_id_local_in_c_FFT, dataspace1D_id_local_out_c_FFT;
    hid_t dataspace1D_id_local_in_c_iFFT, dataspace1D_id_local_out_c_iFFT;

	/* Grab the amount of data allocated by local_size routines */
    N0 = Nx;

    alloc_local_FFT = fftw_mpi_local_size_1d(N0, MPI_COMM_WORLD,
                                    FFTW_FORWARD, FFTW_ESTIMATE,
                                    &local_ni_FFT, &local_i_start_FFT,
                                    &local_no_FFT, &local_o_start_FFT);

    alloc_local_iFFT = fftw_mpi_local_size_1d(N0, MPI_COMM_WORLD,
                                    FFTW_BACKWARD, FFTW_ESTIMATE,
                                    &local_ni_iFFT, &local_i_start_iFFT,
                                    &local_no_iFFT, &local_o_start_iFFT);

	/* Create in/out dataspaces for FFT/iFFT */
    dims1D_c[0] = local_ni_FFT;
    dataspace1D_id_local_in_c_FFT = H5Screate_simple(Rank, dims1D_c, NULL);
    dims1D_c[0] = local_no_FFT;
    dataspace1D_id_local_out_c_FFT = H5Screate_simple(Rank, dims1D_c, NULL);

    dims1D_c[0] = local_ni_iFFT;
    dataspace1D_id_local_in_c_iFFT = H5Screate_simple(Rank, dims1D_c, NULL);
    dims1D_c[0] = local_no_iFFT;
    dataspace1D_id_local_out_c_iFFT = H5Screate_simple(Rank, dims1D_c, NULL);

	/* Allocate memory */
	x_arr_local_FFT = (double *) fftw_malloc(sizeof(double) * alloc_local_FFT);
    kx_arr_local_FFT = (double *) fftw_malloc(sizeof(double) * alloc_local_FFT);
    x_arr_local_iFFT = (double *) fftw_malloc(sizeof(double) * alloc_local_iFFT);
    fx_arr_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_FFT);
    FFT_c2c_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_FFT);
    FFT_analytic_c2c_local = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * alloc_local_FFT);
    iFFT_c2c_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_iFFT);

	/* Create plans */
    plan_FFT_c2c = fftw_mpi_plan_dft_1d(N0, fx_arr_local, FFT_c2c_local, MPI_COMM_WORLD,
                                        FFTW_FORWARD, FFTW_ESTIMATE);
    plan_iFFT_c2c = fftw_mpi_plan_dft_1d(N0, FFT_c2c_local, iFFT_c2c_local, MPI_COMM_WORLD,
                                         FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Fill in x for FFT and iFFT */
    for (i = 0; i < local_ni_FFT; i++)
    {
        x_arr_local_FFT[i] = xmin + (local_i_start_FFT * dx) + i * dx;
    }

	for (i = 0; i < local_no_iFFT; i++)
    {
        x_arr_local_iFFT[i] = xmin + (local_o_start_iFFT * dx) + i * dx;
    }

	/* Fill in kx info */
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

	/* Write the xarray domain (input of FFT plan), kxarray domain (output of FFT plan), xarray domain (output of iFFT plan) */
	/* Expect both x_arr_local to be the same */
	Write_HDF5_dataset(grp_1D_id, "x_arr_local_FFT", dataspace1D_id_local_in_c_FFT, &x_arr_local_FFT[0]);
    Write_HDF5_dataset(grp_1D_id, "kx_arr_local_FFT", dataspace1D_id_local_in_c_FFT, &kx_arr_local_FFT[0]);
    Write_HDF5_dataset(grp_1D_id, "x_arr_local_iFFT", dataspace1D_id_local_in_c_iFFT, &x_arr_local_iFFT[0]);

    /* Run Test 1 */
    a = 2.;
    for (i = 0; i < local_ni_FFT; i++)
    {
        fx_arr_local[i] = TestFunctionOne(x_arr_local_FFT[i], a);
    }

    for (i = 0; i < local_no_FFT; i++)
    {
        FFT_analytic_c2c_local[i] = TestFunctionOne_FFT(kx_arr_local_FFT[i], a);
    }

	/* Create group for this test */
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
        fx_arr_local[i] = TestFunctionTwo(x_arr_local_FFT[i], a);
    }
    for (i = 0; i < local_no_FFT; i++)
    {
        FFT_analytic_c2c_local[i] = TestFunctionTwo_FFT(kx_arr_local_FFT[i], a);
    }

    grp_test_id = H5Gcreate(grp_1D_id, "TestTwo", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    Write_FFTWarr_1Dgrouptest(grp_test_id, "fx_arr", dataspace1D_id_local_in_c_FFT, &fx_arr_local[0], local_ni_FFT);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace1D_id_local_out_c_FFT, &FFT_analytic_c2c_local[0], local_no_FFT);
    fftw_execute(plan_FFT_c2c);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_c2c", dataspace1D_id_local_out_c_FFT, &FFT_c2c_local[0], local_no_FFT);
    fftw_execute(plan_iFFT_c2c);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "iFFT_c2c", dataspace1D_id_local_out_c_iFFT, &iFFT_c2c_local[0], local_no_iFFT);


	/* Run Test 3 */
    a = 17.;
    for (i = 0; i < local_ni_FFT; i++)
    {
        fx_arr_local[i] = TestFunctionThree(x_arr_local_FFT[i], a);
    }
    for (i = 0; i < local_no_FFT; i++)
    {
        FFT_analytic_c2c_local[i] = TestFunctionThree_FFT(kx_arr_local_FFT[i], a);
    }

    grp_test_id = H5Gcreate(grp_1D_id, "TestThree", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    Write_FFTWarr_1Dgrouptest(grp_test_id, "fx_arr", dataspace1D_id_local_in_c_FFT, &fx_arr_local[0], local_ni_FFT);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace1D_id_local_out_c_FFT, &FFT_analytic_c2c_local[0], local_no_FFT);
    fftw_execute(plan_FFT_c2c);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_c2c", dataspace1D_id_local_out_c_FFT, &FFT_c2c_local[0], local_no_FFT);
    fftw_execute(plan_iFFT_c2c);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "iFFT_c2c", dataspace1D_id_local_out_c_iFFT, &iFFT_c2c_local[0], local_no_iFFT);


	/* Destroy plans */
    fftw_destroy_plan(plan_FFT_c2c);
    fftw_destroy_plan(plan_iFFT_c2c);

    /* Free memory */
    free(x_arr_local_FFT);
    free(x_arr_local_iFFT);
    fftw_free(fx_arr_local);
    fftw_free(FFT_c2c_local);
    fftw_free(FFT_analytic_c2c_local);
    fftw_free(iFFT_c2c_local);

	/* Close testing group id */
	status = H5Gclose(grp_test_id);
}



fftw_complex TestFunctionFour(double x, double y)
{
    double x2, y2;

	if (x == 0 && y == 0)
	{
		/* Avoid dividing by zero */
		return 0. + 0. * I;
	}
	else
	{
		x2 = x * x;
		y2 = y * y;
		return (1. / sqrt(x2 + y2)) + 0. * I;
	}
}


fftw_complex TestFunctionFour_FFT(double kx, double ky)
{
    double kx2, ky2;

	if (kx == 0 && ky == 0)
	{
		/* Avoid dividing by zero */
		return 0. + 0. * I;
	}
	else
	{
		kx2 = kx * kx;
		ky2 = ky * ky;
		return (1. / sqrt(kx2 + ky2) ) + 0. * I;
	}
}


extern void RunTwoDimensionalTests(hid_t grp_2D_id, int Nx, int Ny, double xmin, double ymin, double dx, double dy)
{
    /* Declare info for HDF5 file */
    hid_t grp_test_id;
    herr_t status;

    /* Declare iterators */
    int i, j;

	/* Declare variable for indexing */
	int indx;

    /* Declare testing info */
    double a;

    /* Declare array of dimensions for datasets */
    hsize_t dims2D_c[2];
	hsize_t dims2D_r[2];
    int Rank = 2;

    /* Declare fftw info*/
    ptrdiff_t N0, N1;
    fftw_plan plan_FFT_c2c, plan_iFFT_c2c;
	ptrdiff_t alloc_local_c, local_n_c, local_n0_start_c;
	ptrdiff_t N1_r2c; // number of complex data for r2c transform
	fftw_plan plan_FFT_r2c, plan_iFFT_c2r;
	ptrdiff_t alloc_local_r, local_n_r, local_n0_start_r;


    /* Declare all FFT-related arrays*/
    double *x_arr_local, *kx_arr_local_c;
    double *y_arr_local, *ky_arr_local_c;
    fftw_complex *fxy_arr_local, *FFT_c2c_local, *FFT_analytic_c2c_local, *iFFT_c2c_local;
	double *fxy_arr_Real_local;
	fftw_complex *FFT_r2c_local, *FFT_analytic_r2c_local;
	double *iFFT_c2r;
    hid_t dataspace2D_id_local_c, dataspace2D_id_local_r;


    /* Grab the amount of data allocated by local_size routines */
    N0 = Nx;
	N1 = Ny;
	N1_r2c = N1/2 + 1;

	alloc_local_c = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
										   &local_n_c, &local_n0_start_c);

	alloc_local_r = fftw_mpi_local_size_2d(N0, N1_r2c, MPI_COMM_WORLD,
											&local_n_r, &local_n0_start_r);


	printf("--- Rank %d : I have %d cells from the global %d cells with %d offsets --- \n", 
			procID, local_n_c*N1, N0*N1, local_n0_start_c*N1);

	/* Create in/out dataspaces for FFT/iFFT */
    dims2D_c[0] = local_n_c;
	dims2D_c[1] = N1;
    dataspace2D_id_local_c = H5Screate_simple(Rank, dims2D_c, NULL);

	dims2D_r[0] = local_n_r;
	dims2D_r[1] = N1_r;
	dataspace2D_id_local_r = H5Screate_simple(Rank, dims2D_r, NULL);


    /* Allocate memory */
    x_arr_local = (double *) fftw_malloc(sizeof(double) * alloc_local_c);
    y_arr_local = (double *) fftw_malloc(sizeof(double) * alloc_local_c);
    kx_arr_local_c = (double *) fftw_malloc(sizeof(double) * alloc_local_c);
	ky_arr_local_c = (double *) fftw_malloc(sizeof(double) * alloc_local_c);
	

    fxy_arr_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_c);
    FFT_c2c_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_c);
    FFT_analytic_c2c_local = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * alloc_local_c);
    iFFT_c2c_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_c);

	fxy_arr_Real_local = (double *) fftw_malloc(sizeof(double) * 2 * alloc_local_r);
	FFT_r2c_local = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local_r);
    FFT_analytic_r2c_local = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * alloc_local_r);
    iFFT_c2r_local = (double *) fftw_malloc(sizeof(double) * 2 * alloc_local_r);

    /* Create plans */
	plan_FFT_c2c = fftw_mpi_plan_dft_2d(N0, N1, fxy_arr_local, FFT_c2c_local, MPI_COMM_WORLD,
										FFTW_FORWARD, FFTW_ESTIMATE);

    plan_iFFT_c2c = fftw_mpi_plan_dft_2d(N0, N1, FFT_c2c_local, iFFT_c2c_local, MPI_COMM_WORLD,
                                         FFTW_BACKWARD, FFTW_ESTIMATE);

	plan_FFT_r2c = fftw_mpi_plan_dft_r2c_2d(N0, N1, fxy_arr_Real_local, FFT_r2c_local, MPI_COMM_WORLD, FFTW_ESTIMATE);


    /* Fill in (x,y) */
    for (i = 0; i < local_n_c; i++)
    {
		for (j = 0; j < N1; j++)
		{
			indx = j + i * N1;
			x_arr_local[indx] = xmin + (local_n0_start_c * dx) + i * dx;
			y_arr_local[indx] = ymin + j * dy;
		}
    }


	/* Fill in (kx,ky) info */
    for (i = 0; i < local_n_c; i++)
    {
		for (j = 0; j < N1; j++)
		{
			indx = j + i * N1;

			/* Assigning kmodes assumes even number of local_n_c */
			if ( (int) (i + local_n0_start_c) > (int) ((N0 / 2) - 1) )
			{
				/* Negative frequencies */
				kx_arr_local_c[indx] = -( N0 - (i + local_n0_start_c)) / (dx * Nx);
			}
			else
			{
				/* Positive frequencies*/
				kx_arr_local_c[indx] = (i + local_n0_start_c) / (dx * Nx);
			}

			if ( j > (int) ((N1 / 2) - 1) )
            {
                /* Negative frequencies */
                ky_arr_local_c[indx] = -( N1 - (j) ) / (dy * Ny);
            }
            else
            {
                /* Positive frequencies*/
                ky_arr_local_c[indx] = j / (dy * Ny);
            }

        }
    }

	printf("--- --- Rank %d : I am responsible for domain (xmin,xmax) = (%f, %f) & (ymin, ymax = (%f, %f) --- \n",
           procID, x_arr_local[0], x_arr_local[local_n_c * N1 - 1],  
				   y_arr_local[0], y_arr_local[local_n_c * N1 - 1]);

	printf("--- --- Rank %d : I am responsible for k-domain (kxmin,kxmax) = (%f, %f) & (kymin, kymax = (%f, %f) --- \n",
           procID, kx_arr_local_c[0], kx_arr_local_c[local_n_c * N1 - 1],
                   ky_arr_local_c[0], ky_arr_local_c[local_n_c * N1 - 1]);

	 /* Write the x,y 2D arrays domain (input of FFT plan), kx,ky 2D arrays domain (output of FFT plan), x,y 2D arrays domain (output of iFFT plan) */
    /* Expect both x,y_arr_local to be the same */
    Write_HDF5_dataset(grp_2D_id, "x_arr_local", dataspace2D_id_local_c, &x_arr_local[0]);
	Write_HDF5_dataset(grp_2D_id, "y_arr_local", dataspace2D_id_local_c, &y_arr_local[0]);
    Write_HDF5_dataset(grp_2D_id, "kx_arr_local_c", dataspace2D_id_local_c, &kx_arr_local_c[0]);
	Write_HDF5_dataset(grp_2D_id, "ky_arr_local_c", dataspace2D_id_local_c, &ky_arr_local_c[0]);


	/* Run Test 4 */
    for (i = 0; i < local_n_c; i++)
    {
		for (j = 0; j < N1; j++)
        {
            indx = j + i * N1;
			fxy_arr_local[indx] = TestFunctionFour(x_arr_local[indx], y_arr_local[indx]);
        }
    }

	for (i = 0; i < local_n_c; i++)
    {
        for (j = 0; j < N1; j++)
        {
            indx = j + i * N1;
            FFT_analytic_c2c_local[indx] = TestFunctionFour_FFT(kx_arr_local_c[indx], ky_arr_local_c[indx]);
        }
    }


    /* Create group for this test */
    grp_test_id = H5Gcreate(grp_2D_id, "TestFour", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    Write_FFTWarr_2Dgrouptest(grp_test_id, "fxy_arr", dataspace2D_id_local_c, &fxy_arr_local[0], local_n_c, N1);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace2D_id_local_c, &FFT_analytic_c2c_local[0], local_n_c, N1);
    fftw_execute(plan_FFT_c2c);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "FFT_c2c", dataspace2D_id_local_c, &FFT_c2c_local[0], local_n_c, N1);
    fftw_execute(plan_iFFT_c2c);
    Write_FFTWarr_2Dgrouptest(grp_test_id, "iFFT_c2c", dataspace2D_id_local_c, &iFFT_c2c_local[0], local_n_c, N1);

	/* Destroy plans */
    fftw_destroy_plan(plan_FFT_c2c);
    fftw_destroy_plan(plan_iFFT_c2c);

    /* Free memory */
    free(x_arr_local);
	free(y_arr_local);
    free(kx_arr_local_c);
	free(ky_arr_local_c);
    fftw_free(fxy_arr_local);
    fftw_free(FFT_c2c_local);
    fftw_free(FFT_analytic_c2c_local);
    fftw_free(iFFT_c2c_local);

    /* Close testing group id */
    status = H5Gclose(grp_test_id);
}









