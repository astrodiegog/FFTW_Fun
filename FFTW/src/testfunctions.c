#include "testfunctions.h"

extern void TestFunctionOne(int Nx, double *x_arr, fftw_complex *y_arr, double a)
{
	int i;
	double x;
	for (i = 0; i < Nx; i++)
	{
		x = x_arr[i];
		if ( fabs(x) > (1. / (2. * a)) )
		{
			y_arr[i] = 0. + 0. * I;
		}
		else
		{
			y_arr[i] = 1. + 0. * I;
		}
	} 

	return;
}


extern void TestFunctionOne_FFT(int Nx, double *kx_arr, fftw_complex *FFT_analytic, double a)
{
    int i;
    double kx;
    for (i = 0; i < Nx; i++)
    {
        kx = kx_arr[i];
		if (kx == 0)
		{
			FFT_analytic[i] = 1. + 0. * I;
		}
		else
		{
			FFT_analytic[i] = sin((M_PI * kx) / a) / ((kx * M_PI) / a)  + 0. * I;
		}
    }

    return;
}


extern void TestFunctionTwo(int Nx, double *x_arr, fftw_complex *y_arr, double a)
{   
    int i; 
    double x, x2;
    for (i = 0; i < Nx; i++)
    {   
        x = x_arr[i];
		x2 = x * x;
		y_arr[i] = exp(-1. * a * x2);
    }
    
    return;
}


extern void TestFunctionTwo_FFT(int Nx, double *kx_arr, fftw_complex *FFT_analytic, double a)
{   
    int i; 
    double kx, kx2, exp_arg;
    for (i = 0; i < Nx; i++)
    {   
        kx = kx_arr[i];
		kx2 = kx * kx;
		exp_arg = -1. * M_PI * M_PI * kx2 / a;
		FFT_analytic[i] = sqrt(M_PI / a) * exp(exp_arg);

    }
    
    return;
}

extern void TestFunctionThree(int Nx, double *x_arr, fftw_complex *y_arr, double a)
{
    int i;
    double x;
    for (i = 0; i < Nx; i++)
    {
        x = x_arr[i];
        y_arr[i] = exp(-1. * a * fabs(x));
    }

    return;
}


extern void TestFunctionThree_FFT(int Nx, double *kx_arr, fftw_complex *FFT_analytic, double a)
{
    int i;
    double kx, kx2, num, den;
    for (i = 0; i < Nx; i++)
    {
        kx = kx_arr[i];
        kx2 = kx * kx;
		num = 2. * a;
        den = (a * a) + (4. * M_PI * M_PI * kx2);
        FFT_analytic[i] = num / den;
    }

    return;
}



extern void RunTestOne(hid_t grp_test_id, double *x_arr, hid_t dataspace_id_c, double *kx_arr_c, hid_t dataspace_id_r, double *kx_arr_r, int Nx, int Nx_r, double a)
{
	/* Declare variables */
	fftw_complex *y_arr;
    fftw_complex *FFT_c2c, *iFFT_c2c, *FFT_analytic_c2c;
    fftw_complex *FFT_r2c, *FFT_analytic_r2c;
    double *iFFT_c2r;

    double y_arr_Real[Nx];

    fftw_plan plan_FFT_c2c, plan_FFT_r2c;
    fftw_plan plan_iFFT_c2c, plan_iFFT_c2r;

	int i;

	/* Allocate memory for arrays+plan */
    y_arr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx_r);
    FFT_analytic_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_analytic_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx_r);
	iFFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    iFFT_c2r = (double *) malloc(sizeof(double) * Nx);

	/* Create FFT c2c plans */
    plan_FFT_c2c = fftw_plan_dft_1d(Nx, y_arr, FFT_c2c, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_iFFT_c2c = fftw_plan_dft_1d(Nx, FFT_c2c, iFFT_c2c, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Create FFT r2c and c2r plans */
    plan_FFT_r2c = fftw_plan_dft_r2c_1d(Nx, y_arr_Real, FFT_r2c, FFTW_ESTIMATE);
    plan_iFFT_c2r = fftw_plan_dft_c2r_1d(Nx, FFT_r2c, iFFT_c2r, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

    /* Evaluate f(x) , grab real values for r2c */
    TestFunctionOne(Nx, &x_arr[0], &y_arr[0], a);
    for (i = 0; i < Nx; i++)
    {
        y_arr_Real[i] = creal(y_arr[i]);
    }

    /* Evalute F(k) */
    TestFunctionOne_FFT(Nx, &kx_arr_c[0], &FFT_analytic_c2c[0], a);
    TestFunctionOne_FFT(Nx_r, &kx_arr_r[0], &FFT_analytic_r2c[0], a);

    /* Execute FFT */
    fftw_execute(plan_FFT_c2c);
    fftw_execute(plan_FFT_r2c);
    fftw_execute(plan_iFFT_c2c);
    fftw_execute(plan_iFFT_c2r);

	/* Write info: f(x), FFT_c2c, iFFT_c2c, FFT_analytic_c2c, FFT_r2c, iFFT_r2c, FFT_analytic_r2c */
	Write_FFTWarr_1Dgrouptest(grp_test_id, "y_arr", dataspace_id_c, &y_arr[0], Nx);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_c2c", dataspace_id_c, &FFT_c2c[0], Nx);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_r2c", dataspace_id_r, &FFT_r2c[0], Nx_r);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace_id_c, &FFT_analytic_c2c[0], Nx);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_r2c", dataspace_id_r, &FFT_analytic_r2c[0], Nx_r);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "iFFT_c2c", dataspace_id_c, &iFFT_c2c[0], Nx);
    Write_HDF5_1Dgrouptest(grp_test_id, "iFFT_c2r", dataspace_id_c, &iFFT_c2r[0]);

	/* Destroy FFT plans */
    fftw_destroy_plan(plan_FFT_c2c);
    fftw_destroy_plan(plan_FFT_r2c);
    fftw_destroy_plan(plan_iFFT_c2c);
    fftw_destroy_plan(plan_iFFT_c2r);

	/* Free memory */
    fftw_free(y_arr);
    fftw_free(FFT_c2c);
    fftw_free(FFT_analytic_c2c);
    fftw_free(iFFT_c2c);
    fftw_free(FFT_r2c);
    fftw_free(FFT_analytic_r2c);
    free(iFFT_c2r);
}


extern void RunTestTwo(hid_t grp_test_id, double *x_arr, hid_t dataspace_id_c, double *kx_arr_c, hid_t dataspace_id_r, double *kx_arr_r, int Nx, int Nx_r, double a)
{
    /* Declare variables */
    fftw_complex *y_arr;
    fftw_complex *FFT_c2c, *iFFT_c2c, *FFT_analytic_c2c;
    fftw_complex *FFT_r2c, *FFT_analytic_r2c;
    double *iFFT_c2r;

    double y_arr_Real[Nx];

    fftw_plan plan_FFT_c2c, plan_FFT_r2c;
    fftw_plan plan_iFFT_c2c, plan_iFFT_c2r;

    int i;

    /* Allocate memory for arrays+plan */
    y_arr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx_r);
    FFT_analytic_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_analytic_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx_r);
    iFFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    iFFT_c2r = (double *) malloc(sizeof(double) * Nx);

    /* Create FFT c2c plans */
    plan_FFT_c2c = fftw_plan_dft_1d(Nx, y_arr, FFT_c2c, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_iFFT_c2c = fftw_plan_dft_1d(Nx, FFT_c2c, iFFT_c2c, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Create FFT r2c and c2r plans */
    plan_FFT_r2c = fftw_plan_dft_r2c_1d(Nx, y_arr_Real, FFT_r2c, FFTW_ESTIMATE);
    plan_iFFT_c2r = fftw_plan_dft_c2r_1d(Nx, FFT_r2c, iFFT_c2r, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

    /* Evaluate f(x) , grab real values for r2c */
    TestFunctionTwo(Nx, &x_arr[0], &y_arr[0], a);
    for (i = 0; i < Nx; i++)
    {
        y_arr_Real[i] = creal(y_arr[i]);
    }

    /* Evalute F(k) */
    TestFunctionTwo_FFT(Nx, &kx_arr_c[0], &FFT_analytic_c2c[0], a);
    TestFunctionTwo_FFT(Nx_r, &kx_arr_r[0], &FFT_analytic_r2c[0], a);

    /* Execute FFT */
    fftw_execute(plan_FFT_c2c);
    fftw_execute(plan_FFT_r2c);
    fftw_execute(plan_iFFT_c2c);
    fftw_execute(plan_iFFT_c2r);

    /* Write info: f(x), FFT_c2c, iFFT_c2c, FFT_analytic_c2c, FFT_r2c, iFFT_r2c, FFT_analytic_r2c */
    Write_FFTWarr_1Dgrouptest(grp_test_id, "y_arr", dataspace_id_c, &y_arr[0], Nx);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_c2c", dataspace_id_c, &FFT_c2c[0], Nx);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_r2c", dataspace_id_r, &FFT_r2c[0], Nx_r);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace_id_c, &FFT_analytic_c2c[0], Nx);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_r2c", dataspace_id_r, &FFT_analytic_r2c[0], Nx_r);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "iFFT_c2c", dataspace_id_c, &iFFT_c2c[0], Nx);
    Write_HDF5_1Dgrouptest(grp_test_id, "iFFT_c2r", dataspace_id_c, &iFFT_c2r[0]);

    /* Destroy FFT plans */
    fftw_destroy_plan(plan_FFT_c2c);
    fftw_destroy_plan(plan_FFT_r2c);
    fftw_destroy_plan(plan_iFFT_c2c);
    fftw_destroy_plan(plan_iFFT_c2r);

    /* Free memory */
    fftw_free(y_arr);
    fftw_free(FFT_c2c);
    fftw_free(FFT_analytic_c2c);
    fftw_free(iFFT_c2c);
    fftw_free(FFT_r2c);
    fftw_free(FFT_analytic_r2c);
    free(iFFT_c2r);
}


extern void RunTestThree(hid_t grp_test_id, double *x_arr, hid_t dataspace_id_c, double *kx_arr_c, hid_t dataspace_id_r, double *kx_arr_r, int Nx, int Nx_r, double a)
{
    /* Declare variables */
    fftw_complex *y_arr;
    fftw_complex *FFT_c2c, *iFFT_c2c, *FFT_analytic_c2c;
    fftw_complex *FFT_r2c, *FFT_analytic_r2c;
    double *iFFT_c2r;

    double y_arr_Real[Nx];

    fftw_plan plan_FFT_c2c, plan_FFT_r2c;
    fftw_plan plan_iFFT_c2c, plan_iFFT_c2r;

    int i;

    /* Allocate memory for arrays+plan */
    y_arr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx_r);
    FFT_analytic_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_analytic_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx_r);
    iFFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    iFFT_c2r = (double *) malloc(sizeof(double) * Nx);

    /* Create FFT c2c plans */
    plan_FFT_c2c = fftw_plan_dft_1d(Nx, y_arr, FFT_c2c, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_iFFT_c2c = fftw_plan_dft_1d(Nx, FFT_c2c, iFFT_c2c, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Create FFT r2c and c2r plans */
    plan_FFT_r2c = fftw_plan_dft_r2c_1d(Nx, y_arr_Real, FFT_r2c, FFTW_ESTIMATE);
    plan_iFFT_c2r = fftw_plan_dft_c2r_1d(Nx, FFT_r2c, iFFT_c2r, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

    /* Evaluate f(x) , grab real values for r2c */
    TestFunctionThree(Nx, &x_arr[0], &y_arr[0], a);
    for (i = 0; i < Nx; i++)
    {
        y_arr_Real[i] = creal(y_arr[i]);
    }

    /* Evalute F(k) */
    TestFunctionThree_FFT(Nx, &kx_arr_c[0], &FFT_analytic_c2c[0], a);
    TestFunctionThree_FFT(Nx_r, &kx_arr_r[0], &FFT_analytic_r2c[0], a);

    /* Execute FFT */
    fftw_execute(plan_FFT_c2c);
    fftw_execute(plan_FFT_r2c);
    fftw_execute(plan_iFFT_c2c);
    fftw_execute(plan_iFFT_c2r);

    /* Write info: f(x), FFT_c2c, iFFT_c2c, FFT_analytic_c2c, FFT_r2c, iFFT_r2c, FFT_analytic_r2c */
    Write_FFTWarr_1Dgrouptest(grp_test_id, "y_arr", dataspace_id_c, &y_arr[0], Nx);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_c2c", dataspace_id_c, &FFT_c2c[0], Nx);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_r2c", dataspace_id_r, &FFT_r2c[0], Nx_r);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_c2c", dataspace_id_c, &FFT_analytic_c2c[0], Nx);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "FFT_analytic_r2c", dataspace_id_r, &FFT_analytic_r2c[0], Nx_r);
    Write_FFTWarr_1Dgrouptest(grp_test_id, "iFFT_c2c", dataspace_id_c, &iFFT_c2c[0], Nx);
    Write_HDF5_1Dgrouptest(grp_test_id, "iFFT_c2r", dataspace_id_c, &iFFT_c2r[0]);

    /* Destroy FFT plans */
    fftw_destroy_plan(plan_FFT_c2c);
    fftw_destroy_plan(plan_FFT_r2c);
    fftw_destroy_plan(plan_iFFT_c2c);
    fftw_destroy_plan(plan_iFFT_c2r);

    /* Free memory */
    fftw_free(y_arr);
    fftw_free(FFT_c2c);
    fftw_free(FFT_analytic_c2c);
    fftw_free(iFFT_c2c);
    fftw_free(FFT_r2c);
    fftw_free(FFT_analytic_r2c);
	free(iFFT_c2r);
}
