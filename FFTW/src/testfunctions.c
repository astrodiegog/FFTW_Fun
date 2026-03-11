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




extern void RunTestOne(hid_t grp_test_id, double *x_arr, hid_t dataspace_id_c, double *kx_arr_c, hid_t dataspace_id_r, double *kx_arr_r, int Nx, int Nx_r, double a)
{
	/* Declare variables */
	fftw_complex *y_arr;
    fftw_complex *FFT_c2c, *iFFT_c2c, *FFT_analytic_c2c;
    fftw_complex *FFT_r2c, *FFT_analytic_r2c;
    double *iFFT_c2r;

    double y_arr_Real[Nx], y_arr_Imag[Nx];

    fftw_plan plan_FFT_c2c, plan_FFT_r2c;
    fftw_plan plan_iFFT_c2c, plan_iFFT_c2r;

	int i;

	/* Allocate memory for arrays+plan */
    y_arr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    FFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    iFFT_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    FFT_analytic_c2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    FFT_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx_r);
    iFFT_c2r = (double *) malloc(sizeof(double) * Nx_r);
    FFT_analytic_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx_r);

	/* Create FFT c2c plans */
    plan_FFT_c2c = fftw_plan_dft_1d(Nx, y_arr, FFT_c2c, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_iFFT_c2c = fftw_plan_dft_1d(Nx, FFT_c2c, iFFT_c2c, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Create FFT r2c and c2r plans */
    plan_FFT_r2c = fftw_plan_dft_r2c_1d(Nx, y_arr_Real, FFT_r2c, FFTW_ESTIMATE);
    plan_iFFT_c2r = fftw_plan_dft_c2r_1d(Nx_r, FFT_r2c, iFFT_c2r, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

    /* Evaluate f(x) , grab real values for r2c */
    TestFunctionOne(Nx, &x_arr[0], &y_arr[0], a);

    for (i = 0; i < Nx; i++)
    {
        y_arr_Real[i] = creal(y_arr[i]);
        y_arr_Imag[i] = cimag(y_arr[i]);
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
    add_yarr_1Dgrouptest(grp_test_id, dataspace_id_c, &y_arr[0], Nx);
    add_FFTc2c_1Dgrouptest(grp_test_id, dataspace_id_c, &FFT_c2c[0], Nx);
    add_iFFTc2c_1Dgrouptest(grp_test_id, dataspace_id_c, &iFFT_c2c[0], Nx);
    add_FFTanalyticc2c_1Dgrouptest(grp_test_id, dataspace_id_c, &FFT_analytic_c2c[0], Nx);
    add_FFTr2c_1Dgrouptest(grp_test_id, dataspace_id_c, &FFT_r2c[0], Nx);
    add_iFFTc2r_1Dgrouptest(grp_test_id, dataspace_id_r, &iFFT_c2r[0]);
    add_FFTanalyticr2c_1Dgrouptest(grp_test_id, dataspace_id_r, &FFT_analytic_r2c[0], Nx_r);

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



