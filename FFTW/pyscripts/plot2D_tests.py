import numpy as np
import h5py

from pathlib import Path

import matplotlib.pyplot as plt

_ = plt.style.use('dstyle')
_ = plt.style.use('dstyle')

fName = "FFTWFun_out.h5"
fObj = h5py.File(fName, 'r')

TwoD_Tests = fObj['TwoDimensionalTests']
TwoD_attrs = TwoD_Tests.attrs

x_arr, y_arr = TwoD_attrs['x_arr'], TwoD_attrs['y_arr']
Nx, Ny = x_arr.size, y_arr.size
dx, dy = x_arr[1] - x_arr[0], y_arr[1] - y_arr[0]

kx_arr_c, ky_arr_c = TwoD_attrs['kx_arr_c'], TwoD_attrs['ky_arr_c']
kx_arr_r, ky_arr_r = TwoD_attrs['kx_arr_r'], TwoD_attrs['ky_arr_r']

kx_np_c, ky_np_c = np.fft.fftfreq(Nx, dx), np.fft.fftfreq(Ny, dy)
kx_np_r, ky_np_r = np.fft.rfftfreq(Nx, dx), np.fft.rfftfreq(Ny, dy)

bad_kx_r = np.sum(~np.isclose(kx_np_r, kx_arr_r))
bad_kx_c = np.sum(~np.isclose(kx_np_c, kx_arr_c))

bad_ky_r = np.sum(~np.isclose(ky_np_r, ky_arr_r))
bad_ky_c = np.sum(~np.isclose(ky_np_c, ky_arr_c))


bad_kx_r = np.sum(~np.isclose(kx_np_r, kx_arr_r))
bad_ky_r = np.sum(~np.isclose(ky_np_r, ky_arr_r))
bad_kx_c = np.sum(~np.isclose(kx_np_c, kx_arr_c))
bad_ky_c = np.sum(~np.isclose(ky_np_c, ky_arr_c))

assert bad_kx_r == 0
assert bad_kx_r == 0
assert bad_kx_c == 0
assert bad_ky_r == 0


cwd = Path.cwd()


TwoD_testkeys = TwoD_Tests.keys()
    
for testkey in TwoD_testkeys:
    TwoD_Test = TwoD_Tests[testkey]
    print(testkey)
   
    fName = Path(testkey + ".png")
    fPath = cwd / fName

    y_arr = TwoD_Test.get('y_arr_Real')[:] + 1.j * TwoD_Test.get('y_arr_Imag')[:]

    FFT_c2c_arr = TwoD_Test.get('FFT_c2c_Real')[:] + 1.j * TwoD_Test.get('FFT_c2c_Imag')[:]
    FFT_r2c_arr = TwoD_Test.get('FFT_r2c_Real')[:] + 1.j * TwoD_Test.get('FFT_r2c_Imag')[:]

    FFT_analytic_c2c_arr = TwoD_Test.get('FFT_analytic_c2c_Real')[:] + 1.j * TwoD_Test.get('FFT_analytic_c2c_Imag')[:]
    FFT_analytic_r2c_arr = TwoD_Test.get('FFT_analytic_r2c_Real')[:] + 1.j * TwoD_Test.get('FFT_analytic_r2c_Imag')[:]

    iFFT_c2c_arr = TwoD_Test.get('iFFT_c2c_Real')[:] + 1.j * TwoD_Test.get('iFFT_c2c_Imag')[:]
    iFFT_c2r_arr = TwoD_Test.get('iFFT_c2r')[:]

    Pk_c2c_arr = np.abs(FFT_c2c_arr)**2
    Pk_r2c_arr = np.abs(FFT_r2c_arr)**2

    Pk_analytic_c2c_arr = np.abs(FFT_analytic_c2c_arr)**2
    Pk_analytic_r2c_arr = np.abs(FFT_analytic_r2c_arr)**2

    fig, ax = plt.subplots(figsize=(20,6), ncols=4)

    ax_phys, ax_config_Real, ax_config_Imag, ax_PS = ax

    ax_configs = [ax_config_Real, ax_config_Imag, ax_PS]

    _ = ax_phys.scatter(x_arr, y_arr.real, s=1)
    _ = ax_phys.scatter(x_arr, iFFT_c2r_arr / y_arr.size, s=1)
    _ = ax_phys.scatter(x_arr, iFFT_c2c_arr.real / y_arr.size, s=1)


    FFT_real_scale = np.max(FFT_c2c_arr.real) / np.max(FFT_analytic_c2c_arr.real)
    _ = ax_config_Real.scatter(kx_arr_c, FFT_c2c_arr.real, s=1, label='c2c')
    _ = ax_config_Real.scatter(kx_arr_r, FFT_r2c_arr.real, s=1, label='r2c')
    _ = ax_config_Real.scatter(kx_arr_c, FFT_real_scale * FFT_analytic_c2c_arr.real, s=1, label='c2c_analytic')
    _ = ax_config_Real.scatter(kx_arr_r, FFT_real_scale * FFT_analytic_r2c_arr.real, s=1, label='r2c_analytic')


    FFT_imag_scale = np.max(FFT_c2c_arr.imag) / np.max(FFT_analytic_c2c_arr.imag)
    _ = ax_config_Imag.scatter(kx_arr_c, FFT_c2c_arr.imag, s=1)
    _ = ax_config_Imag.scatter(kx_arr_r, FFT_r2c_arr.imag, s=1)
    _ = ax_config_Imag.scatter(kx_arr_c, FFT_imag_scale * FFT_analytic_c2c_arr.imag, s=1)
    _ = ax_config_Imag.scatter(kx_arr_r, FFT_imag_scale * FFT_analytic_r2c_arr.imag, s=1)


    _ = ax_PS.scatter(kx_arr_c, Pk_c2c_arr / FFT_real_scale / FFT_real_scale, s=1)
    _ = ax_PS.scatter(kx_arr_r, Pk_r2c_arr / FFT_real_scale / FFT_real_scale, s=1)
    _ = ax_PS.scatter(kx_arr_c, Pk_analytic_c2c_arr, s=1)
    _ = ax_PS.scatter(kx_arr_r, Pk_analytic_r2c_arr, s=1)


    _ = ax_phys.set_xlim(-2.,2.)
    _ = ax_config_Real.set_xlim(-5.,5.)
    _ = ax_config_Imag.set_xlim(-5.,5.)
    _ = ax_PS.set_xlim(-5.,5.)


    _ = ax_phys.set_xlabel(r'$x$')
    _ = ax_phys.set_ylabel(r'$f(x)$')

    _ = ax_config_Real.set_ylabel(r'$Re[\widetilde{f}(k)]$')
    _ = ax_config_Imag.set_ylabel(r'$Im[\widetilde{f}(k)]$')

    _ = ax_PS.set_ylabel(r'$P(k)$')


    for sub_ax in ax_configs:
        _ = sub_ax.set_xlabel(r'$k_x$')

    for sub_ax in ax:
        _ = sub_ax.grid(which='both', alpha=0.3)
        _ = sub_ax.axvline(x=0., c='k', ls='--', lw=1.)
        _ = sub_ax.axhline(y=0., c='k', ls='--', lw=1.)

        
    _ = ax_config_Real.legend()
    
    plt.tight_layout()
    
    plt.savefig(fPath, dpi=256)
    #plt.show()



