# FFTW_Fun
Basic examples of using FFTW with HDF5 files

The general goal is to get comfortable with the Fourier Transform routines in four ways building on one another:

1. Regular FFTW
2. FFTW with openmp (distributed computing / multi-threading)
3. FFTW with MPI (distributed memory)
4. FFTW with MPI and openmp


For each of these building blocks, we will compute the Fourier Transform of a couple functions to get comfortable with all routines

Test 1:

$$
f(x) = \textrm{rect}(ax)
$$

with Fourier Transform

$$
\widetilde{f}(k) = \frac{1}{|a|} \sinc{\frac{k}{a}}
$$

Test 2:

$$
f(x) = \exp\left(-\alpha x^2 \right)
$$

with Fourier Transform

$$
\widetilde{f}(k) = \sqrt{\frac{\pi}{\alpha}} \exp\left[- \left(\pi k\right)^2 / \alpha\right]
$$


Test 3:

$$
f(x) = \exp(-a|x|)
$$

with Fourier Transform

$$
\widetilde{f}(k) = \frac{2a}{a^2 + 4 \pi^2 k^2}
$$


Test 4:

$$
f(x) = \cos(ax)
$$

with Fourier Transform

$$
\widetilde{f}(k) = \frac{1}{2} \left[\delta\left(k - \frac{a}{2\pi} \right) + \delta\left(k + \frac{a}{2\pi} \right)\right]
$$

where $\delta(x)$ is the Dirac-Delta distribution

Test 5:

$$
f(x) = \sin(ax^2)
$$

with Fourier Transform

$$
\widetilde{f}(k) = - \sqrt{\frac{\pi}{a}} \sin \left(\frac{\pi^2 k^2}{a} - \frac{\pi}{4}\right)
$$


Test 6:

$$
f(x,y) = \exp\left[-\pi \left(a^2 x^2 + b^2 y^2 \right) \right]
$$

with Fourier Transform

$$
\widetilde{f}(k_x, k_y) = \frac{1}{|ab|} \exp\left[-\pi \left(\frac{k_x^2}{a^2} + \frac{k_y^2}{b^2} \right) \right]
$$

Test 7:

$$
f(x,y) = \frac{1}{\sqrt{x^2 + y^2}}
$$

with Fourier Transform

$$
\widetilde{f}(k_x, k_y) = \frac{1}{\sqrt{k_x^2 + k_y^2}}
$$


Test 8:

$$
f(x,y) = \frac{i}{x + iy}
$$

with Fourier Transform

$$
\widetilde{f}(k_x, k_y) = \frac{1}{k_x + i k_y}
$$


For the one-dimensional cases, we will use array from $x \in (-5,5)$ with 2048 grid points.
For the two-dimensional cases, we will use domain with same x values and y values $y \in (-5,5)$ with (1024,1024) grid points.


In the future, would like to also test Heaviside function (equations 205 and 315), Bessel function (equation 317), and log|x| (equation 319), as well as exponential of absolute value (equation 504) and Reisz potential (equation 502). Source of functions and their Fourier Transforms: https://en.wikipedia.org/wiki/Fourier_transform


