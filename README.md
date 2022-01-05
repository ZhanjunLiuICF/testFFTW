# testFFTW
使用FFTW_Plan_many_dft的例子，针对三维数组中的任意两维做2D FFT。

* 1. 数组比较小的情况下，采用FFTW_plan_dft经过多次变换和采用FFTW_Plan_many_dft的结果相同。
* 2. 数组比较大的情况下，采用FFTW_plan_dft经过多次变换和采用FFTW_Plan_many_dft的结果略微不同。舍入误差的因素？

