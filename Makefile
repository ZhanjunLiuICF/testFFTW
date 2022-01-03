all:
	gfortran -ffree-line-length-200  tstfftw.f90  -lfftw3  -lm -I/usr/include
