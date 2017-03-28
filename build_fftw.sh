#!/bin/sh

cd fftw-2.1.5/sourceAndDoc

# Make the fftw library
make clean
make

cd ..

# Copy the compiled library and headers to the place where smoke includes it.
cp sourceAndDoc/fftw/.libs/libfftw.a lib
cp sourceAndDoc/fftw/libfftw.la lib
cp sourceAndDoc/rfftw/.libs/librfftw.a lib
cp sourceAndDoc/rfftw/librfftw.la lib
cp sourceAndDoc/fftw/fftw.h include
cp sourceAndDoc/rfftw/rfftw.h include
