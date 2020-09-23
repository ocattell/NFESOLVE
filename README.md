# NFESOLVE
NFESOLVE is a suite of differential equation solvers with a focus on the efficient computation of delay differential equations. 


1. Make sure Armadillo is installed with an OpenBLAS backend (OpenBLAS isn't strictly necessary but provides optimised BLAS routines that will make code faster).

2. If you are wanting to run the parallel version of the code then you must have OpenMP 3.1 or later installed, as detailed in the Armadillo documentation.

3. Open terminal and cd to NFESOLVE directory. To compile library type "make NFESOLVE", or for the parallel version type "make NFESOLVE_PAR". If you have OpenMP installed and wish to compile both the serial and parallel versions just type "make".

4. To include the NFESOLVE library type '#include "NFESOLVE.hpp"' in the header definitions of your files.

4. See 'Examples' folder for a variety of ODE and DDE example problems. For each problem a Makefile is supplied that follows the same template. These can be copied and edited to suit your own usage. The DelayNFE_Example1 and SparseDelayNFE_Example1 makefiles contain both serial and parallel versions which can be made individually by adding or removing the "_PAR" suffix to the make target when making.
