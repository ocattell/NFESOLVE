#ifndef NFESOLVEHEADERDEF
#define NFESOLVEHEADERDEF

#include <armadillo>
#include <iostream>
#include <stdexcept>
#include <cmath>

#include "NFESOLVE/CompilerSetup.hpp"

#if defined(NFESOLVE_USE_OPENMP)
  #include <omp.h>
#endif

#include "NFESOLVE/Debug.hpp"
#include "NFESOLVE/DE/DE.hpp"
#include "NFESOLVE/DDE/DDE.hpp"
#include "NFESOLVE/MeshHandler/MeshHandler.hpp"
#include "NFESOLVE/ODE/ODE.hpp"

#endif
