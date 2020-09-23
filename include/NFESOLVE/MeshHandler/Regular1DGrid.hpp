#ifndef REGULAR1DGRIDHEADERDEF
#define REGULAR1DGRIDHEADERDEF

#include "Mesh.hpp"
#include <armadillo>

class Regular1DGrid : public Mesh
{
public:
    // Constructor for setting up regular 1D grid mesh
    Regular1DGrid(const double x0, const double xn, const arma::uword numXPoints);
    Regular1DGrid(const arma::vec x);
};

#endif
