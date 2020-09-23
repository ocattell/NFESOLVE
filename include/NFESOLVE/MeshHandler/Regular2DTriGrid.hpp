#ifndef REGULAR2DTRIGRIDHEADERDEF
#define REGULAR2DTRIGRIDHEADERDEF

#include "Mesh.hpp"
#include <armadillo>

class Regular2DTriGrid : public Mesh
{
public:
    // Constructors for setting up regular 2D triangulated grid mesh
    Regular2DTriGrid(const double x0, const double xn, const double y0, const double yn, const arma::uword numXYPoints);
    Regular2DTriGrid(const double x0, const double xn, const double y0, const double yn, const arma::uword numXPoints, const arma::uword numYPoints);
    Regular2DTriGrid(const arma::vec x, const arma::vec y);
};

#endif
