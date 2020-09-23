#ifndef REGULAR2DQUADGRIDHEADERDEF
#define REGULAR2DQUADGRIDHEADERDEF

#include "Mesh.hpp"
#include <armadillo>

class Regular2DQuadGrid : public Mesh
{
public:
    // Constructors for setting up regular 2D quadrilateral grid mesh
    Regular2DQuadGrid(const double x0, const double xn, const double y0, const double yn, const arma::uword numXYPoints);
    Regular2DQuadGrid(const double x0, const double xn, const double y0, const double yn, const arma::uword numXPoints, const arma::uword numYPoints);
    Regular2DQuadGrid(const arma::vec x, const arma::vec y);
};

#endif
