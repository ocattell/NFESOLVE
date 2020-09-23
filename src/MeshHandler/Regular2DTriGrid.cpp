#include "Debug.hpp"
#include "Regular2DTriGrid.hpp"

Regular2DTriGrid::Regular2DTriGrid(const double x0, const double xn, const double y0, const double yn, const arma::uword numXYPoints)
{
    DebugCheck(x0 >= xn, "Regular2DTriGrid(): x0 must be less than xn");
    DebugCheck(y0 >= yn, "Regular2DTriGrid(): y0 must be less than yn");
    DebugCheck(numXYPoints <= 0, "Regular2DTriGrid(): numXYPoints must be positive");

    mDimension = 2;
    mNumGridPoints = numXYPoints*numXYPoints;
    mNumElementNodes = 3;
    mNumElements = (numXYPoints-1)*(numXYPoints-1)*2;

    arma::vec xPoints = arma::linspace(x0, xn, numXYPoints);
    arma::vec yPoints = arma::linspace(y0, yn, numXYPoints);

    mGridPoints.zeros(mNumGridPoints,mDimension);
    mElementConnectivity.zeros(mNumElements,mNumElementNodes);

    arma::uword elementNum = 0;
    for (arma::uword j=0; j<numXYPoints; ++j)
    {
        for (arma::uword i=0; i<numXYPoints; ++i)
        {
            mGridPoints(j*numXYPoints + i,0) = xPoints(i);
            mGridPoints(j*numXYPoints + i,1) = yPoints(j);

            if (i != numXYPoints-1 && j != numXYPoints-1)
            {
                mElementConnectivity(elementNum,0) = (j*numXYPoints)+i;
                mElementConnectivity(elementNum,1) = (j*numXYPoints)+i+1;
                mElementConnectivity(elementNum,2) = ((j+1)*numXYPoints)+i;

                mElementConnectivity(elementNum+1,0) = (j*numXYPoints)+i+1;
                mElementConnectivity(elementNum+1,1) = ((j+1)*numXYPoints)+i+1;
                mElementConnectivity(elementNum+1,2) = ((j+1)*numXYPoints)+i;
                elementNum += 2;
            }
        }
    }
}

Regular2DTriGrid::Regular2DTriGrid(const double x0, const double xn, const double y0, const double yn, const arma::uword numXPoints, const arma::uword numYPoints)
{
    DebugCheck(x0 >= xn, "Regular2DTriGrid(): x0 must be less than xn");
    DebugCheck(y0 >= yn, "Regular2DTriGrid(): y0 must be less than yn");
    DebugCheck(numXPoints <= 0, "Regular2DTriGrid(): numXPoints must be positive");
    DebugCheck(numYPoints <= 0, "Regular2DTriGrid(): numYPoints must be positive");

    mDimension = 2;
    mNumGridPoints = numXPoints*numYPoints;
    mNumElementNodes = 3;
    mNumElements = (numXPoints-1)*(numYPoints-1)*2;

    arma::vec xPoints = arma::linspace(x0, xn, numXPoints);
    arma::vec yPoints = arma::linspace(y0, yn, numYPoints);

    mGridPoints.zeros(mNumGridPoints,mDimension);
    mElementConnectivity.zeros(mNumElements,mNumElementNodes);

    arma::uword elementNum = 0;
    for (arma::uword j=0; j<numYPoints; ++j)
    {
        for (arma::uword i=0; i<numXPoints; ++i)
        {
            mGridPoints(j*numXPoints + i,0) = xPoints(i);
            mGridPoints(j*numXPoints + i,1) = yPoints(j);

            if (i != numXPoints-1 && j != numYPoints-1)
            {
                mElementConnectivity(elementNum,0) = (j*numXPoints)+i;
                mElementConnectivity(elementNum,1) = (j*numXPoints)+i+1;
                mElementConnectivity(elementNum,2) = ((j+1)*numXPoints)+i;

                mElementConnectivity(elementNum+1,0) = (j*numXPoints)+i+1;
                mElementConnectivity(elementNum+1,1) = ((j+1)*numXPoints)+i+1;
                mElementConnectivity(elementNum+1,2) = ((j+1)*numXPoints)+i;
                elementNum += 2;
            }
        }
    }
}

Regular2DTriGrid::Regular2DTriGrid(const arma::vec x, const arma::vec y)
{
    DebugCheck(!(x.is_sorted("strictascend")), "Regular2DTriGrid(): x vector must be in ascending order and cannot have repeated values");
    DebugCheck(!(y.is_sorted("strictascend")), "Regular2DTriGrid(): y vector must be in ascending order and cannot have repeated values");
    DebugCheck(x.n_elem < 2, "Regular2DTriGrid(): x must have two or more elements");
    DebugCheck(y.n_elem < 2, "Regular2DTriGrid(): y must have two or more elements");

    mDimension = 2;
    mNumGridPoints = x.n_elem*y.n_elem;
    mNumElementNodes = 3;
    mNumElements = (x.n_elem-1)*(y.n_elem-1)*2;

    mGridPoints.zeros(mNumGridPoints,mDimension);
    mElementConnectivity.zeros(mNumElements,mNumElementNodes);

    arma::uword elementNum = 0;
    for (arma::uword j=0; j<y.n_elem; ++j)
    {
        for (arma::uword i=0; i<x.n_elem; ++i)
        {
            mGridPoints(j*x.n_elem + i,0) = x(i);
            mGridPoints(j*x.n_elem + i,1) = y(j);

            if (i != x.n_elem-1 && j != y.n_elem-1)
            {
                mElementConnectivity(elementNum,0) = (j*x.n_elem)+i;
                mElementConnectivity(elementNum,1) = (j*x.n_elem)+i+1;
                mElementConnectivity(elementNum,2) = ((j+1)*x.n_elem)+i;

                mElementConnectivity(elementNum+1,0) = (j*x.n_elem)+i+1;
                mElementConnectivity(elementNum+1,1) = ((j+1)*x.n_elem)+i+1;
                mElementConnectivity(elementNum+1,2) = ((j+1)*x.n_elem)+i;
                elementNum += 2;
            }
        }
    }
}
