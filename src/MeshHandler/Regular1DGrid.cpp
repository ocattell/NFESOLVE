#include "Debug.hpp"
#include "Regular1DGrid.hpp"

Regular1DGrid::Regular1DGrid(const double x0, const double xn, const arma::uword numXPoints)
{
    DebugCheck(x0 >= xn, "Regular1DGrid(): x0 must be less than xn");
    DebugCheck(numXPoints <= 0, "Regular1DGrid(): numXPoints must be positive");

    mDimension = 1;
    mNumGridPoints = numXPoints;
    mNumElementNodes = 2;
    mNumElements = (numXPoints - 1);

    mGridPoints.zeros(mNumGridPoints,mDimension);
    mElementConnectivity.zeros(mNumElements,mNumElementNodes);

    mGridPoints.col(0) = arma::linspace(x0, xn, numXPoints);
    for (arma::uword i=0; i<mElementConnectivity.n_rows; ++i)
    {
        mElementConnectivity(i,0) = i;
        mElementConnectivity(i,1) = i+1;
    }
}

Regular1DGrid::Regular1DGrid(const arma::vec x)
{
    DebugCheck(!(x.is_sorted("strictascend")), "Regular1DGrid(): x vector must be in ascending order and cannot have repeated values");
    DebugCheck(x.n_elem < 2, "Regular1DGrid(): x must have two or more elements");

    mDimension = 1;
    mNumGridPoints = x.n_elem;
    mNumElementNodes = 2;
    mNumElements = (x.n_elem - 1);

    mGridPoints.zeros(mNumGridPoints,mDimension);
    mElementConnectivity.zeros(mNumElements,mNumElementNodes);

    mGridPoints.col(0) = x;
    for (arma::uword i=0; i<mElementConnectivity.n_rows; ++i)
    {
        mElementConnectivity(i,0) = i;
        mElementConnectivity(i,1) = i+1;
    }
}
