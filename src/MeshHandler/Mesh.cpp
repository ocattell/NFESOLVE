#include "Debug.hpp"
#include "Mesh.hpp"

#include <iomanip>
#include <fstream>

void Mesh::SaveToFile(const std::string fileName)
{
    std::ofstream meshFile;
    meshFile.open(fileName);
    DebugCheck(!meshFile.is_open(), "Mesh file is not open");

    meshFile << mDimension << std::endl
             << mNumGridPoints << std::endl
             << mNumElementNodes << std::endl
             << mNumElements << std::endl;

    meshFile << std::endl;

    std::streamsize p = meshFile.precision();
    meshFile.precision(10);

    for (arma::uword i = 0; i < mNumGridPoints; ++i)
    {
        for (arma::uword j = 0; j < mDimension; ++j)
        {
            meshFile << mGridPoints(i, j) << " ";
        }
        meshFile << std::endl;
    }
    meshFile << std::endl;

    meshFile.precision(p);
    for (arma::uword i = 0; i < mNumElements; ++i)
    {
        for (arma::uword j = 0; j < mNumElementNodes; ++j)
        {
            meshFile << mElementConnectivity(i, j) << " ";
        }
        meshFile << std::endl;
    }
    meshFile.close();
}

void Mesh::SetGridPoints(const arma::mat& gridPoints)
{
    mGridPoints = gridPoints;
    mDimension = gridPoints.n_cols;
    mNumGridPoints = gridPoints.n_rows;
}

void Mesh::SetElementConnectivity(const arma::umat& elementConnectivity)
{
    mElementConnectivity = elementConnectivity;
    mNumElementNodes = elementConnectivity.n_cols;
    mNumElements = elementConnectivity.n_rows;
}

const arma::mat& Mesh::GetGridPoints() const
{
    return mGridPoints;
}

arma::subview_row<double>  Mesh::GetGridPoint(const arma::uword gridPointNumber) const
{
    DebugCheck(gridPointNumber < 0 || gridPointNumber >= mNumGridPoints, "Mesh(): gridPointNumber not within correct bounds");
    return mGridPoints.row(gridPointNumber);
}

arma::mat Mesh::GetElementGridPoints(const arma::uword elementNumber) const
{
    DebugCheck(elementNumber < 0 || elementNumber >= mNumElements, "Mesh(): elementNumber not within correct bounds");
    arma::mat elementGridPoints(mNumElementNodes,mDimension);
    for (arma::uword i=0; i<mNumElementNodes; ++i)
    {
        elementGridPoints.row(i) = mGridPoints.row(mElementConnectivity(elementNumber,i));
    }
    return elementGridPoints;
}

const arma::umat& Mesh::GetElementConnectivity() const
{
    return mElementConnectivity;
}

arma::urowvec Mesh::GetElementConnectivity(const int elementNumber) const
{
    DebugCheck(elementNumber < 0 || elementNumber >= mNumElements, "Mesh(): elementNumber not within correct bounds");
    return mElementConnectivity.row(elementNumber);
}

arma::uword Mesh::GetDimension() const
{
  return mDimension;
}

arma::uword Mesh::GetNumGridPoints() const
{
  return mNumGridPoints;
}

arma::uword Mesh::GetNumElementNodes() const
{
  return mNumElementNodes;
}

arma::uword Mesh::GetNumElements() const
{
  return mNumElements;
}
