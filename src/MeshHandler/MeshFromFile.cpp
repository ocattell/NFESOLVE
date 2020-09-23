#include "Debug.hpp"
#include "MeshFromFile.hpp"
#include <fstream>

MeshFromFile::MeshFromFile(const std::string fileName)
{
    DebugCheck(fileName.substr(fileName.length() - 4, 4) != ".dat", "File type should be .dat");

    std::ifstream meshFile;
    meshFile.open(fileName);
    DebugCheck(!meshFile.is_open(), "Mesh file is not open");

    meshFile >> mDimension >> mNumGridPoints >> mNumElementNodes >> mNumElements;
    mGridPoints.zeros(mNumGridPoints, mDimension);
    mElementConnectivity.zeros(mNumElements, mNumElementNodes);

    for (arma::uword i = 0; i < mNumGridPoints; ++i)
    {
        for (arma::uword j = 0; j < mDimension; ++j)
        {
            meshFile >> mGridPoints(i, j);
        }
    }
    for (arma::uword i = 0; i < mNumElements; ++i)
    {
        for (arma::uword j = 0; j < mNumElementNodes; ++j)
        {
            meshFile >> mElementConnectivity(i, j);
        }
    }
    meshFile.close();
}