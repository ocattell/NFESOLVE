#ifndef MESHHEADERDEF
#define MESHHEADERDEF

#include <armadillo>
#include <string>

class Mesh
{
public:
  // Saves mesh data to a file in the same format as reading in with MeshFromFile
  void SaveToFile(const std::string fileName);

  // Sets grid points
  void SetGridPoints(const arma::mat& gridPoints);

  // Sets element connectivity array
  void SetElementConnectivity(const arma::umat& elementConnectivity);

  // Getters for all grid points and specific grid points
  const arma::mat& GetGridPoints() const;
  arma::subview_row<double> GetGridPoint(const arma::uword gridPointNumber) const;
  arma::mat GetElementGridPoints(const arma::uword elementNumber) const;

  // Getters for whole element connectivity array and specific element connectivity array
  const arma::umat& GetElementConnectivity() const;
  arma::urowvec GetElementConnectivity(const int elementNumber) const;

  // Getters for dimension, number of grid points, number of element nodes and number of elements
  arma::uword GetDimension() const;
  arma::uword GetNumGridPoints() const;
  arma::uword GetNumElementNodes() const;
  arma::uword GetNumElements() const;

protected:

  // Storage for dimension, number of grid points, number of element nodes and number of elements
  arma::uword mDimension;
  arma::uword mNumGridPoints;
  arma::uword mNumElementNodes;
  arma::uword mNumElements;

  // Grid points and element connectivity array
  arma::mat mGridPoints;
  arma::umat mElementConnectivity;
};

#endif
