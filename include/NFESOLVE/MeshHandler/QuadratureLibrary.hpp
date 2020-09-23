#ifndef QUADRATURELIBRARYHEADERDEF
#define QUADRATURELIBRARYHEADERDEF

#include "Mesh.hpp"
#include "Regular2DQuadGrid.hpp"
#include <armadillo>

namespace QuadratureLibrary
{
    // Generate weights for vertex quadrature rule
  void GenerateVQWeights(const Mesh& mesh, arma::vec& weights);
  void GenerateVQWeights(const Regular2DQuadGrid& mesh, arma::vec& weights);

  // Generate points and weights for n point (per element) Gauss quadrature rule
  void GenerateGQPointsWeights(const Mesh& mesh, const int n, arma::mat& GQPoints, arma::vec& weights);
};

#endif
