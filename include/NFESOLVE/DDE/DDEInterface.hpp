#ifndef DDEINTERFACEHEADERDEF
#define DDEINTERFACEHEADERDEF

#include <armadillo>

class DDEInterface
{
public:

  // Compute right-hand side
  virtual void ComputeF(const double t, const arma::vec& u, const arma::mat& Z, arma::vec& F) const = 0;

  // Compute history
  virtual void ComputeHistory(const double t, arma::vec& history) const = 0;

  // Compute analytical solution
  virtual void ComputeAnalyticSolution(const double t, arma::vec& u) const;

};

#endif
