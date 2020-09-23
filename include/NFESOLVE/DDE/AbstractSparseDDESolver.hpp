#ifndef ABSTRACTSPARSEDDESOLVERHEADERDEF
#define ABSTRACTSPARSEDDESOLVERHEADERDEF

#include "../DE/AbstractDESolver.hpp"
#include "SparseDDEInterface.hpp"
#include <armadillo>

class AbstractSparseDDESolver : public AbstractDESolver
{
public:

  // Set DDE system
  void SetDDESystem(SparseDDEInterface& aSparseDDESystem);

  // Set sparsity pattern of Z
  void SetZLocations(const arma::umat& ZLocations);

  // Set number of delay equations
  void SetNumDelayEqs(const int numDelayEqs);

  // Set delays
  void SetDelays(const arma::vec& delays);

protected:

  // DDE system
  SparseDDEInterface* mpSparseDDESystem;

  // ZLocations
  arma::umat mZLocations;

  // Number of equations containing delays
  int mNumDelayEqs;

  // Delay vector
  arma::vec mDelays;
};

#endif
