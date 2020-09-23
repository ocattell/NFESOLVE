#ifndef ABSTRACTDDESOLVERHEADERDEF
#define ABSTRACTDDESOLVERHEADERDEF

#include "../DE/AbstractDESolver.hpp"
#include "DDEInterface.hpp"
#include <armadillo>

class AbstractDDESolver : public AbstractDESolver
{
public:

  // Set DDE system
  void SetDDESystem(DDEInterface& aDDESystem);

  // Set number of delay equations
  void SetNumDelayEqs(int numDelayEqs);

  // Set delays
  void SetDelays(arma::vec& delays);

protected:

  // DDE system
  DDEInterface* mpDDESystem;

  // Number of equations containing delays
  int mNumDelayEqs;

  // Delay vector
  arma::vec mDelays;
};

#endif
