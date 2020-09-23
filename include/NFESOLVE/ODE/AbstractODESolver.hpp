#ifndef ABSTRACTODESOLVERHEADERDEF
#define ABSTRACTODESOLVERHEADERDEF

#include "../DE/AbstractDESolver.hpp"
#include "ODEInterface.hpp"
#include <armadillo>

class AbstractODESolver : public AbstractDESolver
{
public:

  // Set ODE system
  void SetODESystem(ODEInterface& anODESystem);

protected:

  // ODE system
  ODEInterface* mpODESystem;
};

#endif
