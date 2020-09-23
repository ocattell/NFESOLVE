#include "AbstractODESolver.hpp"

// Set ODE system
void AbstractODESolver::SetODESystem(ODEInterface& anODESystem)
{
  mpODESystem = &anODESystem;
}