#include "AbstractDDESolver.hpp"

// Set ODE system
void AbstractDDESolver::SetDDESystem(DDEInterface& aDDESystem)
{
  mpDDESystem = &aDDESystem;
}

void AbstractDDESolver::SetNumDelayEqs(int numDelayEqs)
{
	mNumDelayEqs = numDelayEqs;
}

void AbstractDDESolver::SetDelays(arma::vec& delays)
{
	mDelays = delays;
}