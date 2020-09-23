#include "AbstractSparseDDESolver.hpp"

// Set DDE system
void AbstractSparseDDESolver::SetDDESystem(SparseDDEInterface& aSparseDDESystem)
{
    mpSparseDDESystem = &aSparseDDESystem;
}

// Set sparsity pattern for Z
void AbstractSparseDDESolver::SetZLocations(const arma::umat& ZLocations)
{
    if (ZLocations.n_rows != 2)
    {
        throw std::length_error("ZLocations should be of size 2xN");
    }
    mZLocations = ZLocations;
}

void AbstractSparseDDESolver::SetNumDelayEqs(const int numDelayEqs)
{
	mNumDelayEqs = numDelayEqs;
}

void AbstractSparseDDESolver::SetDelays(const arma::vec& delays)
{
	mDelays = delays;
}