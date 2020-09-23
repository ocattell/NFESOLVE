#ifndef SPARSEDDEINTERFACEHEADERDEF
#define SPARSEDDEINTERFACEHEADERDEF

#include "SparseDelayMatrix.hpp"
#include <armadillo>

class SparseDDEInterface
{
public:

	// Compute right-hand side using sparse matrix
    virtual void ComputeF(const double t, const arma::vec& u, const SparseDelayMatrix& Z, arma::vec& F) const = 0;

	// Compute history
	virtual void ComputeHistory(const double t, arma::vec& history) const = 0;

	// Compute analytical solution
	virtual void ComputeAnalyticSolution(const double t, arma::vec& u) const;
};

#endif
