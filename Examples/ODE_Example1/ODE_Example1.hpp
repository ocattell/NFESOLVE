#ifndef ODE_EXAMPLE1HEADERDEF
#define ODE_EXAMPLE1HEADERDEF

// Include NFESOLVE library and armadillo
#include "NFESOLVE.hpp"
#include <armadillo>

// Class name is ODE_Example1 and it inherits from ODEInterface
class ODE_Example1 : public ODEInterface
{
public:
	// Constructor takes in parameters vector
	ODE_Example1(const arma::vec& parameters);

	// ComputeF method deriving from ODEInterface
	void ComputeF(const double t, const arma::vec& u, arma::vec& F) const;

private:
	// Private storage for parameters vector
	arma::vec mParameters;
};

#endif
