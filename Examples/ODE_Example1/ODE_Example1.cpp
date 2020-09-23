// Include ODE_Example1 header
#include "ODE_Example1.hpp"

// Constructor
ODE_Example1::ODE_Example1(const arma::vec& parameters)
{
	// Set private parameters vector
	mParameters = parameters;
}

// ComputeF
void ODE_Example1::ComputeF(const double t, const arma::vec& u, arma::vec& F) const
{
	// Clearly define parameter names
	double a = mParameters(0);
	double b = mParameters(1);

	// Right-hand side definition
	F(0) = -a*u(0);
	F(1) = b*u(0) + u(1);
	F(2) = -u(2) + b*t;
}
