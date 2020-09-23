#ifndef DDE_EXAMPLE1HEADERDEF
#define DDE_EXAMPLE1HEADERDEF

// Include NFESOLVE library and armadillo
#include "NFESOLVE.hpp"
#include <armadillo>

// Class name is DDE_Example1 and it inherits from DDEInterface
class DDE_Example1 : public DDEInterface
{
public:
	// Constructor takes in parameters vector and Amplitudes Am and Ap
	DDE_Example1(const arma::vec& parameters, const double Am, const double Ap);

	// ComputeF method deriving from DDEInterface
	void ComputeF(const double t, const arma::vec& u, const arma::mat& Z, arma::vec& F) const;

	// ComputeHistory method deriving from DDEInterface
	void ComputeHistory(const double t, arma::vec& history) const;

private:
	// Private storage for parameters vector
	arma::vec mParameters;
	double mAm, mAp;
};

#endif
