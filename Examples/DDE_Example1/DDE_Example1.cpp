// Include DDE_Example1 header
#include "DDE_Example1.hpp"

// Constructor
DDE_Example1::DDE_Example1(const arma::vec& parameters, const double Am, const double Ap)
{
	// Set private parameters vector and amplitudes
	mParameters = parameters;
	mAm = Am;
	mAp = Ap;
}

// ComputeF
void DDE_Example1::ComputeF(const double t, const arma::vec& u, const arma::mat& Z, arma::vec& F) const
{
	double f = 0.0;
	if (Z(0,0) >= -mParameters(2) && Z(0,0) <= mParameters(2))
	{
		f = mParameters(1);
	}
	F = f - mParameters(0)*u;
}

// ComputeHistory
void DDE_Example1::ComputeHistory(const double t, arma::vec& history) const
{
	history = { mAp };
}
