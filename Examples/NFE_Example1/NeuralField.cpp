// Include NeuralField
#include "NeuralField.hpp"

// Constructor
NeuralField::NeuralField(const arma::vec& parameters, arma::mat& connectivity, Mesh& mesh)
{
		mParameters = parameters;
  	mpConnectivity = &connectivity;
  	mpMesh = &mesh;
}

// ComputeF
void NeuralField::ComputeF(const double t, const arma::vec& u, arma::vec& F) const
{
  	arma::vec firingRate;
  	ComputeFiringRate(u, firingRate);
  	F = -u + (*mpConnectivity)*firingRate;
}

// Method to compute firing rate
void NeuralField::ComputeFiringRate(const arma::vec& u, arma::vec& f) const
{
	double mu = mParameters(0);
	double theta = mParameters(1);
 	f = 1.0/(1.0 + exp(-(mu*(u - theta))));
}

// Set parameters
void NeuralField::SetParameters(const arma::vec& parameters)
{
  mParameters = parameters;
}

// Set connectivity matrix reference
void NeuralField::SetConnectivityKernel(arma::mat& connectivity)
{
  mpConnectivity = &connectivity;
}

// Set mesh reference
void NeuralField::SetMesh(Mesh& mesh)
{
  mpMesh = &mesh;
}
