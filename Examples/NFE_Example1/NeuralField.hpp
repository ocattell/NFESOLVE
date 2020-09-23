#ifndef NEURALFIELDHEADERDEF
#define NEURALFIELDHEADERDEF

// Include NFESOLVE library and armadillo
#include "NFESOLVE.hpp"
#include <armadillo>

// Class name is NeuralField and it inherits from ODEInterface
class NeuralField : public ODEInterface
{
public:
 	// Constructor takes in parameters, connectivity matrix and mesh
  	NeuralField(const arma::vec& parameters, arma::mat& connectivity, Mesh& mesh);

	// ComputeF method deriving from ODEInterface
	void ComputeF(const double t, const arma::vec& u, arma::vec& F) const;

  	// Method to compute firing rate
  	void ComputeFiringRate(const arma::vec& u, arma::vec& f) const;

  	// Set mParameters
  	void SetParameters(const arma::vec& parameters);

 	// Set connectivity matrix
  	void SetConnectivityKernel(arma::mat& connectivity);

	// Set mesh
	void SetMesh(Mesh& mesh);

private:
  	// Private storage for parameters
 	arma::vec mParameters;
	// Private pointers to connectivity matrix and mesh
  	arma::mat* mpConnectivity;
  	Mesh* mpMesh;
};

#endif
