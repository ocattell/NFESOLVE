#ifndef DELAYNEURALFIELDHEADERDEF
#define DELAYNEURALFIELDHEADERDEF

// Include NFESOLVE library and armadillo
#include "NFESOLVE.hpp"
#include <armadillo>

// Class name is DelayNeuralField and it inherits from DDEInterface
class DelayNeuralField : public DDEInterface
{
public:
 	// Constructor takes in parameters, connectivity matrix and mesh
  	DelayNeuralField(const arma::vec& parameters, arma::mat& connectivity, Mesh& mesh);

	// ComputeF method deriving from DDEInterface
	void ComputeF(const double t, const arma::vec& u, const arma::mat& Z, arma::vec& F) const;

    // ComputeHistory method deriving from DDEInterface
    void ComputeHistory(const double t, arma::vec& history) const;

  	// Methods to compute firing rate
  	void ComputeFiringRate(const arma::vec& u, arma::vec& f) const;
    void ComputeFiringRate(const double u, double& f) const;

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
