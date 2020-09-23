// Include DelayNeuralField
#include "DelayNeuralField.hpp"
#include <cmath>

// Constructor
DelayNeuralField::DelayNeuralField(const arma::vec& parameters, arma::mat& connectivity, Mesh& mesh)
{
	mParameters = parameters;
  	mpConnectivity = &connectivity;
  	mpMesh = &mesh;
}

// ComputeF
void DelayNeuralField::ComputeF(const double t, const arma::vec& u, const arma::mat& Z, arma::vec& F) const
{
	double firingRate;

	#pragma omp parallel for default(shared) private(firingRate)
	for (arma::uword i = 0; i < Z.n_rows; ++i)
	{
		double sum = 0;
		for (arma::uword j = 0; j < Z.n_rows; ++j)
		{
			if (i == j)
			{
				ComputeFiringRate(Z(j, 0), firingRate);
			}
			else if (i > j)
			{
				int k = round((2 * j * Z.n_rows - pow(j, 2.0) + 2 * i - 3 * j - 2) / 2) + 1;
				ComputeFiringRate(Z(j, k), firingRate);
			}
			else
			{
				int k = round((2 * i * Z.n_rows - pow(i, 2.0) + 2 * j - 3 * i - 2) / 2) + 1;
				ComputeFiringRate(Z(j, k), firingRate);
			}

			sum += (*mpConnectivity)(i, j) * firingRate;
		}
		F(i) = -u(i) + sum;
	}
}

// ComputeHistory
void DelayNeuralField::ComputeHistory(const double t, arma::vec& history) const
{
    history.zeros();
}

// Methods to compute firing rate
void DelayNeuralField::ComputeFiringRate(const arma::vec& u, arma::vec& f) const
{
	double mu = mParameters(0);
	double theta = mParameters(1);
 	f = 1.0/(1.0 + exp(-(mu*(u - theta))));
}
void DelayNeuralField::ComputeFiringRate(const double u, double& f) const
{
	double mu = mParameters(0);
	double theta = mParameters(1);
	f = 1.0 / (1.0 + exp(-(mu * (u - theta))));
}

// Set parameters
void DelayNeuralField::SetParameters(const arma::vec& parameters)
{
  mParameters = parameters;
}

// Set connectivity matrix reference
void DelayNeuralField::SetConnectivityKernel(arma::mat& connectivity)
{
  mpConnectivity = &connectivity;
}

// Set mesh reference
void DelayNeuralField::SetMesh(Mesh& mesh)
{
  mpMesh = &mesh;
}
