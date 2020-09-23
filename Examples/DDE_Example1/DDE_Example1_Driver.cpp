//// DDE_Example1 solves a DDE with a periodic solution using a fixed step third order Runge-Kutta solver ////

// Include relevant headers
#include "DDE_Example1.hpp"
#include "NFESOLVE.hpp"

#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	// Set up solver parameters
	double initialTime = 0.0;
	double finalTime = 10;
	double stepSize = 0.001;
	int numDelayEqs = 1;
	std::string outputFileName = "DDE_Example1.dat";
	int saveGap = 10;
	int printGap = 1;

	// Set up DDE parameters
	double alpha = 1.0;
	double gamma = 5.0;
	double tau = 10.0;
	double theta = 1.0;
	arma::vec P = {alpha, gamma, theta};

	// Set up amplitudes
	double Am = exp(-alpha*tau)*theta;
	double Ap = (exp(-alpha * tau)*(-gamma+ exp(alpha * tau)*gamma + alpha*theta))/alpha;

	// Set up delay vector
	arma::vec delays = { tau };

	// Set up initial condition
	arma::vec IC = { Ap };

	// Set up DDE
	DDE_Example1 prob(P, Am, Ap);

	// Set up solver
	DelayRungeKutta3Solver solver(prob, IC, delays, numDelayEqs, initialTime, finalTime, stepSize, outputFileName, saveGap, printGap);

	// Solve
	solver.Solve();

	return 0;
}
