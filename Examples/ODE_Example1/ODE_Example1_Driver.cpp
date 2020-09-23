//// ODE_Example1 solves a system of 3 ODEs using a fixed step third order Runge-Kutta solver ////

// Include relevant headers
#include "NFESOLVE.hpp"
#include "ODE_Example1.hpp"

#include <armadillo>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	// Set up solver parameters
	double initialTime = 0.0;
	double finalTime = 10;
	double stepSize = 0.01;
	std::string outputFileName = "ODE_Example1.dat";
	int saveGap = 1;
	int printGap = 1;

	// Set up ODE parameters
	double a = 0.5;
	double b = 1.0;
	arma::vec P = {a, b};

	// Set up initial condition
	arma::vec IC = {1, 2, 1};

	// Set up ODE
	ODE_Example1 prob(P);

	// Set up solver
	RungeKutta3Solver solver(prob, IC, initialTime, finalTime, stepSize, outputFileName, saveGap, printGap);
	
	// Solve
	solver.Solve();

	return 0;
}
