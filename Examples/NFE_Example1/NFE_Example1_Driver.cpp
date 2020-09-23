//// NFE_Example1 solves a neural field equation using an adaptive step third order Runge-Kutta solver ////

// Include relevant headers
#include "NFESOLVE.hpp"
#include "NeuralField.hpp"

#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
    // Set up solver parameters
    double startTime = 0.0;
    double finalTime = 50.0;
    double ATol = 1e-6;
    double RTol = 1e-6;
    std::string outputFileName = "NFE_Example1.dat";
    int saveGap = 5;
    int printGap = 5;

    // Set up NFE parameters
    double mu = 25.0; // firing rate steepness
    double theta = 0.1; // firing rate threshold
    arma::vec P = {mu, theta};

    // Set up 1D domain
    Regular1DGrid grid(-10, 10, 101);

    // Set up initial condition
    double IC_val = 0.0; // homogeneous state value
    arma::vec IC(grid.GetNumGridPoints() ,arma::fill::ones);
    IC *= IC_val;

    // Set up connectivity kernel
    arma::mat connectivity(grid.GetNumGridPoints(), grid.GetNumGridPoints());
    for (arma::uword i=0; i<grid.GetNumGridPoints() ; ++i)
    {
        arma::rowvec x = grid.GetGridPoint(i);
        for (arma::uword j=0; j<grid.GetNumGridPoints(); ++j)
        {
            arma::rowvec xp = grid.GetGridPoint(j);
            double dist = fabs(x(0)-xp(0));

            connectivity(i,j) = (1.0-dist)*exp(-dist);
        }
    }

    // Instantiate quadrature library and apply weights to connectivity matrix
    arma::vec weights;
    QuadratureLibrary::GenerateVQWeights(grid, weights);
    connectivity.each_row() %= weights.t();

    // Set up problem
    NeuralField NFE(P, connectivity, grid);

    // Set up solver
    RungeKutta32Solver solver(NFE, IC, startTime, finalTime, ATol, RTol, outputFileName, saveGap, printGap);

    // Solve
    solver.Solve();

    return 0;
}
