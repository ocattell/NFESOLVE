/* SparseDelayNFE_Example1 solves a space-dependent delay neural field equation
   using an adaptive step third order sparse delay Runge-Kutta solver */

// Include relevant headers
#include "NFESOLVE.hpp"
#include "SparseDelayNeuralField.hpp"

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
    double RTol = 1e-3;
    std::string outputFileName = "SparseDelayNFE_Example1.dat";
    int saveGap = 5;
    int printGap = 5;

    // Set up NFE parameters
    double mu = 25.0; // firing rate steepness
    double theta = 0.1; // firing rate threshold
    arma::vec P = {mu, theta};

    // Set up 2D triangulated domain
    Regular2DTriGrid grid(-10, 10, -10, 10, 51);

    // Set up initial condition
    double IC_val = 0.0; // homogeneous state value
    arma::vec IC(grid.GetNumGridPoints() ,arma::fill::ones);
    IC *= IC_val;

    // Set up connectivity kernel, delays and delay matrix locations
    double tau_c = 3.0; // constant delay
    double c = 0.8; // transmission delay speed
    arma::vec delays(grid.GetNumGridPoints() * (grid.GetNumGridPoints() - 1) / 2 + 1, arma::fill::zeros);
    delays(0) = tau_c;
    arma::uword delayIndex = 1;
    arma::umat ZLocations(2, grid.GetNumGridPoints() * grid.GetNumGridPoints());
  	arma::uword ZLocationsIndex = grid.GetNumGridPoints();

    arma::mat connectivity(grid.GetNumGridPoints(), grid.GetNumGridPoints());
    for (arma::uword i=0; i<grid.GetNumGridPoints() ; ++i)
    {
        arma::uvec ZLoc = { i, 0 }; // first n entries make up the whole first column of Z
        ZLocations.col(i) = ZLoc;

        arma::rowvec x = grid.GetGridPoint(i);
        for (arma::uword j=0; j<grid.GetNumGridPoints(); ++j)
        {
            arma::rowvec xp = grid.GetGridPoint(j);
            double dist = arma::norm(x-xp);

            connectivity(i,j) = (1.0-dist)*exp(-dist);

            if (j > i)
            {
                delays(delayIndex) = tau_c + dist / c;

                arma::umat ZLocs = { {i,j}, {delayIndex, delayIndex} }; // after first n entries add in the i and jth entry for each tau_ij column
        				ZLocations.cols(ZLocationsIndex, ZLocationsIndex + 1) = ZLocs;
        				ZLocationsIndex += 2;

                ++delayIndex;
            }
        }
    }

    // Instantiate quadrature library and apply weights to connectivity matrix
    arma::vec weights;
    QuadratureLibrary::GenerateVQWeights(grid, weights);
    connectivity.each_row() %= weights.t();

    // Set up problem
    SparseDelayNeuralField NFE(P, connectivity, grid);

    // Set up solver
    DelaySparseRungeKutta32Solver solver(NFE, IC, delays, grid.GetNumGridPoints(), ZLocations, startTime, finalTime, ATol, RTol, outputFileName, saveGap, printGap);

    // Solve
    solver.Solve();

    return 0;
}
