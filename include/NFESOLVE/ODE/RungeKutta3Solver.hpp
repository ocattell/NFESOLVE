#ifndef RUNGEKUTTA3SOLVERHEADERDEF
#define RUNGEKUTTA3SOLVERHEADERDEF

#include "RungeKuttaSolver.hpp"
#include "ODEInterface.hpp"
#include <armadillo>

class RungeKutta3Solver : public RungeKuttaSolver
{
  public:

    // Constructor
    RungeKutta3Solver(ODEInterface& anODESystem,
	const arma::vec& initialState,
        const double initialTime,
	const double finalTime,
        const double stepSize,
	const std::string outputFileName="output.dat",
        const int saveGap=1,
        const int printGap=1);

    RungeKutta3Solver(ODEInterface& anODESystem,
	const arma::vec& initialState,
        const double initialTime,
	const double finalTime,
        const double stepSize,
	const std::string outputFileName,
        const int saveGap,
        const int printGap,
        const arma::uvec& outputIndices);

    // Solution
    void Solve();

    // Destructor
    ~RungeKutta3Solver();

  private:

    // Butcher tableau
    const arma::mat ma = { {0.0, 0.0, 0.0}, {1.0 / 2.0, 0.0, 0.0}, {-1.0, 2.0, 0.0} };
    const arma::vec mb = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
    const arma::vec mc = { 0.0, 1.0 / 2.0, 1.0 };
};

#endif
