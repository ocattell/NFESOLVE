#ifndef RUNGEKUTTA32SOLVERHEADERDEF
#define RUNGEKUTTA32SOLVERHEADERDEF

#include "RungeKuttaSolver.hpp"
#include "ODEInterface.hpp"
#include <armadillo>

class RungeKutta32Solver:public RungeKuttaSolver
{
  public:

    // Constructor
    RungeKutta32Solver(ODEInterface& anODESystem,
	const arma::vec& initialState,
        const double initialTime, const double finalTime, const double ATol, const double RTol,
        const std::string outputFileName="output.dat",
        const int saveGap=1,
        const int printGap=1);

    // Solution
    void Solve();

    // Set tolerances
    void SetATol(const double aTol);
    void SetRTol(const double rTol);

    // Destructor
    ~RungeKutta32Solver();

  private:

    // Butcher tableau
    const arma::mat ma = { {0.0, 0.0, 0.0}, {1.0 / 2.0, 0.0, 0.0}, {-1.0, 2.0, 0.0} };
    const arma::vec mb = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
    const arma::vec mc = { 0.0, 1.0 / 2.0, 1.0 };
    const arma::vec mbHat = { 1.0 - (1.0 / (2.0 * mc(1))), 1.0 / (2.0 * mc(1)), 0.0 };

    // Adaptation tolerances
    double mATol;
    double mRTol;
};

#endif
