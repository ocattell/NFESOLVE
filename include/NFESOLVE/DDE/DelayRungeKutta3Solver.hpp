#ifndef DELAYRUNGEKUTTA3SOLVERHEADERDEF
#define DELAYRUNGEKUTTA3SOLVERHEADERDEF

#include "DelayRungeKuttaSolver.hpp"
#include "DDEInterface.hpp"
#include <armadillo>

class DelayRungeKutta3Solver : public DelayRungeKuttaSolver
{
  public:

    // Constructor
    DelayRungeKutta3Solver(DDEInterface& aDDESystem,
	const arma::vec& initialState, const arma::vec& delays, const int numDelayEqs,
        const double initialTime,
	const double finalTime,
        const double stepSize,
	const std::string outputFileName="output.dat",
        const int saveGap=1,
        const int printGap=1);

    // Constructor
    DelayRungeKutta3Solver(DDEInterface& aDDESystem,
	const arma::vec& initialState, const arma::vec& delays, const int numDelayEqs,
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
    ~DelayRungeKutta3Solver();

  private:

      // Butcher tableau
      const arma::mat ma = { {0.0, 0.0, 0.0}, {1.0 / 2.0, 0.0, 0.0}, {-1.0, 2.0, 0.0} };
      const arma::vec mb = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
      const arma::vec mc = { 0.0, 1.0 / 2.0, 1.0 };
};

#endif
