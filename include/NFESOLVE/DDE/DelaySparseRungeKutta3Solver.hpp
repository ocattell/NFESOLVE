#ifndef DELAYSPARSERUNGEKUTTA3SOLVERHEADERDEF
#define DELAYSPARSERUNGEKUTTA3SOLVERHEADERDEF

#include "DelaySparseRungeKuttaSolver.hpp"
#include "SparseDDEInterface.hpp"
#include <armadillo>

class DelaySparseRungeKutta3Solver : public DelaySparseRungeKuttaSolver
{
  public:

    // Constructor
    DelaySparseRungeKutta3Solver(SparseDDEInterface& aSparseDDESystem,
	    const arma::vec& initialState, const arma::vec& delays, const int numDelayEqs, const arma::umat& ZLocations,
        const double initialTime,
	    const double finalTime,
        const double stepSize,
	    const std::string outputFileName="output.dat",
        const int saveGap=1,
        const int printGap=1);

    // Solution
    void Solve();

    // Destructor
    ~DelaySparseRungeKutta3Solver();

  private:

      // Butcher tableau
      const arma::mat ma = { {0.0, 0.0, 0.0}, {1.0 / 2.0, 0.0, 0.0}, {-1.0, 2.0, 0.0} };
      const arma::vec mb = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
      const arma::vec mc = { 0.0, 1.0 / 2.0, 1.0 };
};

#endif
