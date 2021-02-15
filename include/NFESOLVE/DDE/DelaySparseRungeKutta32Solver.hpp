#ifndef DELAYSPARSERUNGEKUTTA32SOLVERHEADERDEF
#define DELAYSPARSERUNGEKUTTA32SOLVERHEADERDEF

#include "DelaySparseRungeKuttaSolver.hpp"
#include "SparseDDEInterface.hpp"
#include <armadillo>

class DelaySparseRungeKutta32Solver : public DelaySparseRungeKuttaSolver
{
  public:

    // Constructor
    DelaySparseRungeKutta32Solver(SparseDDEInterface& aSparseDDESystem,
        const arma::vec& initialState, const arma::vec& delays, const int numDelayEqs, const arma::umat& ZLocations,
        const double initialTime,
        const double finalTime,
        const double ATol,
        const double RTol,
        const std::string outputFileName = "output.dat",
        const int saveGap = 1,
        const int printGap = 1);

    // Constructor
    DelaySparseRungeKutta32Solver(SparseDDEInterface& aSparseDDESystem,
        const arma::vec& initialState, const arma::vec& delays, const int numDelayEqs, const arma::umat& ZLocations,
        const double initialTime,
        const double finalTime,
        const double ATol,
        const double RTol,
        const std::string outputFileName,
        const int saveGap,
        const int printGap,
        const arma::uvec& outputIndices);

    // Solution
    void Solve();

    // Set tolerances
    void SetATol(const double aTol);
    void SetRTol(const double rTol);

    // Destructor
    ~DelaySparseRungeKutta32Solver();

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
