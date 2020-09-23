#include "RungeKutta32Solver.hpp"

#include <iostream>

// Constructor
RungeKutta32Solver::
RungeKutta32Solver(ODEInterface& anODESystem,
    const arma::vec& initialState,
    const double initialTime, const double finalTime, const double ATol, const double RTol,
    const std::string outputFileName,
    const int saveGap,
    const int printGap)
{
  // Initialise
  SetTimeInterval(initialTime, finalTime);
  SetStepSize(0.01);
  mpODESystem = &anODESystem;
  mpState = new arma::vec(initialState);
  mATol = ATol;
  mRTol = RTol;
  mOutputFileName = outputFileName;
  mSaveGap = saveGap;
  mPrintGap = printGap;
}

// Destructor
RungeKutta32Solver::~RungeKutta32Solver()
{
  delete mpState;
}

// Solution
void RungeKutta32Solver::Solve()
{
  // Print info on screen
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Calling Runge-Kutta 3(2) order solver using the following parameters: " << std::endl;
  std::cout << "Initial Time " << mInitialTime << std::endl;
  std::cout << "Final Time " << mFinalTime << std::endl;
  std::cout << "Absolute Tolerance " << mATol << std::endl;
  std::cout << "Relative Tolerance " << mRTol << std::endl;
  std::cout << "Output File " << mOutputFileName << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  // Run algorithm
  AdaptiveRungeKuttaAlgorithm(ma, mb, mc, mbHat, mATol, mRTol);
}

// Set absolute tolerance
void RungeKutta32Solver::SetATol(const double aTol)
{
    mATol = aTol;
}

// Set relative tolerance
void RungeKutta32Solver::SetRTol(const double rTol)
{
    mRTol = rTol;
}
