#include "Debug.hpp"
#include "DelayRungeKutta3Solver.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

DelayRungeKutta3Solver::
DelayRungeKutta3Solver(DDEInterface& aDDESystem,
    const arma::vec& initialState, const arma::vec& delays, const int numDelayEqs,
    const double initialTime, const double finalTime,
    const double stepSize,
    const std::string outputFileName,
    const int saveGap,
    const int printGap)
{
    DebugCheck(numDelayEqs > initialState.n_elem, "DelayRungeKutta3Solver(): numDelayEquations is greater than the system size");
    DebugCheck(delays.min() < 1e-9, "DelayRungeKutta3Solver(): Minimum delay value too small");
    DebugCheck(stepSize < mStepSizeMin, "DelayRungeKutta3Solver(): stepSize should be greater than mStepSizeMin = " + std::to_string(mStepSizeMin));

    // Initialise
    SetTimeInterval(initialTime, finalTime);
    SetStepSize(stepSize);
    mpDDESystem = &aDDESystem;
    mpState = new arma::vec(initialState);
    mDelays = delays;
    mNumDelayEqs = numDelayEqs;
    mOutputFileName = outputFileName;
    mSaveGap = saveGap;
    mPrintGap = printGap;
}

DelayRungeKutta3Solver::~DelayRungeKutta3Solver()
{
  delete mpState;
}

void DelayRungeKutta3Solver::Solve()
{

  // Print info on screen
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Calling Delay Runge-Kutta 3rd order solver using the following parameters: " << std::endl;
  std::cout << "Initial Time " << mInitialTime << std::endl;
  std::cout << "Final Time " << mFinalTime << std::endl;
  std::cout << "Step Size " << mStepSize << std::endl;
  std::cout << "Output File " << mOutputFileName << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  DelayRungeKuttaAlgorithm(ma, mb, mb);
}
