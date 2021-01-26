#include "Debug.hpp"
#include "DelaySparseRungeKutta3Solver.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <thread>

DelaySparseRungeKutta3Solver::
DelaySparseRungeKutta3Solver(SparseDDEInterface& aSparseDDESystem,
    const arma::vec& initialState, const arma::vec& delays,  const int numDelayEqs, const arma::umat& ZLocations,
    const double initialTime, const double finalTime,
    const double stepSize,
    const std::string outputFileName,
    const int saveGap,
    const int printGap)
{
    DebugCheck(ZLocations.n_rows != 2, "DelaySparseRungeKutta3Solver(): ZLocations should be of size 2xN");
    DebugCheck(numDelayEqs > initialState.n_elem, "DelaySparseRungeKutta3Solver(): numDelayEquations is greater than the system size");
    DebugCheck(delays.min() < 1e-9, "DelaySparseRungeKutta3Solver(): Minimum delay value too small");
    DebugCheck(stepSize < mStepSizeMin, "DelaySparseRungeKutta3Solver(): stepSize should be greater than mStepSizeMin = " + std::to_string(mStepSizeMin));

    // Initialise
    SetTimeInterval(initialTime, finalTime);
    mStepSize = stepSize;
    while (mStepSize > delays.min())
    {
      mStepSize /= 10.0;
    }
    if (mStepSize < mStepSizeMin)
    {
      mStepSizeMin = mStepSize/10.0;
      if (mStepSize <= 1e-8)
      {
          mStepSizeMin = 1e-9;
      }
    }
    if (mStepSize != stepSize)
    {
        std::cout << std::endl;
        std::cout << "DelaySparseRungeKutta3Solver: Step size must be smaller than min(delays)" << std::endl;
        std::cout << "Automatically adjusting step size from " << stepSize << " to " << mStepSize << std::endl;
        std::cout << std::endl;
	std::this_thread::sleep_for(std::chrono::seconds(2));
    }
    mpSparseDDESystem = &aSparseDDESystem;
    mpState = new arma::vec(initialState);
    mDelays = delays;
    mNumDelayEqs = numDelayEqs;
    mZLocations = ZLocations;
    mOutputFileName = outputFileName;
    mSaveGap = saveGap;
    mPrintGap = printGap;
}

DelaySparseRungeKutta3Solver::~DelaySparseRungeKutta3Solver()
{
  delete mpState;
}

void DelaySparseRungeKutta3Solver::Solve()
{

  // Print info on screen
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Calling Delay Runge-Kutta 3rd order sparse solver using the following parameters: " << std::endl;
  std::cout << "Initial Time " << mInitialTime << std::endl;
  std::cout << "Final Time " << mFinalTime << std::endl;
  std::cout << "Step Size " << mStepSize << std::endl;
  std::cout << "Output File " << mOutputFileName << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  DelaySparseRungeKuttaAlgorithm(ma, mb, mc);
}
