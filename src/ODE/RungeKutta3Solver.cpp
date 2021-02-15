#include "Debug.hpp"
#include "RungeKutta3Solver.hpp"

#include <iostream>

// Constructor
RungeKutta3Solver::
RungeKutta3Solver(ODEInterface& anODESystem,
    const arma::vec& initialState,
    const double initialTime, const double finalTime,
    const double stepSize,
    const std::string outputFileName,
    const int saveGap,
    const int printGap)
{
    DebugCheck(stepSize < mStepSizeMin, "RungeKutta3Solver(): stepSize should be greater than mStepSizeMin = " + std::to_string(mStepSizeMin));

    // Initialise
    SetTimeInterval(initialTime, finalTime);
    SetStepSize(stepSize);
    mpODESystem = &anODESystem;
    mpState = new arma::vec(initialState);
    mOutputFileName = outputFileName;
    mSaveGap = saveGap;
    mPrintGap = printGap;
    mOutputIndices = arma::regspace<arma::uvec>(0,initialState.n_elem - 1);
}

// Constructor
RungeKutta3Solver::
RungeKutta3Solver(ODEInterface& anODESystem,
    const arma::vec& initialState,
    const double initialTime, const double finalTime,
    const double stepSize,
    const std::string outputFileName,
    const int saveGap,
    const int printGap,
    const arma::uvec& outputIndices)
{
    DebugCheck(stepSize < mStepSizeMin, "RungeKutta3Solver(): stepSize should be greater than mStepSizeMin = " + std::to_string(mStepSizeMin));
    DebugCheck(outputIndices.max() >= initialState.n_elem, "RungeKutta3Solver(): outputIndices value out of bounds. It should only contain values between 0 and " + std::to_string(initialState.n_elem - 1));
    DebugCheck(outputIndices.is_empty(), "RungeKutta3Solver(): outputIndices must be non-empty");
    DebugCheck(arma::any(outputIndices != arma::unique(outputIndices)), "RungeKutta3Solver(): outputIndices must only contain unique elements");

    // Initialise
    SetTimeInterval(initialTime, finalTime);
    SetStepSize(stepSize);
    mpODESystem = &anODESystem;
    mpState = new arma::vec(initialState);
    mOutputFileName = outputFileName;
    mSaveGap = saveGap;
    mPrintGap = printGap;
    mOutputIndices = outputIndices;
}

// Destructor
RungeKutta3Solver::~RungeKutta3Solver()
{
  delete mpState;
}

// Solution
void RungeKutta3Solver::Solve()
{

  // Print info on screen
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Calling Runge-Kutta 3rd order solver using the following parameters: " << std::endl;
  std::cout << "Initial Time " << mInitialTime << std::endl;
  std::cout << "Final Time " << mFinalTime << std::endl;
  std::cout << "Step Size " << mStepSize << std::endl;
  std::cout << "Output File " << mOutputFileName << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  // Run algorithm
  RungeKuttaAlgorithm(ma, mb, mc);
}
