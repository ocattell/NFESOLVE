#include "Debug.hpp"
#include "DelaySparseRungeKutta32Solver.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

DelaySparseRungeKutta32Solver::
DelaySparseRungeKutta32Solver(SparseDDEInterface& aSparseDDESystem,
    const arma::vec& initialState, const arma::vec& delays, const int numDelayEqs, const arma::umat& ZLocations,
    const double initialTime, const double finalTime, const double ATol, const double RTol,
    const std::string outputFileName,
    const int saveGap,
    const int printGap)
{
    DebugCheck(ZLocations.n_rows != 2, "DelaySparseRungeKutta32Solver(): ZLocations should be of size 2xN");
    DebugCheck(numDelayEqs > initialState.n_elem, "DelaySparseRungeKutta32Solver(): numDelayEquations is greater than the system size");
    DebugCheck(delays.min() < 1e-9, "DelaySparseRungeKutta32Solver(): Minimum delay value too small");

    // Initialise
    SetTimeInterval(initialTime, finalTime);
    mStepSize = 0.01;
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
    mpSparseDDESystem = &aSparseDDESystem;
    mpState = new arma::vec(initialState);
    mDelays = delays;
    mNumDelayEqs = numDelayEqs;
    mZLocations = ZLocations;
    mATol = ATol;
    mRTol = RTol;
    mOutputFileName = outputFileName;
    mSaveGap = saveGap;
    mPrintGap = printGap;
}

DelaySparseRungeKutta32Solver::~DelaySparseRungeKutta32Solver()
{
    delete mpState;
}

void DelaySparseRungeKutta32Solver::Solve()
{

    // Print info on screen
    std::cout << "*****************************************************************" << std::endl;
    std::cout << "Calling Delay Runge-Kutta 3(2) order sparse solver using the following parameters: " << std::endl;
    std::cout << "Initial Time " << mInitialTime << std::endl;
    std::cout << "Final Time " << mFinalTime << std::endl;
    std::cout << "Absolute Tolerance " << mATol << std::endl;
    std::cout << "Relative Tolerance " << mRTol << std::endl;
    std::cout << "Output File " << mOutputFileName << std::endl;
    std::cout << "*****************************************************************" << std::endl;

    AdaptiveDelaySparseRungeKuttaAlgorithm(ma, mb, mc, mbHat, mATol, mRTol);
}

void DelaySparseRungeKutta32Solver::SetATol(const double aTol)
{
    mATol = aTol;
}

void DelaySparseRungeKutta32Solver::SetRTol(const double rTol)
{
    mRTol = rTol;
}
