#include "DelaySparseRungeKuttaSolver.hpp"
#include "SparseDelayMatrix.hpp"
#include "StaticHistoryBuffer.hpp"
#include "DynamicHistoryBuffer.hpp"

void DelaySparseRungeKuttaSolver::DelaySparseRungeKuttaAlgorithm(const arma::mat& a, const arma::vec& b, const arma::vec& c)
{
    // Set current time
    double time = mInitialTime;

    // Set up delays
    SparseDelayMatrix delayStates(mZLocations, arma::ones<arma::vec>(mZLocations.n_cols), mNumDelayEqs, mDelays.n_elem);
    double delayTime;
    arma::vec historyTemp(mNumDelayEqs);
    for (arma::uword i = 0; i < mDelays.n_elem; ++i)
    {
        delayTime = time - mDelays(i);
        mpSparseDDESystem->ComputeHistory(delayTime, historyTemp);

        arma::uvec rowIndices = delayStates.Col_RowIndices(i);

        delayStates.Col(i) = historyTemp.elem(rowIndices);
    }

    // History buffer
    int bufferLength = ceil((mDelays.max()) / mStepSize) + 1;
    StaticHistoryBuffer historyBuffer(mNumDelayEqs, bufferLength);

    // Print and save initial condition
    if (mPrintGap > 0)
    {
        PrintSolution(time, *mpState, true);
    }
    if (mSaveGap > 0)
    {
        SaveSolution(time, *mpState, true, mOutputFileName);
    }

    // Runge kutta stages and auxiliary variables
    arma::mat k(mpState->n_elem, a.n_rows);
    arma::vec temp_sum(mpState->n_elem);
    arma::vec temp_sum2(mpState->n_elem);
    arma::vec temp_F(mpState->n_elem);

    // Step counter
    int num_steps = 0;

    // Take the current time point to monitor elapsed wall time
    std::chrono::high_resolution_clock::time_point t0
        = std::chrono::high_resolution_clock::now();

    // Main loop
    bool isFinalRun = false;
    while (mFinalTime - time >= 1e-10)
    {
        // Check if the code is on its final run
        if (fabs((time + mStepSize) - mFinalTime) <= 1e-10) isFinalRun = true;

        temp_sum2.zeros();
        // Loop through stages of RK
        for (arma::uword i = 0; i < a.n_rows; ++i)
        {
            temp_sum.zeros();
            for (arma::uword j = 0; j < i; ++j)
            {
                temp_sum += a(i, j) * k.col(j);
            }

            // Compute delayStates matrix
            #pragma omp parallel for default(shared) private(delayTime)
            for (arma::uword j = 0; j < mDelays.n_elem; ++j)
            {
                delayTime = time + mStepSize * c(i) - mDelays(j);
                arma::uvec rowIndices = delayStates.Col_RowIndices(j);

                if (delayTime < mInitialTime)
                {
                    arma::vec histTemp(mNumDelayEqs);
                    mpSparseDDESystem->ComputeHistory(delayTime, histTemp);
                    delayStates.Col(j) = histTemp.elem(rowIndices);
                }
                else
                {
                    delayStates.Col(j) = historyBuffer.ComputeDelayState(delayTime, rowIndices);
                }
            }

            mpSparseDDESystem->ComputeF(time + mStepSize * c(i), *mpState + mStepSize * temp_sum, delayStates, temp_F);

            if (i == 0)
            {
                // Add to history buffer
                historyBuffer.AddEntry((*mpState).subvec(0, mNumDelayEqs - 1), temp_F.subvec(0, mNumDelayEqs - 1), time);
            }

            k.col(i) = temp_F;

            temp_sum2 += b(i) * k.col(i);
        }

        // Update state
        *mpState = *mpState + mStepSize * temp_sum2;

        // Update time
        time += mStepSize;

        // Update number of steps counter
        num_steps++;

        // Eventually print and save new state
        if (mPrintGap > 0 && num_steps % mPrintGap == 0 || isFinalRun)
        {
            PrintSolution(time, *mpState);
        }
        if (mSaveGap > 0 && num_steps % mSaveGap == 0 || isFinalRun)
        {
            SaveSolution(time, *mpState, false, mOutputFileName);
        }

    }

    // Take the current time point to monitor elapsed wall time
    std::chrono::high_resolution_clock::time_point t1
        = std::chrono::high_resolution_clock::now();

    // Print Elapsed time
    PrintElapsedTime(t0, t1);
}

void DelaySparseRungeKuttaSolver::AdaptiveDelaySparseRungeKuttaAlgorithm(const arma::mat& a, const arma::vec& b, const arma::vec& c, const arma::vec& bHat, const double ATol, const double RTol)
{
    // Set current time
    double time = mInitialTime;

    // Set up delays
    SparseDelayMatrix delayStates(mZLocations, arma::zeros<arma::vec>(mZLocations.n_cols), mNumDelayEqs, mDelays.n_elem);
    double delayTime;
    arma::vec historyTemp(mNumDelayEqs);

    for (arma::uword i = 0; i < mDelays.n_elem; ++i)
    {
        delayTime = time - mDelays(i);
        mpSparseDDESystem->ComputeHistory(delayTime, historyTemp);

        arma::uvec rowIndices = delayStates.Col_RowIndices(i);
        delayStates.Col(i) = historyTemp.elem(rowIndices);
    }

    // History buffer
    int bufferLength = ceil((mDelays.max()) / mStepSize) + 1;
    DynamicHistoryBuffer historyBuffer(mNumDelayEqs, bufferLength, mDelays.max());

    // Print and save initial condition
    if (mPrintGap > 0)
    {
        PrintSolution(time, *mpState, true);
    }
    if (mSaveGap > 0)
    {
        SaveSolution(time, *mpState, true, mOutputFileName);
    }

    // Runge kutta stages and auxiliary variables
    arma::mat k(mpState->n_elem, a.n_rows);
    arma::vec temp_sum(mpState->n_elem);
    arma::vec temp_sum2(mpState->n_elem);
    arma::vec temp_F(mpState->n_elem);
    arma::vec temp_state(*mpState);
    arma::vec y_err(mpState->n_elem);

    // Step counter
    int num_steps = 0;

    // Take the current time point to monitor elapsed wall time
    std::chrono::high_resolution_clock::time_point t0
        = std::chrono::high_resolution_clock::now();

    // Main loop
    bool isFinalRun = false;
    while (mFinalTime - time >= 1e-10)
    {
        // Check if the code is on its final run
        if (fabs((time + mStepSize) - mFinalTime) <= 1e-10) isFinalRun = true;

        temp_sum2.zeros();
        y_err.zeros();
        // Loop through stages of RK
        for (arma::uword i = 0; i < a.n_rows; ++i)
        {
            temp_sum.zeros();
            for (arma::uword j = 0; j < i; ++j)
            {
                temp_sum += a(i, j) * k.col(j);
            }

            // Compute delayStates matrix
            #pragma omp parallel for default(shared) private(delayTime)
            for (arma::uword j = 0; j < mDelays.n_elem; ++j)
            {
                delayTime = time + mStepSize * c(i) - mDelays(j);
                arma::uvec rowIndices = delayStates.Col_RowIndices(j);

                if (delayTime < mInitialTime)
                {
                    arma::vec histTemp(mNumDelayEqs);
                    mpSparseDDESystem->ComputeHistory(delayTime, histTemp);

                    delayStates.Col(j) = histTemp.elem(rowIndices);
                }
                else
                {
                    delayStates.Col(j) = historyBuffer.ComputeDelayState(delayTime, rowIndices);
                }
            }

            mpSparseDDESystem->ComputeF(time + mStepSize * c(i), *mpState + mStepSize * temp_sum, delayStates, temp_F);

            if (i == 0)
            {
                // Add to history buffer
                historyBuffer.AddEntry((*mpState).subvec(0, mNumDelayEqs - 1), temp_F.subvec(0, mNumDelayEqs - 1), time);
            }

            k.col(i) = temp_F;

            temp_sum2 += b(i) * k.col(i);
            y_err += mStepSize * (b(i) - bHat(i)) * k.col(i);
        }
        // Update state
        temp_state = *mpState;
        *mpState = *mpState + mStepSize * temp_sum2;

        // Update time
        time += mStepSize;

        // Step size adjustment
        arma::vec sc = ATol + arma::max(arma::abs(temp_state), arma::abs(*mpState)) * RTol;
        double temp_sum3 = 0;
        for (arma::uword i = 0; i < mpState->n_elem; ++i)
        {
            temp_sum3 += pow((y_err(i) / sc(i)), 2.0);
        }

        double err = sqrt((1.0 / mpState->n_elem) * temp_sum3);
        mStepSize = std::max(std::min(mDelays.min(), mStepSize * std::min(mFacMax, std::max(mFacMin, mFac * pow(1.0 / err, 1.0 / a.n_rows)))), mStepSizeMin);

        if (err > 1.0)
        {
            mStepSize = std::max(std::min(mDelays.min(), mStepSize * std::min(1.0, std::max(mFacMin, mFac * pow(1.0 / err, 1.0 / a.n_rows)))), mStepSizeMin);
        }

        if (time + mStepSize > mFinalTime)
        {
            mStepSize = mFinalTime - time;
        }

        // Update number of steps counter
        num_steps++;

        // Eventually print and save new state
        if (mPrintGap > 0 && num_steps % mPrintGap == 0 || isFinalRun)
        {
            PrintSolution(time, *mpState);
        }
        if (mSaveGap > 0 && num_steps % mSaveGap == 0 || isFinalRun)
        {
            SaveSolution(time, *mpState, false, mOutputFileName);
        }

    }

    // Take the current time point to monitor elapsed wall time
    std::chrono::high_resolution_clock::time_point t1
        = std::chrono::high_resolution_clock::now();

    // Print Elapsed time
    PrintElapsedTime(t0, t1);
}

void DelaySparseRungeKuttaSolver::SetStepSizeMin(const double stepSizeMin)
{
    mStepSizeMin = stepSizeMin;
}

void DelaySparseRungeKuttaSolver::SetFacMin(const double facMin)
{
    mFacMin = facMin;
}

void DelaySparseRungeKuttaSolver::SetFac(const double fac)
{
    mFac = fac;
}

void DelaySparseRungeKuttaSolver::SetFacMax(const double facMax)
{
    mFacMax = facMax;
}
