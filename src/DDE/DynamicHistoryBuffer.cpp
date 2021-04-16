#include "Debug.hpp"
#include "DynamicHistoryBuffer.hpp"
#include <chrono>
#include <thread>
#include <iostream>
#include <iomanip>

DynamicHistoryBuffer::DynamicHistoryBuffer(const int sysSize, const int bufferLength, const double maxDelay)
{
    mBufferLength = bufferLength;
    mSysSize = sysSize;
    mMaxDelay = maxDelay;

    mSolutionBuffer.zeros(mSysSize, mBufferLength);
    mDerivativeBuffer.zeros(mSysSize, mBufferLength);
    mTimeBuffer.zeros(mBufferLength);
}

void DynamicHistoryBuffer::AddEntry(const arma::vec& solution, const arma::vec& derivative, const double time)
{
	if (time - mMaxDelay < (mTimeBuffer((mColStart+1)%mBufferLength)+1e-10) && mBufferFull)
	{
		double prevStepSize = (time - mTimeBuffer((((mColStart - 1) % mBufferLength) + mBufferLength) % mBufferLength));
		double targetTime = mTimeBuffer((mColStart+1)%mBufferLength) + mMaxDelay;
		int expandAmount = ceil((targetTime - time) / prevStepSize) + 10;
		ExpandBuffer(expandAmount);
	}

	int entryIndex = (mColStart + mNumEntries) % mBufferLength;

	mSolutionBuffer.col(entryIndex) = solution;
	mDerivativeBuffer.col(entryIndex) = derivative;
	mTimeBuffer(entryIndex) = time;

	if (mBufferFull)
	{
		++mColStart;
		mColStart = mColStart % mBufferLength;
	}
	else
	{
		++mNumEntries;
		if (mNumEntries >= mBufferLength)
		{
			mBufferFull = true;
		}
	}
}

arma::vec DynamicHistoryBuffer::ComputeDelayState(double time) const
{
	arma::vec delayState(mSysSize);
	int nearestInd = arma::index_min(arma::abs(mTimeBuffer - time));

	if (fabs(mTimeBuffer(nearestInd) - time) <= 1e-10)
	{
		delayState = mSolutionBuffer.col(nearestInd);
	}
	else
	{
		int lowerInd = nearestInd;
		if (mTimeBuffer(lowerInd) == 0 && mTimeBuffer(mColStart) <= 1e-10 )
		{
			lowerInd = mColStart;
		}
		else if (mTimeBuffer(nearestInd) > time)
		{
			--lowerInd;
			lowerInd = ((lowerInd % mBufferLength) + mBufferLength) % mBufferLength;
		}
		int upperInd = (lowerInd + 1) % mBufferLength;

		double h = mTimeBuffer(upperInd) - mTimeBuffer(lowerInd);
        	DebugCheck(h<=1e-13, "DynamicHistoryBuffer: ComputeDelayState(double): Consecutive time buffer values are the same");
		double ht = time - mTimeBuffer(lowerInd);
		double theta = ht/h;

		delayState = (1.0 - theta)*mSolutionBuffer.col(lowerInd) + theta * mSolutionBuffer.col(upperInd) + theta * (theta - 1.0)*((1.0 - (2.0 * theta))*(mSolutionBuffer.col(upperInd) - mSolutionBuffer.col(lowerInd)) + (theta - 1.0)*(mTimeBuffer(upperInd) - mTimeBuffer(lowerInd))*mDerivativeBuffer.col(lowerInd) + theta * (mTimeBuffer(upperInd) - mTimeBuffer(lowerInd))*mDerivativeBuffer.col(upperInd));
	}

	return delayState;
}

double DynamicHistoryBuffer::ComputeDelayState(double time, arma::uword solIndex) const
{
	double delayState;
	int nearestInd = arma::index_min(arma::abs(mTimeBuffer - time));

	if (fabs(mTimeBuffer(nearestInd) - time) <= 1e-10)
	{
		delayState = mSolutionBuffer(solIndex, nearestInd);
	}
	else
	{
		int lowerInd = nearestInd;
		if (mTimeBuffer(lowerInd) == 0 && mTimeBuffer(mColStart) <= 1e-10)
		{
			lowerInd = mColStart;
		}
		else if (mTimeBuffer(nearestInd) > time)
		{
			--lowerInd;
			lowerInd = ((lowerInd % mBufferLength) + mBufferLength) % mBufferLength;
		}
		int upperInd = (lowerInd + 1) % mBufferLength;

		double h = mTimeBuffer(upperInd) - mTimeBuffer(lowerInd);
        	DebugCheck(h<=1e-13, "DynamicHistoryBuffer: ComputeDelayState(double, arma::uword): Consecutive time buffer values are the same");
		double ht = time - mTimeBuffer(lowerInd);
		double theta = ht / h;

		delayState = (1.0 - theta) * mSolutionBuffer(solIndex, lowerInd) + theta * mSolutionBuffer(solIndex, upperInd) + theta * (theta - 1.0) * ((1.0 - (2.0 * theta)) * (mSolutionBuffer(solIndex, upperInd) - mSolutionBuffer(solIndex, lowerInd)) + (theta - 1.0) * (mTimeBuffer(upperInd) - mTimeBuffer(lowerInd)) * mDerivativeBuffer(solIndex, lowerInd) + theta * (mTimeBuffer(upperInd) - mTimeBuffer(lowerInd)) * mDerivativeBuffer(solIndex, upperInd));
	}

	return delayState;
}

arma::vec DynamicHistoryBuffer::ComputeDelayState(double time, arma::uvec solIndices) const
{
	arma::vec delayState;
	int nearestInd = arma::index_min(arma::abs(mTimeBuffer - time));

	if (fabs(mTimeBuffer(nearestInd) - time) <= 1e-10)
	{
		delayState = mSolutionBuffer.col(nearestInd);
		delayState = delayState.elem(solIndices);
	}
	else
	{
		int lowerInd = nearestInd;
		if (mTimeBuffer(lowerInd) == 0 && mTimeBuffer(mColStart) <= 1e-10)
		{
			lowerInd = mColStart;
		}
		else if (mTimeBuffer(nearestInd) > time)
		{
			--lowerInd;
			lowerInd = ((lowerInd % mBufferLength) + mBufferLength) % mBufferLength;
		}
		int upperInd = (lowerInd + 1) % mBufferLength;

		double h = mTimeBuffer(upperInd) - mTimeBuffer(lowerInd);

        	DebugCheck(h<=1e-13, "DynamicHistoryBuffer: ComputeDelayState(double, arma::uvec): Consecutive time buffer values are the same");
		double ht = time - mTimeBuffer(lowerInd);
		double theta = ht / h;

		arma::vec lowerIndSol = mSolutionBuffer.col(lowerInd);
		arma::vec upperIndSol = mSolutionBuffer.col(upperInd);
		arma::vec lowerIndDer = mDerivativeBuffer.col(lowerInd);
		arma::vec upperIndDer = mDerivativeBuffer.col(upperInd);

		delayState = (1.0 - theta) * lowerIndSol.elem(solIndices) + theta * upperIndSol.elem(solIndices) + theta * (theta - 1.0) * ((1.0 - (2.0 * theta)) * (upperIndSol.elem(solIndices) - lowerIndSol.elem(solIndices)) + (theta - 1.0) * (mTimeBuffer(upperInd) - mTimeBuffer(lowerInd)) * lowerIndDer.elem(solIndices) + theta * (mTimeBuffer(upperInd) - mTimeBuffer(lowerInd)) * upperIndDer.elem(solIndices));
	}

	return delayState;
}

void DynamicHistoryBuffer::SetMaxDelay(const double maxDelay)
{
	mMaxDelay = maxDelay;
}

const arma::mat& DynamicHistoryBuffer::GetSolutionBuffer() const
{
  return mSolutionBuffer;
}

const arma::mat& DynamicHistoryBuffer::GetDerivativeBuffer() const
{
	return mDerivativeBuffer;
}

const arma::rowvec& DynamicHistoryBuffer::GetTimeBuffer() const
{
  return mTimeBuffer;
}

int DynamicHistoryBuffer::GetColStart() const
{
	return mColStart;
}

int DynamicHistoryBuffer::GetBufferLength() const
{
	return mBufferLength;
}

bool DynamicHistoryBuffer::IsBufferFull() const
{
	return mBufferFull;
}

void DynamicHistoryBuffer::ExpandBuffer(const int numCols)
{
    DebugCheck(numCols<=0, "DynamicHistoryBuffer: ExpandBuffer(): numCols must be positive");
	int insertIndex = mColStart;

	mSolutionBuffer.insert_cols(insertIndex, numCols);
	mDerivativeBuffer.insert_cols(insertIndex, numCols);
	mTimeBuffer.insert_cols(insertIndex, numCols);

	mBufferLength += numCols;
	mColStart += numCols;
	mBufferFull = false;
}
