#include "Debug.hpp"
#include "StaticHistoryBuffer.hpp"
#include <iomanip>

StaticHistoryBuffer::StaticHistoryBuffer(const int sysSize, const int bufferLength)
{
  mBufferLength = bufferLength;
  mSysSize = sysSize;
  mSolutionBuffer.zeros(mSysSize, mBufferLength);
  mDerivativeBuffer.zeros(mSysSize, mBufferLength);
  mTimeBuffer.zeros(mBufferLength);
}

void StaticHistoryBuffer::AddEntry(const arma::vec& solution, const arma::vec& derivative, const double time)
{
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

arma::vec StaticHistoryBuffer::ComputeDelayState(double time) const
{
	arma::vec delayState(mSysSize);
	int nearestInd = arma::index_min(arma::abs(mTimeBuffer - time));

	if (fabs(mTimeBuffer(nearestInd) - time) <= 1e-12)
	{
		delayState = mSolutionBuffer.col(nearestInd);
	}
	else
	{
        int lowerInd = nearestInd;
		if (mTimeBuffer(lowerInd) == 0 && mTimeBuffer(mColStart) <= 1e-12 )
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
        DebugCheck(h<=1e-13, "StaticHistoryBuffer: ComputeDelayState(): h too small");
		double ht = time - mTimeBuffer(lowerInd);
		double theta = ht/h;

		delayState = (1.0 - theta)*mSolutionBuffer.col(lowerInd) + theta * mSolutionBuffer.col(upperInd) + theta * (theta - 1.0)*((1.0 - (2.0 * theta))*(mSolutionBuffer.col(upperInd) - mSolutionBuffer.col(lowerInd)) + (theta - 1.0)*(mTimeBuffer(upperInd) - mTimeBuffer(lowerInd))*mDerivativeBuffer.col(lowerInd) + theta * (mTimeBuffer(upperInd) - mTimeBuffer(lowerInd))*mDerivativeBuffer.col(upperInd));
	}

	return delayState;
}

double StaticHistoryBuffer::ComputeDelayState(double time, arma::uword solIndex) const
{
	double delayState;
	int nearestInd = arma::index_min(arma::abs(mTimeBuffer - time));

	if (fabs(mTimeBuffer(nearestInd) - time) <= 1e-12)
	{
		delayState = mSolutionBuffer(solIndex, nearestInd);
	}
	else
	{
        int lowerInd = nearestInd;
		if (mTimeBuffer(lowerInd) == 0 && mTimeBuffer(mColStart) <= 1e-12)
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
        DebugCheck(h<=1e-13, "StaticHistoryBuffer: ComputeDelayState(): h too small");
		double ht = time - mTimeBuffer(lowerInd);
		double theta = ht / h;

		delayState = (1.0 - theta) * mSolutionBuffer(solIndex, lowerInd) + theta * mSolutionBuffer(solIndex, upperInd) + theta * (theta - 1.0) * ((1.0 - (2.0 * theta)) * (mSolutionBuffer(solIndex, upperInd) - mSolutionBuffer(solIndex, lowerInd)) + (theta - 1.0) * (mTimeBuffer(upperInd) - mTimeBuffer(lowerInd)) * mDerivativeBuffer(solIndex, lowerInd) + theta * (mTimeBuffer(upperInd) - mTimeBuffer(lowerInd)) * mDerivativeBuffer(solIndex, upperInd));
	}

	return delayState;
}

arma::vec StaticHistoryBuffer::ComputeDelayState(double time, arma::uvec solIndices) const
{
	arma::vec delayState;
	arma::uword nearestInd = arma::index_min(arma::abs(mTimeBuffer - time));

	if (fabs(mTimeBuffer(nearestInd) - time) <= 1e-12)
	{
        delayState = mSolutionBuffer.col(nearestInd);
		delayState = delayState.elem(solIndices);	}
	else
	{
		arma::uword lowerInd = nearestInd;
        if (mTimeBuffer(lowerInd) == 0 && mTimeBuffer(mColStart) <= 1e-12)
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
        DebugCheck(h<=1e-13, "StaticHistoryBuffer: ComputeDelayState(): h too small");
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

const arma::mat& StaticHistoryBuffer::GetSolutionBuffer() const
{
  return mSolutionBuffer;
}

const arma::mat& StaticHistoryBuffer::GetDerivativeBuffer() const
{
	return mDerivativeBuffer;
}

const arma::rowvec& StaticHistoryBuffer::GetTimeBuffer() const
{
  return mTimeBuffer;
}

int StaticHistoryBuffer::GetColStart() const
{
	return mColStart;
}

bool StaticHistoryBuffer::IsBufferFull() const
{
	return mBufferFull;
}
