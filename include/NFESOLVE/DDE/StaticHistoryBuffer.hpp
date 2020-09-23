#ifndef STATICHISTORYBUFFERHEADERDEF
#define STATICHISTORYBUFFERHEADERDEF

#include <armadillo>

class StaticHistoryBuffer
{
public:
  // Constructor
  StaticHistoryBuffer(const int sysSize, const int bufferLength);

  void AddEntry(const arma::vec& solution, const arma::vec& derivative, const double time);

  arma::vec ComputeDelayState(double time) const;
  double ComputeDelayState(double time, arma::uword solIndex) const;
  arma::vec ComputeDelayState(double time, arma::uvec solIndices) const;

  const arma::mat& GetSolutionBuffer() const;
  const arma::mat& GetDerivativeBuffer() const;
  const arma::rowvec& GetTimeBuffer() const;
  int GetColStart() const;
  bool IsBufferFull() const;

private:
  int mBufferLength;
  int mSysSize;

  arma::mat mSolutionBuffer;
  arma::mat mDerivativeBuffer;
  arma::rowvec mTimeBuffer;

  int mColStart = 0;
  int mNumEntries = 0;
  bool mBufferFull = false;
};

#endif
