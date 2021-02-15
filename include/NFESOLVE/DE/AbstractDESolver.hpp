#ifndef ABSTRACTDESOLVERHEADERDEF
#define ABSTRACTDESOLVERHEADERDEF

#include <armadillo>
#include <chrono>
#include <string>

//! Base class for differential solvers
class AbstractDESolver
{
public:

	//! Solve initial-value problem (pure virtual)
	virtual void Solve() = 0;

	//! Print solution on screen
	void PrintSolution(const double time, const arma::vec& state, const bool writeHeader = false) const;

	//! Save solution to file
	void SaveSolution(const double time, const arma::vec& state, const bool writeHeader = false, const std::string fileName = "output.dat") const;

	//! Print elapsed time
	void PrintElapsedTime(const std::chrono::high_resolution_clock::time_point t0, const std::chrono::high_resolution_clock::time_point t1) const;

	//! Accessor functions for initial time
	double GetInitialTime() const;

	//! Accessor functions for final time
	double GetFinalTime() const;

	//! Return current solution
	arma::vec GetState() const;

	//! Set step size
	void SetStepSize(const double h);

	//! Set time interval
	void SetTimeInterval(const double t0, const double t1);

	//! Set current state
	void SetState(arma::vec& initialState);

	//! Set print gap
	void SetPrintGap(const int printGap);

	//! Set save gap
	void SetSaveGap(const int saveGap);

	//! Set file name
	void SetOutputFileName(const std::string outputFileName);

protected:

	//! Initial and finial time
	double mInitialTime;
	double mFinalTime;

	//! Current state
	arma::vec* mpState;

	//! Step size
	double mStepSize;

	//! Save every saveGap iterations
	int mSaveGap;

	//! Print every printGap iterations
	int mPrintGap;

	//! Vector of variable indices to output
	arma::uvec mOutputIndices;

	//! Output file name
	std::string mOutputFileName;
};

#endif
