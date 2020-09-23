#ifndef SPARSEDELAYMATRIXHEADERDEF
#define SPARSEDELAYMATRIXHEADERDEF

#include <armadillo>

class SparseDelayMatrix
{
public:
	SparseDelayMatrix(const arma::umat& locations, const arma::vec& values, const arma::uword numRows, const arma::uword numCols);

    double& operator()(const arma::uword row, const arma::uword col);
	const double operator()(const arma::uword row, const arma::uword col) const;
	arma::subview_col<double> Col(const arma::uword col);
	const arma::subview_col<double> Col(const arma::uword col) const;
	const arma::subview_col<arma::uword> Col_RowIndices(const arma::uword col) const;

	void SetValue(const arma::uword row, const arma::uword col, const double value, const bool checkEntryExists = true);
	void SetValues(const arma::vec& values);

	const arma::uword GetNumRows() const;
	const arma::uword GetNumCols() const;
	const arma::uword GetNumNonZeros() const;

	const arma::vec& GetValues() const;
	const arma::uvec& GetRowIndex() const;
	const arma::uvec& GetColPtr() const;

private:
	arma::uword mNumRows;
	arma::uword mNumCols;
	arma::uword mNumNonZeros;

	arma::vec mValues;
	arma::uvec mCol_ptr;
	arma::uvec mRow_index;

	bool isSorted(const arma::umat& locations);
};

#endif
