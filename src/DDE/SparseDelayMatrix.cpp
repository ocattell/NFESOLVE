#include "Debug.hpp"
#include "SparseDelayMatrix.hpp"

SparseDelayMatrix::SparseDelayMatrix(const arma::umat& locations, const arma::vec& values, const arma::uword numRows, const arma::uword numCols)
{
	DebugCheck(locations.n_rows != 2, "SparseMatrix::SparseMatrix(): Locations matrix must be 2xN");
	DebugCheck(locations.n_cols != values.n_elem, "SparseMatrix::SparseMatrix(): Locations matrix row length must be the same as values array length");
	DebugCheck(!isSorted(locations), "SparseMatrix::SparseMatrix(): Locations matrix must be sorted in column-major ordering");

	mNumRows = numRows;
	mNumCols = numCols;
	mNumNonZeros = values.n_elem;

	mValues = values;
	mRow_index.zeros(mNumNonZeros);
	mCol_ptr.zeros(numCols + 1);

	mRow_index(0) = locations(0,0);
	mCol_ptr(locations(1,0)) = 0;

	arma::uword col_ptrCounter = locations(1,0)+1;
	for (arma::uword i = 1; i < locations.n_cols; ++i)
	{
		arma::uword row_i = locations(0, i);

		arma::uword col_i = locations(1, i);
		arma::uword col_im1 = locations(1, i - 1);

		if (col_i > col_im1)
		{
			mCol_ptr(col_ptrCounter) = i;
			++col_ptrCounter;
		}
		mRow_index(i) = row_i;
	}
	mCol_ptr(col_ptrCounter) = mNumNonZeros;
}

double& SparseDelayMatrix::operator()(const arma::uword row, const arma::uword col)
{
	DebugCheck(row < 0 || col < 0 || row >= mNumRows || col >= mNumCols, "SparseMatrix::operator(): Index is out of range");

	for (arma::uword k = mCol_ptr(col); k < mCol_ptr(col + 1); ++k)
	{
		if (mRow_index(k) == row)
		{
			return mValues(k);
		}
	}
	DebugCheck(true, "SparseDelayMatrix::operator(): No element exists at index (" + std::to_string(row) + ", " + std::to_string(col) + ")");
}

const double SparseDelayMatrix::operator()(const arma::uword row, const arma::uword col) const
{
	DebugCheck(row < 0 || col < 0 || row >= mNumRows || col >= mNumCols, "SparseMatrix::operator(): Index is out of range");

	for (arma::uword k = mCol_ptr(col); k < mCol_ptr(col + 1); ++k)
	{
		if (mRow_index(k) == row)
		{
			return mValues(k);
		}
	}
	DebugCheck(true, "SparseDelayMatrix::operator(): No element exists at index (" + std::to_string(row) + ", " + std::to_string(col) + ")");
}

arma::subview_col<double> SparseDelayMatrix::Col(const arma::uword col)
{
	DebugCheck(col < 0 || col >= mNumCols, "SparseMatrix::Col(): Index is out of range");
	return mValues.subvec(mCol_ptr(col), mCol_ptr(col + 1)-1);
}

const arma::subview_col<double> SparseDelayMatrix::Col(const arma::uword col) const
{
	DebugCheck(col < 0 || col >= mNumCols, "SparseMatrix::Col(): Index is out of range");
	return mValues.subvec(mCol_ptr(col), mCol_ptr(col + 1) - 1);
}

const arma::subview_col<arma::uword> SparseDelayMatrix::Col_RowIndices(const arma::uword col) const
{
	DebugCheck(col < 0 || col >= mNumCols, "SparseMatrix::GetColRowIndices(): Index is out of range");
	return mRow_index.subvec(mCol_ptr(col), mCol_ptr(col + 1) - 1);
}

void SparseDelayMatrix::SetValue(const arma::uword row, const arma::uword col, const double value, const bool checkEntryExists)
{
	DebugCheck(row < 0 || col < 0 || row >= mNumRows || col >= mNumCols, "SparseMatrix::SetValue(): Index is out of range");
	bool entryExists = false;

	for (arma::uword k = mCol_ptr(col); k < mCol_ptr(col + 1); ++k)
	{
		if (mRow_index(k) == row)
		{
			mValues(k) = value;
			entryExists = true;
			break;
		}
	}

	if (checkEntryExists)
	{
		DebugCheck(!entryExists , "SparseMatrix::SetValue(): No element exists at index (" + std::to_string(row) + ", " + std::to_string(col) + ")");
	}
}

void SparseDelayMatrix::SetValues(const arma::vec& values)
{
	DebugCheck(values.n_elem != mValues.n_elem, "SparseMatrix::SetValues(): Incorrect values array length");
	mValues = values;
}

const arma::uword SparseDelayMatrix::GetNumRows() const
{
	return mNumRows;
}

const arma::uword SparseDelayMatrix::GetNumCols() const
{
	return mNumCols;
}

const arma::uword SparseDelayMatrix::GetNumNonZeros() const
{
	return mNumNonZeros;
}

const arma::vec& SparseDelayMatrix::GetValues() const
{
	return mValues;
}

const arma::uvec& SparseDelayMatrix::GetRowIndex() const
{
	return mRow_index;
}

const arma::uvec& SparseDelayMatrix::GetColPtr() const
{
	return mCol_ptr;
}

bool SparseDelayMatrix::isSorted(const arma::umat& locations)
{
	for (arma::uword i = 1; i < locations.n_cols; ++i)
	{
		arma::uword row_i = locations(0, i);
		arma::uword col_i = locations(1, i);
		arma::uword row_im1 = locations(0, i - 1);
		arma::uword col_im1 = locations(1, i - 1);

		if ((col_i < col_im1) || ((col_i == col_im1) && (row_i <= row_im1)))
		{
			return false;
		}
	}

	return true;
}
