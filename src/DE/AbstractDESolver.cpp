#include "AbstractDESolver.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iomanip>

// Print solution on screen
void AbstractDESolver::
PrintSolution(const double time,
    const arma::vec& state,
    const bool writeHeader) const
{
    // Format output
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout.precision(6);

    // Determine size of output on screen
    int state_size = state.n_elem;
    int max_num_columns = 3;
    int num_columns = std::min(state_size, max_num_columns);
    bool trim_output = state_size > num_columns;

    // Eventually write a header
    std::string label;
    if (writeHeader)
    {
        label = "t";
        std::cout << std::setw(15) << label;
        for (int i = 0; i < num_columns; i++)
        {
            label = "u" + std::to_string(i + 1);
            std::cout << std::setw(15) << label;
        }
        if (trim_output)
        {
            std::cout << "    (displaying " << max_num_columns
                << " out of " << state_size << " components)";
        }
        std::cout << std::endl;

    }

    // Write data
    std::cout << std::setw(15) << time;
    for (int i = 0; i < num_columns; i++)
    {
        std::cout << std::setw(15) << state(i);
    }
    if (trim_output)
    {
        std::cout << std::setw(15) << "...";
    }
    std::cout << std::endl;

}

// Save solution to file
void AbstractDESolver::
SaveSolution(const double time,
    const arma::vec& state,
    const bool writeHeader,
    const std::string fileName) const
{

    // Format output
    std::ofstream output_file;
    output_file.setf(std::ios::scientific, std::ios::floatfield);
    output_file.precision(14);

    // Open file (and perform a check)
    if (writeHeader)
    {
        output_file.open(fileName);
    }
    else
    {
        output_file.open(fileName, std::ios::app);
    }
    assert(output_file.is_open());

    // Eventually write a header
    std::string label;
    if (writeHeader)
    {
        label = "t";
        output_file << "#" << std::setw(24) << label;
        for (int i = 0; i < state.n_elem; i++)
        {
            label = "u" + std::to_string(i + 1);
            output_file << std::setw(25) << label;
        }
        output_file << std::endl;

    }

    // Write data
    output_file << std::setw(25) << time;
    for (int i = 0; i < state.n_elem; i++)
    {
        output_file << std::setw(25) << state(i);
    }
    output_file << std::endl;

    output_file.close();

}

// Print Elapsed time
void AbstractDESolver::PrintElapsedTime(
    const std::chrono::high_resolution_clock::time_point t0,
    const std::chrono::high_resolution_clock::time_point t1) const
{
    std::cout << std::endl
        << "Time elapsed for time stepping = "
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
        << " milliseconds"
        << std::endl << std::endl;
}

// Accessor for initial time
double AbstractDESolver::GetInitialTime() const
{
    return mInitialTime;
}

// Accessor for final time
double AbstractDESolver::GetFinalTime() const
{
    return mFinalTime;
}

// Accessor for current state
arma::vec AbstractDESolver::GetState() const
{
    return *mpState;
}

// Set step size
void AbstractDESolver::SetStepSize(const double h)
{
    mStepSize = h;
}

// Set time interval
void AbstractDESolver::SetTimeInterval(double t0, double t1)
{
    mInitialTime = t0;
    mFinalTime = t1;
}

// Set current state
void AbstractDESolver::SetState(arma::vec& state)
{
    *mpState = state;
}

// Set Print Gap
void AbstractDESolver::SetPrintGap(const int printGap)
{
    mPrintGap = printGap;
}

// Set Save Gap
void AbstractDESolver::SetSaveGap(const int saveGap)
{
    mSaveGap = saveGap;
}

void AbstractDESolver::
SetOutputFileName(const std::string outputFileName)
{
    mOutputFileName = outputFileName;
}