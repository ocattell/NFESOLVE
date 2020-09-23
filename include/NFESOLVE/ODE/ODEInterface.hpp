#ifndef ODEINTERFACEHEADERDEF
#define ODEINTERFACEHEADERDEF

#include <armadillo>

// Interface class for ODE problems of the type
//  du/dt = f(t,u)
//  where t is a real number
//        u is a vector of state variables

class ODEInterface
{

  public:

    // Compute right-hand side
    virtual void ComputeF(const double t, const arma::vec& u, arma::vec& F) const = 0;

    // Compute analytical solution
    virtual void ComputeAnalyticSolution(const double t, arma::vec& u) const;

};

#endif
