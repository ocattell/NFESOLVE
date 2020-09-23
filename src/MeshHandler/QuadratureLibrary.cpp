#include "Debug.hpp"
#include "QuadratureLibrary.hpp"
#include <cmath>

namespace QuadratureLibrary
{
    void GenerateVQWeights(const Mesh& mesh, arma::vec& weights)
    {
        arma::uword NumElements = mesh.GetNumElements();
        arma::uword dim = mesh.GetDimension();
        weights.zeros(mesh.GetNumGridPoints());
        if (dim == 1)
        {
            for (arma::uword k=0; k<NumElements; ++k)
            {
                arma::urowvec elementConnectivity = mesh.GetElementConnectivity(k);
                arma::mat elementGridPoints = mesh.GetElementGridPoints(k);
                double deltaX = fabs(elementGridPoints(1,0) - elementGridPoints(0,0));
                for (arma::uword i=0; i<mesh.GetNumElementNodes(); ++i)
                {
                    weights(elementConnectivity(i)) += deltaX/2.0;
                }
            }
        }
        else if (dim == 2)
        {
            if (mesh.GetNumElementNodes()==3)
            {
                for (arma::uword k=0; k<NumElements; ++k)
                {
                    arma::urowvec elementConnectivity = mesh.GetElementConnectivity(k);
                    arma::mat elementGridPoints = mesh.GetElementGridPoints(k);
                    double A_k = fabs((elementGridPoints(0,0)*(elementGridPoints(1,1) - elementGridPoints(2,1)) +
                        elementGridPoints(1,0)*(elementGridPoints(2,1) - elementGridPoints(0,1)) +
                        elementGridPoints(2,0)*(elementGridPoints(0,1) - elementGridPoints(1,1)))/2.0);
                    for (arma::uword i=0; i<mesh.GetNumElementNodes(); ++i)
                    {
                        weights(elementConnectivity(i)) += 2*A_k*(1.0/6.0);
                    }
                }
            }
            else if (mesh.GetNumElementNodes()==4)
            {

            }
            else
            {
                DebugCheck(true, "GenerateVQWeights(): Quadrature rule not defined for this element type");
            }
        }
        else if (dim == 3)
        {
            if (mesh.GetNumElementNodes()==3)
            {
                for (arma::uword k=0; k<NumElements; ++k)
                {
                    arma::urowvec elementConnectivity = mesh.GetElementConnectivity(k);
                    arma::mat elementGridPoints = mesh.GetElementGridPoints(k);
                    double A_k = 0.5*sqrt(pow(elementGridPoints(0,1)*(elementGridPoints(1,0) - elementGridPoints(2,0)) + elementGridPoints(1,1)*(elementGridPoints(2,0) - elementGridPoints(0,0)) + elementGridPoints(2,1)*(elementGridPoints(0,0) - elementGridPoints(1,0)),2.0) +
                        pow(elementGridPoints(0,2)*(elementGridPoints(1,0) - elementGridPoints(2,0)) + elementGridPoints(1,2)*(elementGridPoints(2,0) - elementGridPoints(0,0)) + elementGridPoints(2,2)*(elementGridPoints(0,0) - elementGridPoints(1,0)),2.0) +
                        pow(elementGridPoints(0,2)*(elementGridPoints(1,1) - elementGridPoints(2,1)) + elementGridPoints(1,2)*(elementGridPoints(2,1) - elementGridPoints(0,1)) + elementGridPoints(2,2)*(elementGridPoints(0,1) - elementGridPoints(1,1)),2.0));
                    for (arma::uword i=0; i<mesh.GetNumElementNodes(); ++i)
                    {
                        weights(elementConnectivity(i)) += 2*A_k*(1.0/6.0);
                    }
                }
            }
            else if (mesh.GetNumElementNodes()==4)
            {

            }
            else
            {
                DebugCheck(true, "GenerateVQWeights(): Quadrature rule not defined for this element type");
            }
        }
        else
        {
            DebugCheck(true, "GenerateVQWeights(): Quadature rule not defined for dimension > 3");
        }
    }

    void GenerateVQWeights(const Regular2DQuadGrid& mesh, arma::vec& weights)
    {
        arma::uword NumElements = mesh.GetNumElements();
        weights.zeros(mesh.GetNumGridPoints());
        for (arma::uword k=0; k<NumElements; ++k)
        {
            arma::urowvec elementConnectivity = mesh.GetElementConnectivity(k);
            arma::mat elementGridPoints = mesh.GetElementGridPoints(k);
            double A_k = arma::norm(elementGridPoints.row(0) - elementGridPoints.row(1)) * arma::norm(elementGridPoints.row(1) - elementGridPoints.row(2));
            for (arma::uword i=0; i<mesh.GetNumElementNodes(); ++i)
            {
                weights(elementConnectivity(i)) += A_k/4.0;
            }
        }
    }

    void GenerateGQPointsWeights(const Mesh& mesh, const int n, arma::mat& GQPoints, arma::vec& weights)
    {

    }
}
