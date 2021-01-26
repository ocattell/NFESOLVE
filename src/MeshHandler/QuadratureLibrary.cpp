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

    void GenerateSimpsonsRuleWeights(const Mesh& mesh, arma::vec& weights)
    {
        DebugCheck(mesh.GetDimension() != 1, "GenerateSimpsonsRuleWeights(): Quadrature rule not defined for dimension > 1");
        DebugCheck(mesh.GetNumElementNodes() != 2, "GenerateSimpsonsRuleWeights(): Quadrature rule not defined for non-interval elements");
        DebugCheck(!((mesh.GetGridPoints()).is_sorted("ascend",0)), "GenerateSimpsonsRuleWeights(): Mesh grid points must be in ascending order for this quadrature rule");

        arma::uword NumElements = mesh.GetNumElements();
        weights.zeros(mesh.GetNumGridPoints());

        bool isNOdd = 0;
        if (NumElements%2==1)
        {
            --NumElements;
            isNOdd = 1;
        }

        for (arma::uword i=0; i<=(NumElements/2)-1; ++i)
        {
            arma::vec h = {((mesh.GetGridPoint(2*i+1))(0)-(mesh.GetGridPoint(2*i))(0)), ((mesh.GetGridPoint(2*i+2))(0)-(mesh.GetGridPoint(2*i+1))(0))};
            double alpha = (2.0*pow(h(1),3)-pow(h(0),3)+3.0*h(0)*pow(h(1),2))/(6.0*h(1)*(h(1)+h(0)));
            double beta = (pow(h(1),3)+pow(h(0),3)+3.0*h(1)*h(0)*(h(1)+h(0)))/(6.0*h(1)*h(0));
            double gamma = (2.0*pow(h(0),3)-pow(h(1),3)+3.0*h(1)*pow(h(0),2))/(6.0*h(0)*(h(1)+h(0)));
            weights(2*i+2) += alpha;
            weights(2*i+1) += beta;
            weights(2*i) += gamma;
        }

        if (isNOdd)
        {
            ++NumElements;
            arma::vec h = {((mesh.GetGridPoint(NumElements-1))(0)-(mesh.GetGridPoint(NumElements-2))(0)), ((mesh.GetGridPoint(NumElements))(0)-(mesh.GetGridPoint(NumElements-1))(0))};
            double alpha = (2.0*pow(h(1),2)+3.0*h(1)*h(0))/(6.0*(h(0)+h(1)));
            double beta = (pow(h(1),2)+3.0*h(1)*h(0))/(6.0*h(0));
            double gamma = (pow(h(1),3))/(6.0*h(0)*(h(0)+h(1)));
            weights(NumElements) += alpha;
            weights(NumElements-1) += beta;
            weights(NumElements-2) += gamma;
        }
    }

    void GenerateGQPointsWeights(const Mesh& mesh, const int n, arma::mat& GQPoints, arma::vec& weights)
    {

    }
}
