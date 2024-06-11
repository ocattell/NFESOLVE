//// NottinghamModel_Delay implements and solves the Nottingham Neural Field Model ////

// Include relevant headers
#include "NFESOLVE.hpp"
#include "NextGenNMM_Delay_Sparse.hpp"

#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>


const double pi = arma::datum::pi;

int main(int argc, char* argv[])
{
    std::string atlas = "DK68";//"MMP360" or "AAL78"
    std::string subject = "Caverage";
    std::string normState = "logrownorm";
    //std::string pSet = subject + "_" + normState + "_P0";
    std::string pSet = "Caverage_logrownorm_P0";

    // Set up solver parameters
    double startTime = 0.0;
    double finalTime = 300.0;
    double ATol = 1e-7;
    double RTol = 1e-4;
    std::string outputFileName = "output_gaps.dat";
    int saveGap = 10;
    int printGap = 10;

    // Read NFE parameters
    arma::rowvec P_in;
    P_in.load("params_gaps.dat", arma::raw_ascii);
    double I_E =P_in[0];
    double tau_E =P_in[1];
    double a_EE = P_in[2];
    double a_EI = P_in[3];
    double a_ext = P_in[4];
    double vs_E = P_in[5];
    double ks_EE = P_in[6];
    double ks_EI = P_in[7];
    double kv_EE = P_in[8];
    double d_E = P_in[9];
    double ks_ext = P_in[10];
    double vs_ext = P_in[11];

    double I_I = P_in[12];
    double tau_I = P_in[13];
    double a_IE = P_in[14];
    double a_II = P_in[15];
    double vs_I = P_in[16];
    double ks_II = P_in[17];
    double ks_IE = P_in[18];
    double kv_II = P_in[19];
    double d_I = P_in[20];
    double tau_c = P_in[21];
    double v = P_in[22];

    // BW parameters
    double BW_epsilon = P_in[23];
    double BW_k = P_in[24];
    double BW_gamma = P_in[25];
    double BW_tau = P_in[26];
    double BW_alpha = P_in[27];
    double BW_rho = P_in[28];
    // Load connectivity matrix
    arma::mat W;
    W.load("C68.dat");
    arma::uword N = W.n_rows;
    arma::uvec M = find(W);
    arma::uword m = M.n_rows;
    arma::cx_vec evals;
    arma::cx_mat evecs;

    arma::eig_gen(evals,evecs,W);
    arma::uvec indices = arma::sort_index(arma::real(evals),"descend");

    // Choose eigenmode to perturb IC by
    arma::uword emode_num = 7;
    arma::vec emode = arma::real(evecs.col(indices(emode_num - 1)));

    // Load distance matrix
    arma::mat D;
    D.load("D68.txt");
    // arma::regspace<arma::uvec>(14*N+2*m,16*N+2*m-1) <- BOLD vars
    arma::uvec outputIndices = arma::join_cols(arma::regspace<arma::uvec>(0,N-1), arma::regspace<arma::uvec>(N,2*N-1), arma::regspace<arma::uvec>(14*N+2*m,15*N+2*m-1), arma::regspace<arma::uvec>(15*N+2*m,16*N+2*m-1) );

    // Load initial conditions
    arma::rowvec IC0;
    IC0.load("initial_conditions_gaps.dat", arma::raw_ascii);
    double RE=IC0[0], RI=IC0[1], VE=IC0[2], VI=IC0[3];
    arma::vec IC(16*N+2*m, arma::fill::ones);
    IC.subvec(0, N - 1) *= RE;
    IC.subvec(N, 2*N - 1) *= VE;
    IC.subvec(2*N, 3*N - 1) *= RI;
    IC.subvec(3*N, 4*N - 1) *= VI;
    IC.subvec(4*N, 5*N - 1) *= (ks_EE*RE);
    IC.subvec(5*N, 6*N - 1) *= (ks_EI*RI);
    IC.subvec(6*N, 7*N - 1) *= (ks_IE*RE);
    IC.subvec(7*N, 8*N - 1) *= (ks_II*RI);
    IC.subvec(8*N,8*N+m-1) = ks_ext*nonzeros(W)*RE;
    IC.subvec(8*N+m, 9*N+m - 1) *= 0.0;
    IC.subvec(9*N+m, 10*N+m - 1) *= 0.0;
    IC.subvec(10*N+m, 11*N+m - 1) *= 0.0;
    IC.subvec(11*N+m, 12*N+m - 1) *= 0.0;
    IC.subvec(12*N+m, 13*N+2*m - 1) *= 0.0;

    double eScale = 0.1;
    IC.subvec(0,N-1) += eScale*emode;
    IC.subvec(N, 2*N - 1) += eScale*emode;
    IC.subvec(2*N, 3*N - 1) += eScale*emode;
    IC.subvec(3*N, 4*N - 1) += eScale*emode;
    IC.subvec(4*N, 5*N - 1) += eScale*emode;
    IC.subvec(5*N, 6*N - 1) += eScale*emode;
    IC.subvec(6*N, 7*N - 1) += eScale*emode;
    IC.subvec(7*N, 8*N - 1) += eScale*emode;
    IC.subvec(8*N+m, 9*N+m - 1) += eScale*emode;
    IC.subvec(9*N+m, 10*N+m - 1) += eScale*emode;
    IC.subvec(10*N+m, 11*N+m - 1) += eScale*emode;
    IC.subvec(11*N+m, 12*N+m - 1) += eScale*emode;
    arma::vec P =
        {tau_c,
        v,
        I_E,
        tau_E,
        a_EE,
        a_EI,
        a_ext,
        vs_E,
        ks_EE,
        ks_EI,
        kv_EE,
        0.0, //sqrt(kv_EE)*sqrt(kv_II),
        d_E,
        ks_ext,
        I_I,
        tau_I,
        a_IE,
        a_II,
        vs_I,
        ks_II,
        ks_IE,
        kv_II,
        0.0, //sqrt(kv_EE)*sqrt(kv_II),
        d_I,
	BW_epsilon,
	BW_k,
	BW_gamma,
	BW_tau,
	BW_alpha,
	BW_rho,
	RE,
	vs_ext
        };
    // Set up delays
    arma::uword numDelays = arma::accu(arma::trimatu(W,1) > 0.0);
    arma::uword numDelayValues = numDelays*2;
    if (tau_c != 0.0)
    {
        numDelays += 1;
    }
    arma::vec delays(numDelays);
    arma::uword delayIndex = 0;
    arma::uword ZLocationsIndex = 0;
    if (tau_c != 0.0)
    {
      delays(0) = tau_c;
      delayIndex = 1;
      numDelayValues += N;
      ZLocationsIndex = N;
    }
    arma::umat ZLocations(2, numDelayValues);

    // Set up delays
    arma::uvec colIndex(numDelays);
    arma::uvec rowPointer(N);
    arma::uword colCounter = 0;
    for (arma::uword i=0; i<N ; ++i)
    {
        if (tau_c != 0.0)
        {
          arma::uvec ZLoc = { i, 0 };
          ZLocations.col(i) = ZLoc;
        }

	rowPointer(i) = colCounter;
        for (arma::uword j=0; j<N; ++j)
        {
		if (W(i,j) == 0.0 || D(i,j) == 0.0)
		{
			continue;
		}
            	if (j > i)
      		{
			delays(delayIndex) = tau_c + D(i,j)/v;

      			arma::umat ZLocs = { {i,j}, {delayIndex, delayIndex} };
      			ZLocations.cols(ZLocationsIndex, ZLocationsIndex + 1) = ZLocs;

      			ZLocationsIndex += 2;
      			++delayIndex;

			colIndex(colCounter) = j;
			++colCounter;
      		}
        }
    }

    // Set up problem
    NextGenNMM NFE(P, W, IC.subvec(0,N-1), colIndex, rowPointer); 
    // Set up solver
    DelaySparseRungeKutta32Solver solver(NFE, IC, delays, N, ZLocations, startTime, finalTime, ATol, RTol, outputFileName, saveGap, printGap, outputIndices);
    // Solve
    solver.Solve();

    return 0;
}
