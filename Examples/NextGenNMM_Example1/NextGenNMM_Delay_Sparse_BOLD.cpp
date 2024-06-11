// Include NextGenNMM
#include "NextGenNMM_Delay_Sparse.hpp"
#include <cmath>

const double pi = arma::datum::pi;

// Constructor
NextGenNMM::NextGenNMM(const arma::vec& parameters, arma::mat& connectivity, const arma::vec& constantHistory, const arma::uvec& colIndex, const arma::uvec& rowPointer)
{
	mParameters = parameters;
  	mpConnectivity = &connectivity;
	mHistoryVal = constantHistory;
	mColIndex = colIndex;
	mRowPointer = rowPointer;
}

// ComputeF
void NextGenNMM::ComputeF(const double t, const arma::vec& u, const SparseDelayMatrix& Z, arma::vec& F) const
{
		// Split up u vector to each variable
		arma::uword n = mpConnectivity->n_rows;
		arma::uword m = n*(n-1);
		
		arma::vec Re = u.subvec(0,(n-1));
		arma::vec Ve = u.subvec(n,(2*n-1));
		arma::vec Ri = u.subvec(2*n,(3*n-1));
		arma::vec Vi = u.subvec(3*n,(4*n-1));
		arma::vec gee = u.subvec(4*n,(5*n-1));
		arma::vec gei = u.subvec(5*n,(6*n-1));
		arma::vec gie = u.subvec(6*n,(7*n-1));
		arma::vec gii = u.subvec(7*n,(8*n-1));
		arma::vec gext = u.subvec(8*n,8*n+m-1);
		arma::vec hee = u.subvec(8*n+m,(9*n+m-1));
		arma::vec hei = u.subvec(9*n+m,(10*n+m-1));
		arma::vec hie = u.subvec(10*n+m,(11*n+m-1));
		arma::vec hii = u.subvec(11*n+m,(12*n+m-1));
		arma::vec hext = u.subvec(12*n+m,12*n+2*m-1);

		arma::vec BW_x = u.subvec(12*n+2*m,(13*n+2*m-1));
		arma::vec BW_f = u.subvec(13*n+2*m,(14*n+2*m-1));
		arma::vec BW_v = u.subvec(14*n+2*m,(15*n+2*m-1));
		arma::vec BW_q = u.subvec(15*n+2*m,(16*n+2*m-1));

		// Clearly define parameters
		double tau_c = mParameters(0);
		double v = mParameters(1);

		double I_E = mParameters(2);
		double tau_E = mParameters(3);
		double a_EE = mParameters(4);
		double a_EI = mParameters(5);
		double a_ext = mParameters(6);
		double vs_E = mParameters(7);
		double ks_EE = mParameters(8);
		double ks_EI = mParameters(9);
		double kv_EE = mParameters(10);
		double kv_EI = mParameters(11);
		double d_E = mParameters(12);
		double ks_ext = mParameters(13);

		double I_I = mParameters(14);
		double tau_I = mParameters(15);
		double a_IE = mParameters(16);
		double a_II = mParameters(17);
		double vs_I = mParameters(18);
		double ks_II = mParameters(19);
		double ks_IE = mParameters(20);
		double kv_II = mParameters(21);
		double kv_IE = mParameters(22);
		double d_I = mParameters(23);

		double BW_epsilon = mParameters(24);
		double BW_k = mParameters(25);
		double BW_gamma = mParameters(26);
		double BW_tau = mParameters(27);
		double BW_alpha = mParameters(28);
		double BW_rho = mParameters(29);

		double RE = mParameters(30);
		double vs_ext = mParameters(31);
		// Compute RHS

		arma::mat gmat=reshape(gext,n-1,n);
        arma::vec DRe=(1.0/tau_E)*(2.0*Re%Ve + d_E/(pi*tau_E) - Re%(gee + gei + kv_EE + kv_EI + arma::trans(arma::sum(gmat))) );
        arma::vec DVe=(1.0/tau_E)*(I_E - pow(pi,2.0)*pow(tau_E,2.0)*arma::pow(Re,2.0) + arma::pow(Ve,2.0) + kv_EI*(Vi-Ve) + gee%(vs_E-Ve) + arma::trans(arma::sum(gmat) )%(vs_ext-Ve) + gei%(vs_I-Ve) );
		arma::vec DRi = (1.0/tau_I)*(2.0*Ri%Vi + d_I/(pi*tau_I) - Ri%(gie + gii + kv_IE + kv_II));
		arma::vec DVi = (1.0/tau_I)*(I_I - pow(pi,2.0)*pow(tau_I,2.0)*arma::pow(Ri,2.0) + arma::pow(Vi,2.0) + kv_IE*(Ve-Vi) + gie%(vs_E-Vi) + gii%(vs_I-Vi));

		arma::vec Dgee = a_EE*(hee-gee);
		arma::vec Dgei = a_EI*(hei-gei);
		arma::vec Dgie = a_IE*(hie-gie);
		arma::vec Dgii = a_II*(hii-gii);
		arma::vec Dgext = a_ext*(hext-gext);

		arma::vec Dhee = a_EE*(ks_EE*Re - hee);
		arma::vec Dhei = a_EI*(ks_EI*Ri - hei);
		arma::vec Dhie = a_IE*(ks_IE*Re - hie);
		arma::vec Dhii = a_II*(ks_II*Ri - hii);
		arma::vec Dhext(m);

		if (tau_c == 0.0)
		{
			#pragma omp parallel for default(shared)
			for (arma::uword i = 0; i < Z.GetNumRows(); ++i)
			{		
				for (arma::uword j = 0; j < Z.GetNumRows(); ++j)
				{
					if ((*mpConnectivity)(i,j) == 0.0)
					{
						continue;
					}
					if (i > j)
					{
						arma::uword k = Compute_k(j,i);
						Dhext(i*(n-1)+j) =a_ext*(ks_ext * (*mpConnectivity)(i, j) * Z(j, k) - hext(i*(n-1)+j));
					}
					else
					{
						arma::uword k = Compute_k(i,j);
						Dhext(i*(n-1)+j-1) =a_ext*(ks_ext * (*mpConnectivity)(i, j) * Z(j, k) - hext(i*(n-1)+j-1));
					}
				}
			}
		}
		else
		{
			#pragma omp parallel for default(shared)
			for (arma::uword i = 0; i < Z.GetNumRows(); ++i)
			{
				double sum = 0;
				for (arma::uword j = 0; j < Z.GetNumRows(); ++j)
				{
                    			if ((*mpConnectivity)(i,j) == 0.0)
                    			{
                        			continue;
                    			}
					if (i > j)
					{
						arma::uword k = Compute_k(j,i);
						Dhext(i*(n-1)+j) =a_ext*(ks_ext * (*mpConnectivity)(i, j) * Z(j, k) - hext(i*(n-1)+j));
					}
					else
					{
						arma::uword k = Compute_k(i,j);
						Dhext(i*(n-1)+j-1) =a_ext*(ks_ext * (*mpConnectivity)(i, j) * Z(j, k) - hext(i*(n-1)+j-1));
					}
				}
			}
		}

		// BW Model
		// Z_E: (arma::pow(arma::pow(1-pi*tau_E*Re,2)+arma::pow(Ve,2),0.5)/arma::pow(arma::pow(1+pi*tau_E*Re,2)+arma::pow(Ve,2),0.5));
		arma::vec DBW_x = BW_epsilon*(Re) - BW_k*BW_x - BW_gamma*(BW_f - 1.0);
		arma::vec DBW_f = BW_x;
		arma::vec DBW_v = (1.0/BW_tau)*(BW_f - arma::pow(BW_v,(1.0/BW_alpha)));

		arma::vec temp = BW_f;
		temp.transform( [&BW_rho](double val) {return pow((1-BW_rho),(1.0/val));} );

		arma::vec DBW_q = (1.0/BW_tau)*(((BW_f/BW_rho) % (1.0 - temp)) - (BW_q % arma::pow(BW_v,(1.0/BW_alpha) - 1.0)));

		arma::vec F1 = arma::join_cols(DRe, DVe, DRi, DVi);
		arma::vec F2 = arma::join_cols(Dgee, Dgei, Dgie, Dgii);
		arma::vec F3 = arma::join_cols(Dgext, Dhee, Dhei, Dhie);
		arma::vec F4 = arma::join_cols(Dhii, Dhext);
		arma::vec F5 = arma::join_cols(F1, F2, F3, F4);

		arma::vec F6 = arma::join_cols(DBW_x, DBW_f, DBW_v, DBW_q);

		F = arma::join_cols(F5,F6);
}

// ComputeHistory
void NextGenNMM::ComputeHistory(const double t, arma::vec& history) const
{
	history = mHistoryVal;
}

// Set parameters
void NextGenNMM::SetParameters(const arma::vec& parameters)
{
	mParameters = parameters;
}

// Set connectivity matrix reference
void NextGenNMM::SetConnectivityKernel(arma::mat& connectivity)
{
	mpConnectivity = &connectivity;
}

// Compute k
arma::uword NextGenNMM::Compute_k(const arma::uword i, const arma::uword j) const
{
	for (arma::uword k = mRowPointer(i); k < mRowPointer(i+1); ++k)
	{
		if (mColIndex(k) == j)
		{
			return k;
		}
	}
	return 0;
}
