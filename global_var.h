
#ifndef GLOBAL_VAR_H_
#define GLOBAL_VAR_H_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <mpi.h>


using namespace std;


class globalvar
{
	public:
	void initialization(ofstream *, ofstream *, double*);
	void allocate(int *);
	void deallocate();
	void boundary_conditions();
	void solve();
	void matrix_assemble();
	void generate_matrix_inverse();
	void poisson_solver();
	void slope_limiter();
	void update();
	void updateurk();
	void updateukrkfinal();
	void updateukrk(int RK);
	void back_up();
	void output();
	void calctimesteps(int=1);
	void mpipassu();
	void mpipassq();



	double Fni(int i);
	double Fne(int i);
	double Qni(int i);
	double Qne(int i);
	double CFni(int i, int j);
	double CFne(int i, int j);

	double T();
	double DT();
	double TMP();
	double DTreset();
	double DTplasma();
	double DTconvection();
	double DTdiffusion();
	double DTpospresv();
	double starttime();
	void updatet();
	double retV();
	double ret_current();
	double ret_conduction_current();
	double ret_total_current();
	void updatedelV();
	void residuals(double* );

	double int_P(double xi, int k);
	double int_gradP(double xi, int k);
	double int_gradE(int i, double xi);
	double int_ni(int i, double xi);
	double int_ne(int i, double xi);
	double int_qi(int i, double xi);
	double int_qe(int i, double xi);
	double int_phi(int, double);
	double int_qphi(int, double);
	double int_E(int i, double xi);
	double aG(int i, int k);
	double bnn(int i, int k);
	double nmuE_i(int i, int k);
	double nmuE_e(int i, int k);
	double Alpha(int i);
	double getalpha(int i, double xi);

	double int_newni(int i, double xi);
	double int_newne(int i, double xi);
	void fluxpencalc();
	void reshapeq();
	void positivity_check();
	void findniqk(int);
	void findneqk(int);
	void positivity_limiter();
	double pospresvtimestep();
	double reconstructqi(int,int);
	double reconstructqe(int,int);


	int N,Nlocal,N1,N2,nallocate;
	int dgaci,dgace,dgacp,RKac;
	int i_limiter_max,e_limiter_max;
	double l;
	double dx;

	int iteration, init_choice,itercount;

	int procid,procsize;
	double XL,XR;


	double *ni,*ne,*qi,*qe,*phi,*qphi;
	double **ni_k,**ne_k,**qi_k,**qe_k,**phi_k,**qphi_k;
	double **new_ni_k,**new_ne_k;

	double *nirk,*nerk;
	double **kiRK, **keRK;
	double **ni_k_RK, **ne_k_RK;


	double *Ef,*alpha;

	double **AA_prototype, **AA_inverse;
	double* limflag,limchoice;
	double pen,pen_bound;
	double start_time;
	double mu_i,mu_e,D_i,D_e,beta,e,eps_0,p,T_i,T_e;
	double N0;
	double dtinput,dt,t;
	double V,EMF,R0;
	double Lambda_D;
	double boltz_const;
	double m_p;
	double t_mp;
	double Time;
	double A,B,gamma_coeff;

	double rkfactor[3],rk2factor[3];
	double ufactor[3];

	double *maxgradni, *maxgradne;
	double *beta_i,*beta_e;

	int ionization;
	int positivity_on,convection_flux_choice;

	double timesteplim_convection,timesteplim_diffusion,timesteplim_pospresv;
	double timesteplim_plasma;

	double delv,current,conduction_current,total_current,Eold,Enew;

	double eps;

	int outputcount,potdropcount;
	ofstream *potdropfile, *potentialfile;

};







#endif /* GLOBAL_VAR_H_ */
