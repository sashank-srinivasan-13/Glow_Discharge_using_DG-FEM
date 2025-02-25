
#include  "initialize.h"

using namespace std;


void globalvar::initialization(ofstream *F1, ofstream *F2, double* in)
{
	int i,j,k;

	//ESSENTIAL CONSTANTS
	boltz_const=1.38064852*pow(10,-23);
	m_p=2.0*pow(10,-3)/(6.022*pow(10,23));
	e=1.6*pow(10,-19);
	eps_0 = 8.854*pow(10,-12);



	//PLASMA ASPECTS
	A=in[0];
	B=in[1];
	beta=in[2];
	N0=in[3];
	p=in[4];
	T_i=in[5];
	T_e=in[6];
	mu_i=in[7]/p;
	mu_e=in[8]/p;
	D_i = mu_i*T_i;
	D_e = mu_e*T_e;
	gamma_coeff=in[9];
	EMF=in[10];
	R0=in[11];
	delv=EMF;
	t_mp=eps_0/(e*N0*mu_i);

	ionization=1;

	Eold=0.0;
	Enew=0.0;

	//TRANSIENT SHEATH ASPECTS
	Lambda_D=sqrt(eps_0*boltz_const*T_e*11600.0/(N0*e*e));

	//COMPUTATIONAL ASPECTS
	l=in[12];
	dx=l/(N-1);
	pen=1.0;
	pen_bound=100000000000.0;
	Time=in[13];
	dtinput=in[14];
	t=0.0;
	itercount=0;
	outputcount=0;
	potdropcount=0;

	eps=pow(10,-10);

	potdropfile=F1;
	potentialfile=F2;


	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



	if(init_choice==2)
	{
		dt=pow(10,-24);
		for(int op=0;op<procsize;op++)
		{
			if(op==procid)
			{
				cout<<"Re initializing in process: "<<procid<<endl;
				ifstream bak2;
				ifstream bak3;
				double c,n2,dgaci2,dgace2,dgacp2,elimmax,ilimmax;
				int cell;
				bak3.open("a_bakdata.dat");
				bak3>>n2;
				bak3>>dgaci2;
				bak3>>dgace2;
				bak3>>dgacp2;
				if(n2!=N ||dgaci!=dgaci2 ||dgace!=dgace2 ||dgacp!=dgacp2)
				{
					cout<<"Preious and current runs dont match. ABORTING"<<endl;
					MPI_Abort(MPI_COMM_WORLD,1);
				}

				bak3>>ilimmax;
				bak3>>elimmax;
				if(ilimmax!=i_limiter_max || elimmax!=e_limiter_max)
					cout<<"WARNING: Limiter minimums different"<<endl;

				bak3>>start_time;
				bak3>>outputcount;
				bak3>>delv;
				t=start_time;

				bak2.open("a_bak.dat");
				for(i=0;i<N1;i++)
				{
					bak2>>cell;
					for(k=0;k<2*(dgaci+dgace+dgacp);k++)
						bak2>>c;
				}
				for(i=1;i<Nlocal;i++)
				{
					bak2>>cell;
					if(cell!=N1+i)
					{
						cout<<"Re initialization error. Reading cell: "<<c<<" instead of cell "<<N1+i<<" in process "<<procid<<endl;
						MPI_Abort(MPI_COMM_WORLD,1);
					}
					for(k=0;k<dgaci;k++)
					{
						bak2>>ni_k[i][k];
						bak2>>qi_k[i][k];
					}
					for(k=0;k<dgace;k++)
					{
						bak2>>ne_k[i][k];
						bak2>>qe_k[i][k];
					}
					for(k=0;k<dgacp;k++)
					{
						bak2>>phi_k[i][k];
						bak2>>qphi_k[i][k];
					}
				}
				bak2.close();
				bak3.close();
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		for(i=1;i<Nlocal;i++)
		{
			limflag[i]=0.0;

			ni[2*i-1]=int_ni(i,-1.0);
			ne[2*i-1]=int_ne(i,-1.0);
			qi[2*i-1]=int_qi(i,-1.0);
			qe[2*i-1]=int_qe(i,-1.0);
			phi[2*i-1]=int_phi(i,-1.0);
			qphi[2*i-1]=int_qphi(i,-1.0);

			ni[2*i]=int_ni(i,1.0);
			ne[2*i]=int_ne(i,1.0);
			qi[2*i]=int_qi(i,1.0);
			qe[2*i]=int_qe(i,1.0);
			phi[2*i]=int_phi(i,1.0);
			qphi[2*i]=int_qphi(i,1.0);


			Ef[2*i-1]=int_E(i,-1.0);
			Ef[2*i]=int_E(i,+1.0);
			alpha[i]=Alpha(i);
		}

		double points[10]={-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,-0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,-0.9739065285171717,0.9739065285171717};
		double xi;

		for(i=1;i<Nlocal;i++)
		{
			for(j=0;j<dgaci;j++)
				ni_k_RK[i][j]=ni_k[i][j];

			for(j=0;j<dgaci;j++)
				ne_k_RK[i][j]=ne_k[i][j];
		}

		updateurk();

	}
	else
	{
		dt=dtinput;
		start_time=0.0;
		for(i=0;i<nallocate;i++)
		{
			ni[2*i+1]=N0;
			ne[2*i+1]=N0;
			qi[2*i+1]=0.0;
			qe[2*i+1]=0.0;
			phi[2*i+1]=-V+V*(N1+i+1)/(N-1);
			qphi[2*i+1]=1;
			Ef[2*i+1]=0.0;

			ni[2*i]=N0;
			ne[2*i]=N0;
			qi[2*i]=0.0;
			qe[2*i]=0.0;
			phi[2*i]=-V+V*(N1+i)/(N-1);
			qphi[2*i]=1;
			Ef[2*i]=0.0;
		}

		//cell values
		for(i=0;i<nallocate;i++)
		{
			if(dgaci>=1)
			{
				ni_k[i][0]=N0;
				qi_k[i][0]=0.0;
			}

			for(j=1;j<dgaci;j++)
			{
				ni_k[i][j]=0.0;
				qi_k[i][j]=0.0;
			}


			if(dgace>=1)
			{
				ne_k[i][0]=N0;
				qe_k[i][0]=0.0;
			}

			for(j=1;j<dgace;j++)
			{
				ne_k[i][j]=0.0;
				qe_k[i][j]=0.0;
			}

			if(dgacp>=1)
			{
				phi_k[i][0]=(phi[2*i+1]+phi[2*i])/2;
				qphi_k[i][0]=(phi[2*i+1]-phi[2*i])/dx;
			}

			if(dgacp>=2)
			{
				phi_k[i][1]=(phi[2*i+1]-phi[2*i])/2;
				qphi_k[i][1]=0.0;
			}

			for(j=2;j<dgacp;j++)
			{
				phi_k[i][j]=0.0;
				qphi_k[i][j]=0.0;
			}


			Ef[2*i+1]=0.0;
			Ef[2*i]=0.0;
			alpha[i]=0.0;
			limflag[i]=0.0;
			//maxgradni[i]=abs(delv/l);
			//maxgradne[i]=abs(delv/l);

		}

	}

	for(i=0;i<10*(Nlocal-1);i++)
		for(j=0;j<10*(Nlocal-1);j++)
			AA_prototype[i][j]=0.0;

	//Initialize RK

	for(i=1;i<Nlocal;i++)
	{
		nirk[2*i-1]=ni[2*i-1];
		nirk[2*i]=ni[2*i];
		nerk[2*i-1]=ne[2*i-1];
		nerk[2*i]=ne[2*i];

		for(k=0;k<dgaci;k++)
			ni_k_RK[i][k]=ni_k[i][k];

		for(k=0;k<dgaci;k++)
			ne_k_RK[i][k]=ne_k[i][k];
	}

	for(i=0;i<nallocate;i++)
	{
		beta_i[i]=0.0;
		beta_e[i]=0.0;
	}


	rkfactor[0]=0.0;
	rkfactor[1]=1.0/4.0;
	rkfactor[2]=2.0/3.0;

	rk2factor[0]=1.0;
	rk2factor[1]=1.0/4.0;
	rk2factor[2]=2.0/3.0;

	ufactor[0]=1.0;
	ufactor[1]=3.0/4.0;
	ufactor[2]=1.0/3.0;

	iteration=10;



}







