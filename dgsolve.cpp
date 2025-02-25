#include "dgsolve.h"

using namespace std;

void globalvar::boundary_conditions()
{
	int i,j;

	double dummy1[2],dummy2[2],dummy3[2],dummy4[2];

	MPI_Request req;
	MPI_Status stat,stat2;

	if(procid>0)
	{

		dummy1[0]=nirk[1];
		dummy1[1]=nerk[1];
		MPI_Isend(&dummy1[0],2,MPI_DOUBLE,procid-1,100+iteration,MPI_COMM_WORLD,&req);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&dummy2[0],2,MPI_DOUBLE,procid-1,200+iteration,MPI_COMM_WORLD,&stat2);
		nirk[0]=dummy2[0];
		nerk[0]=dummy2[1];

		MPI_Wait(&req,&stat);
	}

	if(procid!=procsize-1)
	{
		dummy3[0]=nirk[2*Nlocal-2];
		dummy3[1]=nerk[2*Nlocal-2];
		MPI_Isend(&dummy3[0],2,MPI_DOUBLE,procid+1,200+iteration,MPI_COMM_WORLD,&req);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&dummy4[0],2,MPI_DOUBLE,procid+1,100+iteration,MPI_COMM_WORLD,&stat2);
		nirk[2*Nlocal-1]=dummy4[0];
		nerk[2*Nlocal-1]=dummy4[1];

		MPI_Wait(&req,&stat);

	}

	if(procid==0)
	{
		nerk[0]=(gamma_coeff*mu_i/mu_e)*nirk[1];
		nirk[0]=nirk[1];
	}

	if(procid==procsize-1)
	{
		i=Nlocal-1;
		nirk[2*i+1]=0.0;
		nerk[2*i+1]=nerk[2*i];
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------------------------


	double nip,nim,nep,nem,xi;

	for(i=1;i<Nlocal;i++)
	{

		nip=Fni(i)+Fni(i-1);
		nim=Fni(i)-Fni(i-1);

		nep=Fne(i)+Fne(i-1);
		nem=Fne(i)-Fne(i-1);



		if(dgaci>0)
			qi_k[i][0]=(1.0/dx)*nim;

		if(dgaci>1)
			qi_k[i][1]=(3.0/dx)*(nip-2.0*ni_k_RK[i][0]);

		if(dgaci>2)
			qi_k[i][2]=(5.0/dx)*(nim-2.0*ni_k_RK[i][1]);

		if(dgaci>3)
			qi_k[i][3]=(7.0/dx)*(nip-2.0*ni_k_RK[i][2]-2.0*ni_k_RK[i][0]);

		if(dgaci>4)
			qi_k[i][4]=(9.0/dx)*(nim-2.0*ni_k_RK[i][3]-2.0*ni_k_RK[i][1]);



		if(dgace>0)
			qe_k[i][0]=(1.0/dx)*nem;

		if(dgace>1)
			qe_k[i][1]=(3.0/dx)*(nep-2.0*ne_k_RK[i][0]);

		if(dgace>2)
			qe_k[i][2]=(5.0/dx)*(nem-2.0*ne_k_RK[i][1]);

		if(dgace>3)
			qe_k[i][3]=(7.0/dx)*(nep-2.0*ne_k_RK[i][2]-2.0*ne_k_RK[i][0]);

		if(dgace>4)
			qe_k[i][4]=(9.0/dx)*(nem-2.0*ne_k_RK[i][3]-2.0*ne_k_RK[i][1]);

	}

	for(i=1;i<Nlocal;i++)
	{
		qi[2*i-1]=int_qi(i,-1.0);
		qi[2*i]=int_qi(i,+1.0);

		qe[2*i-1]=int_qe(i,-1.0);
		qe[2*i]=int_qe(i,+1.0);
	}


	for(i=1;i<Nlocal;i++)
	{
		Ef[2*i-1]=int_E(i,-1.0);
		Ef[2*i]=int_E(i,+1.0);
	}

	if(positivity_on)
		reshapeq();

	//-----------------------------------------------------------------------------------------------------------------------------------------------------------
	double Dummy1[3], Dummy2[3], Dummy3[3], Dummy4[3];


	if(procid>0)
	{
		Dummy1[0]=qi[1];
		Dummy1[1]=qe[1];
		Dummy1[2]=Ef[1];

		MPI_Isend(&Dummy1[0],3,MPI_DOUBLE,procid-1,300+iteration,MPI_COMM_WORLD,&req);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&Dummy2[0],3,MPI_DOUBLE,procid-1,400+iteration,MPI_COMM_WORLD,&stat2);
		qi[0]=Dummy2[0];
		qe[0]=Dummy2[1];
		Ef[0]=Dummy2[2];

		MPI_Wait(&req,&stat);

	}

	if(procid!=procsize-1)
	{
		Dummy3[0]=qi[2*Nlocal-2];
		Dummy3[1]=qe[2*Nlocal-2];
		Dummy3[2]=Ef[2*Nlocal-2];
		MPI_Isend(&Dummy3[0],3,MPI_DOUBLE,procid+1,400+iteration,MPI_COMM_WORLD,&req);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&Dummy4[0],3,MPI_DOUBLE,procid+1,300+iteration,MPI_COMM_WORLD,&stat2);
		qi[2*Nlocal-1]=Dummy4[0];
		qe[2*Nlocal-1]=Dummy4[1];
		Ef[2*Nlocal-1]=Dummy4[2];

		MPI_Wait(&req,&stat);

	}

	if(procid==0)
		Ef[0]=Ef[1];

	if(procid==procsize-1)
		Ef[2*Nlocal-1]=Ef[2*Nlocal-2];


	//-----------------------------------------------------------------------------------------------------------------------------------------------------------


	if(procid==0)
	{
		i=1;
		qi[0]=0.0;
		qe[0]=qe[1];
		phi[0]=-delv;
	}


	if(procid==procsize-1)
	{
		i=Nlocal-1;
		qi[2*i+1]=qi[2*i];
		qe[2*i+1]=0.0;
		phi[2*i+1]=0.0;
	}


}

void globalvar::solve()
{
	MPI_Barrier(MPI_COMM_WORLD);
	int i,j,RK=-5;

	double points[13]={-1.0,+1.0,0.0,-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,-0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,-0.9739065285171717,0.9739065285171717};
	double xi;

	for(i=1;i<Nlocal;i++)
	{
		for(j=0;j<dgaci;j++)
			ni_k_RK[i][j]=ni_k[i][j];

		for(j=0;j<dgaci;j++)
			ne_k_RK[i][j]=ne_k[i][j];
	}

	updateurk();

	double maxpen=0;

	for(i=0;i<Nlocal+1;i++)
	{
		maxgradni[i]=0.0;
		for(j=0;j<13;j++)
		{
			xi=points[j];
			if(abs(int_E(i,xi))>maxgradni[i])
				maxgradni[i]=abs(int_E(i,xi));
		}
		if(maxpen<maxgradni[i])
			maxpen=maxgradni[i];
	}

	for(i=0;i<Nlocal;i++)
		maxgradne[i]=max(maxgradni[i],maxgradni[i+1]);


	for(i=0;i<Nlocal;i++)
		maxgradni[i]=maxgradne[i];


	if(convection_flux_choice==1)
	{
		if(procid>0)
		{
			MPI_Status rq;
			MPI_Send(&maxpen,1,MPI_DOUBLE,0,2311,MPI_COMM_WORLD);
			MPI_Recv(&maxpen,1,MPI_DOUBLE,0,2311,MPI_COMM_WORLD,&rq);

			for(i=1;i<Nlocal;i++)
				maxgradne[i]=maxgradni[i]=maxpen;
		}
		else
		{
			double dtpc;MPI_Status rq;
			for(int pc=1;pc<procsize;pc++)
			{
				MPI_Recv(&dtpc,1,MPI_DOUBLE,pc,2311,MPI_COMM_WORLD,&rq);
				if(dtpc>maxpen)
					maxpen=dtpc;
			}
			for(int pc=1;pc<procsize;pc++)
				MPI_Send(&maxpen,1,MPI_DOUBLE,pc,2311,MPI_COMM_WORLD);

			for(i=1;i<Nlocal;i++)
				maxgradne[i]=maxgradni[i]=maxpen;

		}

	}

	for(RK=0;RK<RKac;RK++)
	{
		iteration=RK;
		boundary_conditions();

		if(positivity_on)
		{
			//reshapeq();   //needs to be done before mpi passes
			fluxpencalc();
			calctimesteps();
		}



		for(i=1;i<Nlocal;i++)
		{

			alpha[i]=0;
			double nmuEi[5],nmuEe[5],AG[6],BNN[6];
			double Fip,Fim,Fep,Fem;
			double Dip,Dim,Dep,Dem;

			for(j=0;j<6;j++)
			{
				AG[j]=0;
				BNN[j]=0;
			}

			if(ionization==1)
			{
				int dgacc=max(dgaci,dgace);
				for(j=0;j<dgacc;j++)
				{
					AG[j]=aG(i,j);
					BNN[j]=bnn(i,j);
				}

				for(j=dgacc;j<6;j++)
				{
					AG[j]=0;
					BNN[j]=0;
				}

				alpha[i]=Alpha(i);
			}


			for(int k=1;k<dgaci;k++)
				nmuEi[k]=nmuE_i(i,k);

			for(int k=1;k<dgace;k++)
				nmuEe[k]=nmuE_e(i,k);


			Fip=CFni(i,i)+CFni(i-1,i);
			Fim=CFni(i,i)-CFni(i-1,i);
			Fep=CFne(i,i)+CFne(i-1,i);
			Fem=CFne(i,i)-CFne(i-1,i);

			Dip=Qni(i)+Qni(i-1);
			Dim=Qni(i)-Qni(i-1);
			Dep=Qne(i)+Qne(i-1);
			Dem=Qne(i)-Qne(i-1);

			if(dgaci>0)
			{
				kiRK[i][0]=(1.0/dx)*(

						-mu_i*Fim
						+D_i*Dim
						+(dx/2.0)*(AG[0])
						-(dx/2.0)*(BNN[0])
				);

			}

			if(dgace>0)
			{
				keRK[i][0]=(1.0/dx)*(

						+mu_e*Fem
						+D_e*Dem
						+(dx/2.0)*(AG[0])
						-(dx/2.0)*BNN[0]
				);
			}

			if(dgaci>1)
			{
				kiRK[i][1]= (3.0/dx)*(

						-mu_i*Fip
						+D_i*Dip
						+(nmuEi[1])
						-(2.0*D_i*qi_k[i][0])

						+(dx/2.0)*(AG[1])
						-(dx/2.0)*BNN[1]
				);
			}

			if(dgace>1)
			{
				keRK[i][1]=(3.0/dx)*(

						+mu_e*Fep
						+D_e*Dep
						-(nmuEe[1])
						-(2.0*D_e*qe_k[i][0])

						+(dx/2.0)*(AG[1])
						-(dx/2.0)*BNN[1]
				);
			}

			if(dgaci>2)
			{
				kiRK[i][2]= (5.0/dx)*(

						-mu_i*Fim
						+D_i*Dim
						+(nmuEi[2])
						-(2.0*D_i*qi_k[i][1])

						+(dx/2.0)*(AG[2])
						-(dx/2.0)*BNN[2]
				);
			}

			if(dgace>2)
			{
				keRK[i][2]=(5.0/dx)*(

						+mu_e*Fem
						+D_e*Dem
						-(nmuEe[2])
						-(2.0*D_e*qe_k[i][1])

						+(dx/2.0)*(AG[2])
						-(dx/2.0)*BNN[2]
				);
			}

			if(dgaci>3)
			{
				kiRK[i][3]= (7.0/dx)*(

						-mu_i*Fip
						+D_i*Dip
						+(nmuEi[3])
						-(2.0*D_i*qi_k[i][2]-2.0*D_i*qi_k[i][0])

						+(dx/2.0)*(AG[3])
						-(dx/2.0)*(BNN[3])
				);
			}

			if(dgace>3)
			{
				keRK[i][3]=(7.0/dx)*(

						+mu_e*Fep
						+D_e*Dep
						-(nmuEe[3])
						-(2.0*D_e*qe_k[i][2]-2.0*D_e*qe_k[i][0])

						+(dx/2.0)*(AG[3])
						-(dx/2.0)*(BNN[3])
				);
			}

			if(dgaci>4)
			{
				kiRK[i][4]= (9.0/dx)*(

						-mu_i*Fim
						+D_i*Dim
						+(nmuEi[4])
						-(2.0*D_i*qi_k[i][3]-2.0*D_i*qi_k[i][1])

						+(dx/2.0)*(AG[4])
						-(dx/2.0)*(BNN[4])
				);
			}

			if(dgace>4)
			{
				keRK[i][4]=(9.0/dx)*(

						+mu_e*Fem
						+D_e*Dem
						-(nmuEe[4])
						-(2.0*D_e*qe_k[i][3]-2.0*D_e*qe_k[i][1])

						+(dx/2.0)*(AG[4])
						-(dx/2.0)*(BNN[4])
				);
			}

			if(dgaci>5)
			{
				kiRK[i][4]= (11.0/dx)*(

						-mu_i*Fip
						+D_i*Dip
						+(nmuEi[5])
						-(2.0*D_i*qi_k[i][4]-2.0*D_i*qi_k[i][2]-2.0*D_i*qi_k[i][0])

						+(dx/2.0)*(AG[5])
						-(dx/2.0)*(BNN[5])
				);
			}

			if(dgace>5)
			{
				keRK[i][5]=(11.0/dx)*(

						+mu_e*Fep
						+D_e*Dep
						-(nmuEe[5])
						-(2.0*D_e*qe_k[i][4]-2.0*D_e*qe_k[i][2]-2.0*D_e*qe_k[i][0])

						+(dx/2.0)*(AG[5])
						-(dx/2.0)*(BNN[5])
				);
			}
		}


		updateukrk(RK);
		updateurk();

		MPI_Barrier(MPI_COMM_WORLD);
	}

	for(i=1;i<Nlocal;i++)
	{
		for(j=0;j<dgaci;j++)
			new_ni_k[i][j]=ni_k_RK[i][j];

		for(j=0;j<dgace;j++)
			new_ne_k[i][j]=ne_k_RK[i][j];
	}

}



void globalvar::updateukrk(int RK)
{
	int i,j;


	for(i=1;i<Nlocal;i++)
	{
		for(j=0;j<dgaci;j++)
		{
			new_ni_k[i][j]=
					+(ufactor[RK]*ni_k[i][j])
					+(rkfactor[RK]*ni_k_RK[i][j])
					+(rk2factor[RK]*dt*kiRK[i][j]);
		}

		for(j=0;j<dgace;j++)
		{
			new_ne_k[i][j]=
					+(ufactor[RK]*ne_k[i][j])
					+(rkfactor[RK]*ne_k_RK[i][j])
					+(rk2factor[RK]*dt*keRK[i][j]);
		}
	}


	if(positivity_on)
	{
		positivity_check();
		positivity_limiter();
	}

	slope_limiter();

	for(i=1;i<Nlocal;i++)
	{
		for(j=0;j<dgaci;j++)
			ni_k_RK[i][j]=new_ni_k[i][j];

		for(j=0;j<dgace;j++)
			ne_k_RK[i][j]=new_ne_k[i][j];
	}

}

void globalvar::updateurk()
{
	int i;
	for(i=1;i<Nlocal;i++)
	{
		nirk[2*i-1]=int_ni(i,-1.0);
		nirk[2*i]=int_ni(i,+1.0);

		nerk[2*i-1]=int_ne(i,-1.0);
		nerk[2*i]=int_ne(i,+1.0);
	}
}

