#include "poisson.h"
using namespace std;

void globalvar::poisson_solver()
{
	int i,j,k;
	double f1,f2,f3,f4,f5;
	int factor=2*dgacp;
	double P=pen_bound;
	int order=dgacp;
	int r= 2*dgacp*(Nlocal-1);
	int RTOT=2*dgacp*(N-1);
	double bc1=0.0,bc2=0.0;

	if(procid==0)
	{
		bc1=1.0;
		if(ionization)updatedelV();
		i=1;
		phi[0]=-delv;
	}
	if(procid==procsize-1)
	{
		//if(ionization)updatedelV();
		bc2=1.0;
		phi[2*Nlocal-1]=0.0;
	}

	double *L = new double[r];
	double *Ltot=new double [RTOT];

	double **nik, **nek;
	nik =new double* [Nlocal];
	nek =new double* [Nlocal];
	for(i=0;i<Nlocal;i++)
	{
		nik[i]=new double [dgacp];
		nek[i]=new double [dgacp];
	}

	for(i=0;i<Nlocal;i++)
	{
		for(k=0;k<dgacp;k++)
		{
			nik[i][k]=0.0;
			nek[i][k]=0.0;
		}
	}
	for(i=1;i<Nlocal;i++)
	{
		for(k=0;k<min(dgaci,dgacp);k++)
			nik[i][k]=ni_k_RK[i][k];

		for(k=0;k<min(dgace,dgacp);k++)
			nek[i][k]=ne_k_RK[i][k];

	}


	if(order>=1)
	{
		i=1;
		f1=(-e/eps_0*dx)*(nik[i][0]-nek[i][0]);
		L[factor*i-factor+0]=bc1*phi[0];
		L[factor*i-factor+1]=f1-bc1*P*phi[0];

		for(i=2;i<Nlocal-1;i++)
		{
			f1=(-e/eps_0*dx)*(nik[i][0]-nek[i][0]);
			L[factor*i-factor+0]=0;
			L[factor*i-factor+1]=f1;
		}

		i=Nlocal-1;

		f1=(-e/eps_0*dx)*(nik[i][0]-nek[i][0]);
		L[factor*i-factor+0]=-bc2*phi[2*i+1];
		L[factor*i-factor+1]=f1-bc2*P*phi[2*i+1];

	}

	if(order>=2)
	{
		i=1;

		f2=(-e/eps_0*dx/3)*(nik[i][1]-nek[i][1]);
		L[factor*i-factor+2]=-bc1*phi[0];
		L[factor*i-factor+3]=f2+bc1*P*phi[0];

		for(i=2;i<Nlocal-1;i++)
		{
			f2=(-e/eps_0*dx/3)*(nik[i][1]-nek[i][1]);
			L[factor*i-factor+2]=0;
			L[factor*i-factor+3]=f2;
		}

		i=Nlocal-1;

		f2=(-e/eps_0*dx/3)*(nik[i][1]-nek[i][1]);
		L[factor*i-factor+2]=-bc2*phi[2*i+1];
		L[factor*i-factor+3]=f2-bc2*P*phi[2*i+1];


	}

	if(order>=3)
	{
		i=1;

		f3=(-e/eps_0*dx/5)*(nik[i][2]-nek[i][2]);
		L[factor*i-factor+4]=bc1*phi[0];
		L[factor*i-factor+5]=f3-bc1*P*phi[0];

		for(i=2;i<Nlocal-1;i++)
		{

			f3=(-e/eps_0*dx/5)*(nik[i][2]-nek[i][2]);
			L[factor*i-factor+4]=0;
			L[factor*i-factor+5]=f3;
		}

		i=Nlocal-1;

		f3=(-e/eps_0*dx/5)*(nik[i][2]-nek[i][2]);
		L[factor*i-factor+4]=-bc2*phi[2*i+1];
		L[factor*i-factor+5]=f3-bc2*P*phi[2*i+1];

	}
	if(order>=4)
	{
		i=1;

		f4=(-e/eps_0*dx/7)*(nik[i][3]-nek[i][3]);
		L[factor*i-factor+6]=-bc1*phi[0];
		L[factor*i-factor+7]=f4+bc1*P*phi[0];

		for(i=2;i<Nlocal-1;i++)
		{
			f4=(-e/eps_0*dx/7)*(nik[i][3]-nek[i][3]);
			L[factor*i-factor+6]=0;
			L[factor*i-factor+7]=f4;
		}

		i=Nlocal-1;

		f4=(-e/eps_0*dx/7)*(nik[i][3]-nek[i][3]);
		L[factor*i-factor+6]=-bc2*phi[2*i+1];
		L[factor*i-factor+7]=f4-bc2*P*phi[2*i+1];
	}

	if(order>=5)
	{
		i=1;

		f5=(-e/eps_0*dx/9)*(nik[i][4]-nek[i][4]);
		L[factor*i-factor+8]=bc1*phi[0];
		L[factor*i-factor+9]=f5-bc1*P*phi[0];

		for(i=2;i<Nlocal-1;i++)
		{

			f5=(-e/eps_0*dx/9)*(nik[i][4]-nek[i][4]);
			L[factor*i-factor+8]=0;
			L[factor*i-factor+9]=f5;
		}

		i=Nlocal-1;

		f5=(-e/eps_0*dx/9)*(nik[i][4]-nek[i][4]);
		L[factor*i-factor+8]=-bc2*phi[2*i+1];
		L[factor*i-factor+9]=f5-bc2*P*phi[2*i+1];
	}



	//------------------------------------------------------------------------------------------------------------------------------------------------
	MPI_Status stat2;

	if(procid!=0)
	{
		MPI_Send(&r,1,MPI_INT,0,200+procid,MPI_COMM_WORLD);
		MPI_Send(L,r,MPI_DOUBLE,0,2000+procid,MPI_COMM_WORLD);
		MPI_Recv(Ltot,RTOT,MPI_DOUBLE,0,20000+procid,MPI_COMM_WORLD,&stat2);
	}

	if(procid==0)
	{
		int counttot;
		for(j=0;j<r;j++)
			Ltot[j]=L[j];
		counttot=r;
		for(i=1;i<procsize;i++)
		{
			int ralloc;
			double *Lnew;
			MPI_Recv(&ralloc,1,MPI_INT,i,200+i,MPI_COMM_WORLD,&stat2);

			Lnew=new double [ralloc];
			MPI_Recv(Lnew,ralloc,MPI_DOUBLE,i,2000+i,MPI_COMM_WORLD,&stat2);

			for(j=0;j<ralloc;j++)
				Ltot[j+counttot]=Lnew[j];
			counttot=counttot+ralloc;

			delete [] Lnew;

		}

		if(counttot!=RTOT)
			cout<<"Poisson calcs wrong"<<endl;

		for(i=1;i<procsize;i++)
			MPI_Send(Ltot,RTOT,MPI_DOUBLE,i,20000+i,MPI_COMM_WORLD);
	}


	double sum;

	for(i=0;i<r;i++)
	{
		sum=0.0;
		for(j=0;j<RTOT;j++)
			sum=sum+AA_inverse[i][j]*Ltot[j];
		L[i]=sum;
	}

	int index=0;
	for(j=1;j<Nlocal;j++)
	{
		for(k=0;k<dgacp;k++)
		{
			phi_k[j][k]=L[index];
			index++;
		}
		for(k=0;k<dgacp;k++)
		{
			qphi_k[j][k]=L[index];
			index++;
		}
		//qphi_k[j][1]=0.0;
	}


	for(i=0;i<Nlocal;i++)
	{
		delete [] nik[i];
		delete [] nek[i];
	}

	delete [] L;
	delete [] Ltot;
	delete [] nik;
	delete [] nek;


}



