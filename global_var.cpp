
#include "global_var.h"

using namespace std;

void globalvar::allocate(int *args)
{

	N=args[0];
	dgaci=args[1];
	dgace=args[2];
	dgacp=args[3];
	i_limiter_max=args[4];
	e_limiter_max=args[5];
	limchoice=args[6];
	init_choice=args[7];
	positivity_on=args[8];
	convection_flux_choice=args[9];
	procid=args[10];
	procsize=args[11];


	if(dgacp>max(dgace,dgaci))
	{
		dgacp=max(dgaci,dgace);
		if(procid==0)
		{
			cout<<"Unwantedly high potential accuracy desired"<<endl;
			cout<<"Changing potential accuracy to: "<<max(dgaci,dgace)<<endl;
		}
	}
	RKac=3;



	Nlocal=(N-1)/procsize;
	if((N-1)%procsize!=0)
		if((N-1)%procsize>procid)
			Nlocal++;
	int i,sum=0;
	for(i=0;i<procid;i++)
	{
		int nn;
		nn=(N-1)/procsize;
		if((N-1)%procsize!=0)
			if((N-1)%procsize>i)
				nn++;
		sum=sum+nn;
	}

	N1=sum;
	N2=Nlocal+sum;
	Nlocal++; //no of points

	nallocate=Nlocal+1;  //for limiter


	ni=new double[2*nallocate];
	ne=new double[2*nallocate];
	qi=new double[2*nallocate];
	qe=new double[2*nallocate];
	phi=new double[2*nallocate];
	qphi=new double[2*nallocate];
	Ef=new double[2*nallocate];
	alpha=new double[nallocate];
	limflag=new double [nallocate];
	maxgradni=new double [nallocate];
	maxgradne=new double [nallocate];
	nirk=new double[2*nallocate];
	nerk=new double[2*nallocate];

	ni_k=new double* [nallocate];
	ne_k=new double* [nallocate];
	qi_k=new double* [nallocate];
	qe_k=new double* [nallocate];
	phi_k=new double* [nallocate];
	qphi_k=new double* [nallocate];
	new_ni_k=new double* [nallocate];
	new_ne_k=new double* [nallocate];
	ni_k_RK=new double* [nallocate];
	ne_k_RK=new double* [nallocate];
	kiRK=new double* [nallocate];
	keRK=new double* [nallocate];

	beta_i=new double[nallocate];
	beta_e=new double[nallocate];

	for(i=0;i<nallocate;i++)
	{
		ni_k[i]=new double[dgaci];
		ne_k[i]=new double[dgaci];
		qi_k[i]=new double[dgace];
		qe_k[i]=new double[dgace];
		phi_k[i]=new double[dgacp];
		qphi_k[i]=new double[dgacp];
		new_ni_k[i]=new double[dgaci];
		new_ne_k[i]=new double[dgace];
		ni_k_RK[i]=new double[dgaci];
		ne_k_RK[i]=new double[dgace];
		kiRK[i]=new double [dgaci];
		keRK[i]=new double [dgace];
	}

	AA_prototype=new double *[10*(N-1)];
	for(int i= 0; i< 10*(N-1); i++)
		AA_prototype[i] = new double[10*(N-1)];


	AA_inverse= new double *[2*dgacp*(Nlocal-1)];
	for(i=0;i<2*dgacp*(Nlocal-1);i++)
		AA_inverse[i]=new double [2*dgacp*(N-1)];



}

void globalvar::deallocate()
{

	int i;

	delete [] ni;
	delete [] ne;
	delete [] qi;
	delete [] qe;
	delete [] phi;
	delete [] qphi;
	delete [] Ef;
	delete [] alpha;
	delete [] limflag;
	delete [] maxgradni;
	delete [] maxgradne;
	delete [] nirk;
	delete [] nerk;

	for(i=0;i<nallocate;i++)
	{
		delete [] ni_k[i];
		delete [] qi_k[i];
		delete [] kiRK[i];
		delete [] new_ni_k[i];
		delete [] ni_k_RK[i];
		delete [] ne_k[i];
		delete [] qe_k[i];
		delete [] keRK[i];
		delete [] new_ne_k[i];
		delete [] ne_k_RK[i];
		delete [] phi_k[i];
		delete [] qphi_k[i];
	}
	delete [] ni_k;
	delete [] qi_k;
	delete [] kiRK;
	delete [] new_ni_k;
	delete [] ni_k_RK;
	delete [] ne_k;
	delete [] qe_k;
	delete [] keRK;
	delete [] new_ne_k;
	delete [] ne_k_RK;
	delete [] phi_k;
	delete [] qphi_k;

	delete [] beta_i;
	delete [] beta_e;


	for(i=0;i<2*dgacp*(Nlocal-1);i++)
		delete [] AA_inverse[i];

	delete [] AA_inverse;

}




double globalvar::T()
{
	return Time;
}
double globalvar::DTplasma()
{
	return timesteplim_plasma;
}
double globalvar::DTconvection()
{
	return timesteplim_convection;
}
double globalvar::DTdiffusion()
{
	return timesteplim_diffusion;
}
double globalvar::DTpospresv()
{
	return timesteplim_pospresv;
}

double globalvar::starttime()
{
	return start_time;
}
void globalvar::updatet()
{
	t=t+dt;
}

double globalvar::retV()
{
	return delv;
}
double globalvar::ret_current()
{
	return current;
}
double globalvar::ret_conduction_current()
{
	return conduction_current;
}
double globalvar::ret_total_current()
{
	return total_current;
}
double globalvar::TMP()
{
	return t_mp;
}
double globalvar::DT()
{
	return dt;
}
double globalvar::DTreset()
{
	if(ionization)
	{
		dt=dtinput;
		if(dt<pow(10,-11))
			if(t>1.5*pow(10,-5))
				dt=pow(10,-11);


		itercount++;
		if(init_choice==2)
		{
			if(itercount>10000)
			{
				dt=dtinput;
				init_choice=1;
			}
			else
				dt=pow(10,-24);
		}
	}
	else
	{
		dt=dtinput;
		itercount++;
		if(init_choice==2)
		{
			if(itercount>10000)
			{
				dt=1;	//some value so that dt goes back to linear stability		
				init_choice=1;
			}			
			else
				dt=pow(10,-24);
		}
	}

	calctimesteps(0);

	double d=min(timesteplim_convection,timesteplim_diffusion);
	d=min(d,timesteplim_plasma);
	dt=min(d,dt);

	if(dt+t>Time)
		dt=Time-t;



	if(procid>0)
		MPI_Send(&dt,1,MPI_DOUBLE,0,231,MPI_COMM_WORLD);
	else
	{
		double dtpc;MPI_Status rq;
		for(int pc=1;pc<procsize;pc++)
		{
			MPI_Recv(&dtpc,1,MPI_DOUBLE,pc,231,MPI_COMM_WORLD,&rq);
			if(dtpc!=dt)
			{
				cout<<"Timesteps dont match"<<endl;
				MPI_Abort(MPI_COMM_WORLD,23);
			}
		}
	}


	if(Time==0)
		dt=dtinput;


	return dt;
}


void globalvar::updatedelV()
{

	Eold=Enew;
	Enew=qphi_k[1][0];

	double E;
	E=(qphi_k[1][0]);

	current=(ne[1]*mu_e*e*E)*1.25*pow(10,-3);
	conduction_current=current+D_e*qe[1];
	total_current=conduction_current+eps_0*(Enew-Eold)/dt;

	delv=(EMF-R0*current);

	if(delv<0)
		delv=0.0;

/*
	if(procid==procsize-1)
	{

		Eold=Enew;
		Enew=qphi_k[Nlocal-1][0];

		double E;
		E=(qphi_k[Nlocal-1][0]);

		current=(ne[2*Nlocal-2]*mu_e*e*E)*1.25*pow(10,-3);
		conduction_current=current+D_e*qe[2*Nlocal-2];
		total_current=conduction_current+eps_0*(Enew-Eold)/dt;

		MPI_Send(&current,1,MPI_DOUBLE,0,2311,MPI_COMM_WORLD);

	}
	if(procid==0)
	{
		double dtpc;MPI_Status rq;
		MPI_Recv(&dtpc,1,MPI_DOUBLE,procsize-1,2311,MPI_COMM_WORLD,&rq);
		current=dtpc;
		delv=(EMF-R0*current);

		if(delv<0)
			delv=0.0;
	}

*/

}

void globalvar::residuals(double* R)
{
	int i,j;
	for(i=0;i<6;i++) R[i]=0.0;
	double s1,s2;
	for(i=1;i<Nlocal;i++)
	{
		s1=0.0;s2=0.0;
		for(j=0;j<dgaci;j++)
		{
			s1=s1+new_ni_k[i][j]*pow(-1.0,j);
			s2=s2+new_ni_k[i][j];
		}
		R[0]=R[0]+pow((s1-ni[2*i-1]),2.0)+pow((s2-ni[2*i]),2.0);
		s1=0.0;s2=0.0;
		for(j=0;j<dgace;j++)
		{
			s1=s1+new_ne_k[i][j]*pow(-1.0,j);
			s2=s2+new_ne_k[i][j];
		}
		R[1]=R[1]+pow((s1-ne[2*i-1]),2.0)+pow((s2-ne[2*i]),2.0);
		s1=0.0;s2=0.0;
		for(j=0;j<dgacp;j++)
		{
			s1=s1+phi_k[i][j]*pow(-1.0,j);
			s2=s2+phi_k[i][j];
		}
		R[2]=R[2]+pow((s1-phi[2*i-1]),2.0)+pow((s2-phi[2*i]),2.0);

		R[3]=R[3]+ni[2*i-1]+ni[2*i];
		R[4]=R[4]+ne[2*i-1]+ne[2*i];
		R[5]=R[5]+abs(phi[2*i-1]+phi[2*i]);
	}

}




