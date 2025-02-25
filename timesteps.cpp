

#include "timesteps.h"

using namespace std;

void globalvar::calctimesteps(int a)
{
	int i,j;
	double points[10]={-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,-0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,-0.9739065285171717,0.9739065285171717};
	//----------------------------------------------------------------------------------------------------------------------------------------
	//timetsep plasma


	double dtnew;
	double dielectric_relaxation_time, courant_limit;
	double xi,max_E=0,max_ne=0,max_ni=0;
	double maxd;


	for(i=1;i<Nlocal;i++)
	{
		for(j=0;j<10;j++)
		{
			xi=points[j];
			if(max_E<abs(int_E(i,xi)))
				max_E=abs(int_E(i,xi));

			if(max_ni<int_ni(i,xi))
				max_ni=int_ni(i,xi);

			if(max_ne<int_ne(i,xi))
				max_ne=int_ne(i,xi);
		}
	}
	courant_limit=dx*dx/(2.0*D_e+mu_e*dx*max_E);
	dielectric_relaxation_time=eps_0/(e*max_ne*mu_e)*pow(10,-2);

	timesteplim_plasma=min(courant_limit,dielectric_relaxation_time);


	
	//----------------------------------------------------------------------------------------------------------------------------------------

	//timestep convection
	maxd=maxgradni[2];

	if(convection_flux_choice==2)
	{
		for(i=1;i<Nlocal;i++)
			if(maxd<maxgradni[i])
				maxd=maxgradni[i];
	}


	timesteplim_convection=dx/(maxd*(2.0*max(dgaci,dgace)-1));
	if(t==0)
		timesteplim_convection=1000;
	//----------------------------------------------------------------------------------------------------------------------------------------
	//timestep diffusion


	double cci,cce;
	if(dgaci==2)
		cci=0.05;
	if(dgaci==3)
		cci=0.01;
	if(dgaci==4)
		cci=0.005;
	if(dgaci==5)
		cci=0.0002;
	if(dgaci==6)
		cci=0.0001;

	if(dgace==2)
		cce=0.05;
	if(dgace==3)
		cce=0.01;
	if(dgace==4)
		cce=0.005;
	if(dgace==5)
		cce=0.0002;
	if(dgace==6)
		cce=0.0001;

	timesteplim_diffusion=dx*dx*min(cci/D_i,cce/D_e);



	//----------------------------------------------------------------------------------------------------------------------------------------
	//timestep pos presv
	if(positivity_on)
	{
		if(a)
			timesteplim_pospresv=pospresvtimestep();
	}
	else
		timesteplim_pospresv=1000-procid;

	//----------------------------------------------------------------------------------------------------------------------------------------
	//final timestep

	double ts[4];
	ts[0]=timesteplim_plasma;
	ts[1]=timesteplim_convection;
	ts[2]=timesteplim_diffusion;
	ts[3]=timesteplim_pospresv;

	MPI_Status stat,stat2;
	MPI_Request req;

	if(procid>0)
	{
		MPI_Send(ts,4,MPI_DOUBLE,0,650,MPI_COMM_WORLD);
		MPI_Recv(ts,4,MPI_DOUBLE,0,698,MPI_COMM_WORLD,&stat2);
	}
	if(procid==0)
	{
		double tss[4];
		for(int pid=1;pid<procsize;pid++)
		{
			MPI_Recv(tss,4,MPI_DOUBLE,pid,650,MPI_COMM_WORLD,&stat2);
			for(int g=0;g<4;g++)
				if(tss[g]<ts[g])
					ts[g]=tss[g];
		}

		for(int pid=1;pid<procsize;pid++)
			MPI_Send(ts,4,MPI_DOUBLE,pid,698,MPI_COMM_WORLD);

	}

	timesteplim_plasma=ts[0];
	timesteplim_convection=ts[1];
	timesteplim_diffusion=ts[2];
	timesteplim_pospresv=ts[3];

	//----------------------------------------------------------------------------------------------------------------------------------------

}


double globalvar::pospresvtimestep()
{

	double iontimestepreq,electrontimestepreq,maxden=-2,ret;
	int i;
	double localden1,localden2;

	for(i=1;i<Nlocal;i++)
	{

		if(nirk[2*i-1]>0.0)
			localden1=D_i*(beta_i[i-1]+qi[2*i-1]/nirk[2*i-1]);
		else
			localden1=D_i*beta_i[i-1];

		if(nirk[2*i]>0.0)
			localden2=D_i*(beta_i[i]-qi[2*i]/nirk[2*i]);
		else
			localden2=D_i*beta_i[i];

		if(localden1<0)
		{
			cout<<"ilocalden1<0 during timetseps calculations at "<<N1+i<<endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		if(localden2<0)
		{
			cout<<"ilocalden2<0 during timetseps calculations at "<<N1+i<<endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		if(maxden<max(localden1,localden2))
			maxden=max(localden1,localden2);

	}

	iontimestepreq=(1.0/3.0)*dx/maxden;

	if(dgaci>4)
		iontimestepreq/=2.0;

	if(dgaci>6)
		iontimestepreq/=2.0;


	//--------------------------------------------------------------------------------------------------------------------------------------------------


	for(i=1;i<Nlocal;i++)
	{

		if(nerk[2*i-1]>0.0)
			localden1=D_e*(beta_e[i-1]+qe[2*i-1]/nerk[2*i-1]);
		else
			localden1=D_e*beta_e[i-1];

		if(nerk[2*i]>0.0)
			localden2=D_e*(beta_e[i]-qe[2*i]/nerk[2*i]);
		else
			localden2=D_e*beta_e[i];

		if(localden1<0)
		{
			cout<<"elocalden1<0 during timetseps calculations at "<<N1+i<<endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		if(localden2<0)
		{
			cout<<"elocalden2<0 during timetseps calculations at "<<N1+i<<endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		if(maxden<max(localden1,localden2))
			maxden=max(localden1,localden2);

	}


	electrontimestepreq=(1.0/3.0)*dx/maxden;

	if(dgace>4)
		electrontimestepreq/=2.0;

	if(dgace>6)
		electrontimestepreq/=2.0;

	//--------------------------------------------------------------------------------------------------------------------------------------------------

	ret=min(iontimestepreq,electrontimestepreq);


	return ret;

}



