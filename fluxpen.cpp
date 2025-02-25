
#include "fluxpen.h"

using namespace std;

void globalvar::fluxpencalc()
{
	int i;
	double a1,b1;


	if(nirk[1]<eps)
		if(qi[1]!=0.0)
			if(t>0.0)
			{
				cout<<"Need to change q at time "<<t<<endl;
				qi[1]=0.0;
				MPI_Abort(MPI_COMM_WORLD,1);
			}

	if(nerk[1]<eps)
		if(qe[1]!=0.0)
			if(t>0.0)
			{
				cout<<"Need to change q at time "<<t<<endl;
				qe[1]=0.0;
				MPI_Abort(MPI_COMM_WORLD,1);
			}


	for(i=1;i<Nlocal-1;i++)
	{
		if(nirk[2*i]<eps)
			if(qi[2*i]!=0.0)
				if(t>0.0)
				{
					cout<<"Need to change q of cell "<<i<<" at time "<<t<<endl;
					qi[2*i]=0.0;
					MPI_Abort(MPI_COMM_WORLD,1);
				}

		if(nirk[2*i+1]<eps)
			if(qi[2*i+1]!=0.0)
				if(t>0.0)
				{
					cout<<"Need to change q of cell "<<i<<" at time "<<t<<endl;
					qi[2*i+1]=0.0;
					MPI_Abort(MPI_COMM_WORLD,1);
				}

		if(nerk[2*i]<eps)
			if(qe[2*i]!=0.0)
				if(t>0.0)
				{
					cout<<"Need to change q of cell "<<i<<" at time "<<t<<endl;
					qe[2*i]=0.0;
					MPI_Abort(MPI_COMM_WORLD,1);
				}

		if(nerk[2*i+1]<eps)
			if(qe[2*i+1]!=0.0)
				if(t>0.0)
				{
					cout<<"Need to change q of cell "<<i<<" at time "<<t<<endl;
					qe[2*i+1]=0.0;
					MPI_Abort(MPI_COMM_WORLD,1);
				}
	}



	if(nirk[2*Nlocal-2]<eps)
		if(qi[2*Nlocal-2]!=0.0)
			if(t>0.0)
			{
				cout<<"Need to change q at time "<<t<<endl;
				qi[2*Nlocal-2]=0.0;
				MPI_Abort(MPI_COMM_WORLD,1);
			}

	if(nerk[2*Nlocal-2]<eps)
		if(qe[2*Nlocal-2]!=0.0)
			if(t>0.0)
			{
				cout<<"Need to change q at time "<<t<<endl;
				qe[2*Nlocal-2]=0.0;
				MPI_Abort(MPI_COMM_WORLD,1);
			}





	//---------------------------------------------------------------------------------------------------------------------------------------------------


	for(i=0;i<Nlocal;i++)
	{
		a1=0.0;
		b1=0.0;


		a1=abs(qi[2*i]/nirk[2*i]);
		b1=abs(qi[2*i+1]/nirk[2*i+1]);

		if(nirk[2*i]==0)
			a1=0.0;

		if(nirk[2*i+1]==0)
			b1=0.0;

		if(i==0 && procid==0)
			a1=0.0;
		if(i==Nlocal-1 && procid==procsize-1)
			b1=0.0;

		beta_i[i]=1.01*max(a1,b1);  //because beta > q/u .


	}

	for(i=0;i<Nlocal;i++)
	{
		a1=0.0;
		b1=0.0;


		a1=abs(qe[2*i]/nerk[2*i]);
		b1=abs(qe[2*i+1]/nerk[2*i+1]);

		if(nerk[2*i]==0)
			a1=0.0;

		if(nerk[2*i+1]==0)
			b1=0.0;

		if(i==0 && procid==0)
			a1=0.0;
		if(i==Nlocal-1 && procid==procsize-1)
			b1=0.0;

		beta_e[i]=1.01*max(a1,b1);  //because beta > q/u .



	}

}




