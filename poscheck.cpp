
#include "poscheck.h"

using namespace std;

void globalvar::positivity_check()
{
	MPI_Status stat2;
	int i,flag=0;

	for(i=1;i<Nlocal;i++)
		if(new_ni_k[i][0]<0.0)
			flag=1;

	for(i=1;i<Nlocal;i++)
		if(new_ne_k[i][0]<0.0)
			flag=1;


	if(procid>0)
	{
		MPI_Send(&flag,1,MPI_INT,0,550,MPI_COMM_WORLD);
		MPI_Recv(&flag,1,MPI_INT,0,598,MPI_COMM_WORLD,&stat2);
	}
	if(procid==0)
	{
		int msg;
		for(int pid=1;pid<procsize;pid++)
		{
			MPI_Recv(&msg,1,MPI_INT,pid,550,MPI_COMM_WORLD,&stat2);
			if(msg)
				flag=msg;
		}
		for(int pid=1;pid<procsize;pid++)
			MPI_Send(&flag,1,MPI_INT,pid,598,MPI_COMM_WORLD);

	}

	if(flag==1)
	{
		dt=dt/2.0;

		//calctimesteps(); 	//making sure the right time step limits are used.

		if(dt<timesteplim_pospresv*0.5 && dt<timesteplim_diffusion*0.5 && dt<timesteplim_convection*0.5)
		{
			if(procid==0)
			{
				cout<<"DANGER: Timestep is smaller than lowest limit required at time t = "<<t<<endl;
				cout<<"Timestep in use: "<<dt*2.0<<endl;
				cout<<"Timestep limits (convection, diffusion, pos presv): "<<timesteplim_convection<<" "<<timesteplim_diffusion<<" "<<timesteplim_pospresv<<endl;
				cout<<"Exiting"<<endl<<endl;
				update();
				outputcount=12056;
				output();
				MPI_Abort(MPI_COMM_WORLD,2);
			}
		}

		solve();
	}

}

