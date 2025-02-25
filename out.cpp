
#include "out.h"

using namespace std;

char Filename[250];

void globalvar::output()
{
	double Dummy1[7],Dummy2[7], Dummy3[7],Dummy4[7];
	MPI_Status stat,stat2;
	MPI_Request req;

	if(procid>0)
	{
		Dummy1[0]=ni[1];
		Dummy1[1]=ne[1];
		Dummy1[2]=phi[1];
		Dummy1[3]=Ef[1];
		Dummy1[4]=alpha[1];
		Dummy1[5]=qi[1];
		Dummy1[6]=qe[1];

		MPI_Isend(&Dummy1[0],7,MPI_DOUBLE,procid-1,699875,MPI_COMM_WORLD,&req);

		//-------------------------------------------------------------------------------------------------------------------------------------------

		MPI_Recv(&Dummy2[0],7,MPI_DOUBLE,procid-1,688975,MPI_COMM_WORLD,&stat2);
		ni[0]=Dummy2[0];
		ne[0]=Dummy2[1];
		phi[0]=Dummy2[2];
		Ef[0]=Dummy2[3];
		alpha[0]=Dummy3[4];
		qi[0]=Dummy2[5];
		qe[0]=Dummy2[6];

		MPI_Wait(&req,&stat);
	}

	if(procid!=procsize-1)
	{
		Dummy3[0]=ni[2*Nlocal-2];
		Dummy3[1]=ne[2*Nlocal-2];
		Dummy3[2]=phi[2*Nlocal-2];
		Dummy3[3]=Ef[2*Nlocal-2];
		Dummy3[4]=alpha[Nlocal-1];
		Dummy3[5]=qi[2*Nlocal-2];
		Dummy3[6]=qe[2*Nlocal-2];

		MPI_Isend(&Dummy3[0],7,MPI_DOUBLE,procid+1,688975,MPI_COMM_WORLD,&req);

		//-------------------------------------------------------------------------------------------------------------------------------------------

		MPI_Recv(&Dummy4[0],7,MPI_DOUBLE,procid+1,699875,MPI_COMM_WORLD,&stat2);
		ni[2*Nlocal-1]=Dummy4[0];
		ne[2*Nlocal-1]=Dummy4[1];
		phi[2*Nlocal-1]=Dummy4[2];
		Ef[2*Nlocal-1]=Dummy4[3];
		alpha[Nlocal]=Dummy4[4];
		qi[2*Nlocal-1]=Dummy4[5];
		qe[2*Nlocal-1]=Dummy4[6];

		MPI_Wait(&req,&stat);

	}

	if(procid==0)
	{
		ni[0]=nirk[0];
		ne[0]=nerk[0];
		Ef[0]=Ef[1];
		alpha[0]=alpha[1];
	}

	if(procid==procsize-1)
	{
		ni[2*Nlocal-1]=nirk[2*Nlocal-1];
		ne[2*Nlocal-1]=nerk[2*Nlocal-1];
		Ef[2*Nlocal-1]=Ef[2*Nlocal-2];
		alpha[Nlocal]=alpha[Nlocal-1];
	}
	//-------------------------------------------------------------------------------------------------------------------------------------------------------------

	if(procid!=0)
	{
		int j;
		double *UK;
		UK=new double[Nlocal*9];

		int count=0;
		for(j=0;j<Nlocal;j++)
		{
			UK[count++]=(ne[2*j]+ne[2*j+1])/2.0/pow(10,15);
			UK[count++]=(ni[2*j]+ni[2*j+1])/2.0/pow(10,15);
			UK[count++]=(phi[2*j]+phi[2*j+1])/2.0;
			UK[count++]=(alpha[j]+alpha[j+1])/2.0;
			UK[count++]=-(Ef[2*j]+Ef[2*j+1])/2.0;
			UK[count++]=-e*(phi[2*j]+phi[2*j+1])/2.0/(boltz_const*T_e*11600);
			UK[count++]=-(Ef[2*j]+Ef[2*j+1])/2.0/p;
			UK[count++]=mu_e*e*(ne[2*j]+ne[2*j+1])/2.0+D_e*(qe[2*j]+qe[2*j+1])/2.0;
			UK[count++]=mu_i*e*(ni[2*j]+ni[2*j+1])/2.0-D_i*(qi[2*j]+qi[2*j+1])/2.0;
		}


		MPI_Send(&Nlocal,1,MPI_INT,0,200+procid,MPI_COMM_WORLD);
		MPI_Send(UK,Nlocal*9,MPI_DOUBLE,0,20000+procid,MPI_COMM_WORLD);


		delete [] UK;
	}

	//-------------------------------------------------------------------------------------------------------------------------------------------------------------

	if(procid==0)
	{
		ofstream File,file;
		int i,j,count;

		char q='"';
		sprintf(Filename,"a_%d.dat",outputcount);
		File.open (Filename);

		File<<"Title="<<q<<"Solutions at time = "<<t<<q<<endl;
		File<<"Variables="<<q
				<<"X"<<q<<","<<q
				<<"n<sub>e</sub>"<<q<<","<<q
				<<"n<sub>i</sub>"<<q<<","<<q
				<<"<greek>f</greek>"<<q<<","<<q
				<<"<greek>a</greek>"<<q<<","<<q
				<<"E<sub>f</sub>"<<q<<","<<q
				<<"e<greek>f</greek>/(k<sub>B</sub>T<sub>e</sub>"<<q<<","<<q
				<<"Eop"<<q
				<<"j<greek>ion</greek>"<<q<<","<<q
				<<"j<greek>electron</greek>"<<q<<","<<q<<endl;


		for(j=0;j<Nlocal;j++)
		{
			File<<(j)*dx<<" "
					<<(ne[2*j]+ne[2*j+1])/2/pow(10,15)<<" "
					<<(ni[2*j]+ni[2*j+1])/2/pow(10,15)<<" "
					<<(phi[2*j]+phi[2*j+1])/2<<" "
					<<(alpha[j]+alpha[j+1])/2.0<<" "
					<<-(Ef[2*j]+Ef[2*j+1])/2<<" "
					<<-e*(phi[2*j]+phi[2*j+1])/2/(boltz_const*T_e*11600)<<" "
					<<-(Ef[2*j]+Ef[2*j+1])/2/p<<" "
					<<mu_e*e*(ne[2*j]+ne[2*j+1])/2.0+D_e*(qe[2*j]+qe[2*j+1])/2.0<<" "
					<<mu_i*e*(ni[2*j]+ni[2*j+1])/2.0-D_i*(qi[2*j]+qi[2*j+1])/2.0<<endl;

		}

		//-------------------------------------------------------------------------------------------------------------------------------------------------------------


		MPI_Status sta;
		count=Nlocal;
		for(i=1;i<procsize;i++)
		{
			int nalloc;
			MPI_Recv(&nalloc,1,MPI_INT,i,200+i,MPI_COMM_WORLD,&sta);
			double *UUK;
			UUK=new double [nalloc*9];
			MPI_Recv(UUK,nalloc*9,MPI_DOUBLE,i,20000+i,MPI_COMM_WORLD,&sta);

			int count2=0;
			for(j=0;j<nalloc;j++)
			{
				File<<(count+j-1)*dx<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<endl;
			}

			count=count+nalloc-1;
			delete [] UUK;
		}
		if(count!=N)
			cout<<"Output calcs are wrong"<<endl;

		file.close();
		File.close();


		outputcount++;

	}
}



