#include "backup.h"

using namespace std;




void globalvar::back_up()
{
	double Dummy1[5],Dummy2[5],Dummy3[5],Dummy4[5];
	MPI_Status stat2;
	if(procid>0)
	{
		Dummy1[0]=ni[1];
		Dummy1[1]=ne[1];
		Dummy1[2]=phi[1];
		Dummy1[3]=Ef[1];
		Dummy1[4]=alpha[1];

		MPI_Send(&Dummy1[0],5,MPI_DOUBLE,procid-1,699875,MPI_COMM_WORLD);

		//-------------------------------------------------------------------------------------------------------------------------------------------

		MPI_Recv(&Dummy2[0],5,MPI_DOUBLE,procid-1,688975,MPI_COMM_WORLD,&stat2);
		ni[0]=Dummy2[0];
		ne[0]=Dummy2[1];
		phi[0]=Dummy2[2];
		Ef[0]=Dummy2[3];
		alpha[0]=Dummy2[4];
	}

	if(procid!=procsize-1)
	{
		Dummy3[0]=ni[2*Nlocal-2];
		Dummy3[1]=ne[2*Nlocal-2];
		Dummy3[2]=phi[2*Nlocal-2];
		Dummy3[3]=Ef[2*Nlocal-2];
		Dummy3[4]=alpha[2*Nlocal-2];

		MPI_Send(&Dummy3[0],5,MPI_DOUBLE,procid+1,688975,MPI_COMM_WORLD);

		//-------------------------------------------------------------------------------------------------------------------------------------------

		MPI_Recv(&Dummy4[0],5,MPI_DOUBLE,procid+1,699875,MPI_COMM_WORLD,&stat2);
		ni[2*Nlocal-1]=Dummy4[0];
		ne[2*Nlocal-1]=Dummy4[1];
		phi[2*Nlocal-1]=Dummy4[2];
		Ef[2*Nlocal-1]=Dummy4[3];
		alpha[2*Nlocal-1]=Dummy4[4];
	}

	if(procid==0)
	{
		ni[0]=ni[1];
		ne[0]=ne[1];
		phi[0]=phi[1];
		Ef[0]=Ef[1];
		alpha[0]=alpha[1];
	}

	if(procid==procsize-1)
	{
		ni[2*Nlocal-1]=ni[2*Nlocal-2];
		ne[2*Nlocal-1]=ne[2*Nlocal-2];
		phi[2*Nlocal-1]=phi[2*Nlocal-2];
		Ef[2*Nlocal-1]=Ef[2*Nlocal-2];
		alpha[2*Nlocal-1]=alpha[2*Nlocal-2];
	}

	//-------------------------------------------------------------------------------------------------------------------------------------------------------------

	if(procid!=0)
	{
		int i,j;

		double *uk;
		double *UK;
		uk=new double[(dgaci+dgace+dgacp)*(Nlocal-1)*2];
		UK=new double[Nlocal*8];

		int count=0;
		for(int i=1;i<Nlocal;i++)
		{
			for(int j=0;j<dgaci;j++)
			{
				uk[count++]=ni_k[i][j];
				uk[count++]=qi_k[i][j];
			}
			for(int j=0;j<dgace;j++)
			{
				uk[count++]=ne_k[i][j];
				uk[count++]=qe_k[i][j];
			}
			for(int j=0;j<dgacp;j++)
			{
				uk[count++]=phi_k[i][j];
				uk[count++]=qphi_k[i][j];
			}
		}

		//-------------------------------------------------------------------------------------------------------------------------------------------------------------


		count=0;
		for(j=0;j<Nlocal;j++)
		{
			UK[count++]=(ni[2*j]+ni[2*j+1])/2/pow(10,15);
			UK[count++]=(ne[2*j]+ne[2*j+1])/2/pow(10,15);
			UK[count++]=(phi[2*j]+phi[2*j+1])/2;
			UK[count++]=(alpha[j]+alpha[j+1])/2.0;
			UK[count++]=-(Ef[2*j]+Ef[2*j+1])/2;
			UK[count++]=-e*(phi[2*j]+phi[2*j+1])/2/(boltz_const*T_e*11600);
			UK[count++]=-(Ef[2*j]+Ef[2*j+1])/2/p;

			if(j>0)
				UK[count++]=e*(phi[2*j]-phi[2*j-1])/(boltz_const*T_e*11600);
			else
				UK[count++]=0;
		}


		MPI_Send(&Nlocal,1,MPI_INT,0,200+procid,MPI_COMM_WORLD);
		MPI_Send(uk,(dgaci+dgace+dgacp)*(Nlocal-1)*2,MPI_DOUBLE,0,2000+procid,MPI_COMM_WORLD);
		MPI_Send(UK,Nlocal*8,MPI_DOUBLE,0,20000+procid,MPI_COMM_WORLD);


		delete [] uk;
		delete [] UK;
	}

	//-------------------------------------------------------------------------------------------------------------------------------------------------------------

	if(procid==0)
	{
		ofstream file;
		ofstream File;
		ofstream nik2,nek2,phik2;

		file.open("a_bakdata.dat");
		file<<N<<endl;
		file<<dgaci<<endl;
		file<<dgace<<endl;
		file<<dgacp<<endl;
		file<<i_limiter_max<<endl;
		file<<e_limiter_max<<endl;
		file<<t<<endl;
		file<<outputcount<<endl;
		file<<delv;
		file.close();

		//-------------------------------------------------------------------------------------------------------------------------------------------------------------


		int i,j,k,count;

		file.open("a_bak.dat");
		for(i=1;i<Nlocal;i++)
		{
			file<<i<<endl;
			for(j=0;j<dgaci;j++)
			{
				file<<ni_k[i][j]<<endl;
				file<<qi_k[i][j]<<endl;
			}

			for(j=0;j<dgace;j++)
			{
				file<<ne_k[i][j]<<endl;
				file<<qe_k[i][j]<<endl;
			}

			for(j=0;j<dgacp;j++)
			{
				file<<phi_k[i][j]<<endl;
				file<<qphi_k[i][j]<<endl;
			}

		}

		//-------------------------------------------------------------------------------------------------------------------------------------------------------------


		char q='"';
		File.open ("a_backup.plt");
		File<<"Variables="<<q
				<<"X"<<q<<","<<q
				<<"n<sub>i</sub>"<<q<<","<<q
				<<"n<sub>e</sub>"<<q<<","<<q
				<<"<greek>f</greek>"<<q<<","<<q
				<<"<greek>a</greek>"<<q<<","<<q
				<<"E<sub>f</sub>"<<q<<","<<q
				<<"e<greek>f</greek>/(k<sub>B</sub>T<sub>e</sub>"<<q<<","<<q
				<<"Eop"<<q<<endl;

		for(i=0;i<Nlocal;i++)
		{
			File<<(i)*dx<<" "
					<<(ni[2*i]+ni[2*i+1])/2/pow(10,15)<<" "
					<<(ne[2*i]+ne[2*i+1])/2/pow(10,15)<<" "
					<<(phi[2*i]+phi[2*i+1])/2<<" "
					<<(alpha[j]+alpha[j+1])/2.0<<" "
					<<-(Ef[2*i]+Ef[2*i+1])/2<<" "
					<<-e*(phi[2*i]+phi[2*i+1])/2/(boltz_const*T_e*11600)<<" "
					<<-(Ef[2*i]+Ef[2*i+1])/2/p<<endl;

		}

		*potentialfile<<"ZONE"<<endl;
		for(i=1;i<Nlocal;i++)
			*potentialfile<<(i-0.5)*dx<<" "
			<<phi_k[i][0]<<endl;

		*potdropfile<<"ZONE"<<endl;
		for(i=1;i<Nlocal;i++)
			*potdropfile<<(i-0.5)*dx<<" "
			<<e*(phi[2*i]-phi[2*i-1])/(boltz_const*T_e*11600)<<" "<<endl;


		//-------------------------------------------------------------------------------------------------------------------------------------------------------------


		MPI_Status sta;
		count=Nlocal;
		for(i=1;i<procsize;i++)
		{
			int nalloc;
			MPI_Recv(&nalloc,1,MPI_INT,i,200+i,MPI_COMM_WORLD,&sta);

			double *uuk;
			double *UUK;
			uuk=new double [(dgaci+dgace+dgacp)*(nalloc-1)*2];
			UUK=new double [nalloc*8];


			MPI_Recv(uuk,(dgaci+dgace+dgacp)*(nalloc-1)*2,MPI_DOUBLE,i,2000+i,MPI_COMM_WORLD,&sta);
			MPI_Recv(UUK,nalloc*8,MPI_DOUBLE,i,20000+i,MPI_COMM_WORLD,&sta);

			int count2=0;
			for(j=1;j<nalloc;j++)
			{
				file<<j+count-1<<endl;
				for(k=0;k<dgaci;k++)
				{
					file<<uuk[count2++]<<endl;		//ni_k
					file<<uuk[count2++]<<endl;		//qi_k
				}

				for(k=0;k<dgace;k++)
				{
					file<<uuk[count2++]<<endl;		//ne_k
					file<<uuk[count2++]<<endl;		//qe_k
				}

				for(k=0;k<dgacp;k++)
				{
					if(k==0)
					{
						*potentialfile<<(count+j-0.5)*dx<<" "
								<<uuk[count2]<<endl;
					}

					file<<uuk[count2++]<<endl;		//phi_k
					file<<uuk[count2++]<<endl;		//qphi_k
				}
			}

			//-------------------------------------------------------------------------------------------------------------------------------------------------------------


			count2=0;
			for(j=0;j<nalloc;j++)
			{
				File<<(count+j-1)*dx<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<" ";
				File<<UUK[count2++]<<endl;

				if(j>0)
					*potdropfile<<(count+j-0.5)*dx<<" "<<UUK[count2]<<" "<<endl;
				count2++;
			}

			count=count+nalloc-1;
			delete [] uuk;
			delete [] UUK;
		}
		if(count!=N)
			cout<<"Output calcs are wrong"<<endl;

		file.close();
		File.close();
	}

}

