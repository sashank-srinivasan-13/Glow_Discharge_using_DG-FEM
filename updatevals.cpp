
#include "updatevals.h"

using namespace std;

void globalvar::update()
{

	int i,j;
	for(i=1;i<Nlocal;i++)
	{
		for(j=0;j<dgaci;j++)
			if(isnan(new_ni_k[i][j]))
			{
				cout<<"n_i Nan detected at cell: "<<N1+i<<", level: "<<j<<endl;
				MPI_Abort(MPI_COMM_WORLD,1);
			}

		for(j=0;j<dgace;j++)
			if(isnan(new_ne_k[i][j]))
			{
				cout<<"n_e Nan detected at cell: "<<N1+i<<", level: "<<j<<endl;
				MPI_Abort(MPI_COMM_WORLD,1);
			}
	}

	for(i=1;i<Nlocal;i++)
	{
		for(j=0;j<dgaci;j++)
			ni_k[i][j]=new_ni_k[i][j];

		for(j=0;j<dgace;j++)
			ne_k[i][j]=new_ne_k[i][j];

		ni[2*i-1]=0;
		qi[2*i-1]=0;
		ni[2*i]=0;
		qi[2*i]=0;

		for(j=0;j<dgaci;j++)
		{
			ni[2*i-1]+=ni_k[i][j]*pow(-1,j);
			qi[2*i-1]+=qi_k[i][j]*pow(-1,j);

			ni[2*i]+=ni_k[i][j];
			qi[2*i]+=qi_k[i][j];
		}

		ne[2*i-1]=0;
		qe[2*i-1]=0;
		ne[2*i]=0;
		qe[2*i]=0;

		for(j=0;j<dgace;j++)
		{
			ne[2*i-1]+=ne_k[i][j]*pow(-1,j);
			qe[2*i-1]+=qe_k[i][j]*pow(-1,j);

			ne[2*i]+=ne_k[i][j];
			qe[2*i]+=qe_k[i][j];
		}

		phi[2*i-1]=0;
		qphi[2*i-1]=0;
		phi[2*i]=0;
		qphi[2*i]=0;

		for(j=0;j<dgacp;j++)
		{
			phi[2*i-1]+=phi_k[i][j]*pow(-1,j);
			qphi[2*i-1]+=qphi_k[i][j]*pow(-1,j);

			phi[2*i]+=phi_k[i][j];
			qphi[2*i]+=qphi_k[i][j];
		}
	}



}


