#include "limiter.h"
using namespace std;

double minmodbeg(double a, double b)
{
	double ret=0;
	if(a>=0 && b>=0 )
		ret=min(a,b);
	if(a<=0 && b<=0 )
		ret=min(a,b);
	else
		ret=0;

	return ret;
}


double minmod(double a, double b, double c)
{
	double ret=0;
	if(a>=0 && b>=0 && c>=0)
	{
		ret=min(a,b);
		ret=min(ret,c);
	}

	else if(a<0 && b<0 && c<0)
	{
		ret=min(a,b);
		ret=min(ret,c);
	}

	else
		ret=0;


	return ret;
}


void globalvar::slope_limiter()
{
	int i,j,flag;
	double Dipk,Dimk,ans,limcoeff;
	for(i=1;i<Nlocal;i++)
	{
		if(abs(phi[2*i-1]-phi[2*i])>1.0*boltz_const*T_e*11600.0/e)
			limflag[i]=1;
		else
		{
			if(limchoice==0)
				limflag[i]=0;
			else
				limflag[i]=1;
		}
	}

	//------------------------------------------------------------------------------------------------------------------------------------------------

	double dummy1[dgaci+dgace], dummy2[dgaci+dgace], dummy3[dgaci+dgace], dummy4[dgaci+dgace];
	MPI_Status stat2;

	if(procid>0)
	{
		for(j=0;j<dgaci;j++)
			dummy1[j]=new_ni_k[1][j];
		for(j=0;j<dgace;j++)
			dummy1[j+dgaci]=new_ne_k[1][j];
		MPI_Send(&dummy1[0],dgaci+dgace,MPI_DOUBLE,procid-1,99875,MPI_COMM_WORLD);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&dummy2[0],dgaci+dgace,MPI_DOUBLE,procid-1,88975,MPI_COMM_WORLD,&stat2);

		for(j=0;j<dgaci;j++)
			new_ni_k[0][j]=dummy2[j];
		for(j=0;j<dgace;j++)
			new_ne_k[0][j]=dummy2[j+dgaci];

	}

	if(procid!=procsize-1)
	{
		for(j=0;j<dgaci;j++)
			dummy3[j]=new_ni_k[Nlocal-1][j];
		for(j=0;j<dgace;j++)
			dummy3[j+dgaci]=new_ne_k[Nlocal-1][j];
		MPI_Send(&dummy3[0],dgaci+dgace,MPI_DOUBLE,procid+1,88975,MPI_COMM_WORLD);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&dummy4[0],dgaci+dgace,MPI_DOUBLE,procid+1,99875,MPI_COMM_WORLD,&stat2);
		for(j=0;j<dgaci;j++)
			new_ni_k[Nlocal][j]=dummy4[j];
		for(j=0;j<dgace;j++)
			new_ne_k[Nlocal][j]=dummy4[j+dgaci];
	}


	//------------------------------------------------------------------------------------------------------------------------------------------------


	for(i=1;i<Nlocal;i++)
	{
		if(i==1 && procid==0)
			for(j=dgaci-1;j>0;j--)
			{
				if(limflag[i]==0)
					limcoeff=1.0;
				else
					limcoeff=0.5/(2*j-1);


				Dipk=(new_ni_k[i+1][j-1]-new_ni_k[i][j-1]);
				Dimk=0;
				new_ni_k[i][j]=minmodbeg(new_ni_k[i][j],Dipk*limcoeff);

			}
		else if(i==Nlocal-1 && procid==procsize-1)
			for(j=dgaci-1;j>0;j--)
			{
				if(limflag[i]==0)
					limcoeff=1.0;
				else
					limcoeff=0.5/(2*j-1);

				Dipk=0;
				Dimk=(new_ni_k[i][j-1]-new_ni_k[i-1][j-1]);
				new_ni_k[i][j]=minmodbeg(new_ni_k[i][j],Dimk*limcoeff);
			}
		else
		{
			flag=0;
			if(flag==0 && dgaci>i_limiter_max)
				for(j=dgaci-1;j>i_limiter_max-1;j--)
				{
					if(limflag[i]==0)
						limcoeff=1.0;
					else
						limcoeff=0.5/(2*j-1);

					Dipk=(new_ni_k[i+1][j-1]-new_ni_k[i][j-1]);
					Dimk=(new_ni_k[i][j-1]-new_ni_k[i-1][j-1]);
					ans=minmod(new_ni_k[i][j],Dipk*limcoeff,Dimk*limcoeff);
					if(ans==new_ni_k[i][j])
						flag=1;
					new_ni_k[i][j]=ans;
				}
		}
	}

	//--------------------------------------------------------------------------------------------------------------------------------------

	for(i=1;i<Nlocal;i++)
	{
		if(i==1 && procid==0)
			for(j=dgace-1;j>0;j--)
			{
				if(limflag[i]==0)
					limcoeff=1.0;
				else
					limcoeff=0.5/(2*j-1);
				Dipk=(new_ne_k[i+1][j-1]-new_ne_k[i][j-1]);
				Dimk=0;
				new_ne_k[i][j]=minmodbeg(new_ne_k[i][j],Dipk*limcoeff);
			}
		else if(i==Nlocal-1 && procid==procsize-1)
			for(j=dgace-1;j>0;j--)
			{
				if(limflag[i]==0)
					limcoeff=1.0;
				else
					limcoeff=0.5/(2*j-1);

				Dipk=0;
				Dimk=(new_ne_k[i][j-1]-new_ne_k[i-1][j-1]);
				new_ne_k[i][j]=minmodbeg(new_ne_k[i][j],Dimk*limcoeff);
			}
		else
		{
			flag=0;
			if(flag==0 && dgace>e_limiter_max)
				for(j=dgace-1;j>e_limiter_max-1;j--)
				{
					if(limflag[i]==0)
						limcoeff=1.0;
					else
						limcoeff=0.5/(2*j-1);

					Dipk=(new_ne_k[i+1][j-1]-new_ne_k[i][j-1]);
					Dimk=(new_ne_k[i][j-1]-new_ne_k[i-1][j-1]);
					ans=minmod(new_ne_k[i][j],Dipk*limcoeff,Dimk*limcoeff);
					if(ans==new_ne_k[i][j])
						flag=1;
					new_ne_k[i][j]=ans;

				}
		}
	}

	//--------------------------------------------------------------------------------------------------------------------------------------

}

void globalvar::positivity_limiter()
{
	double Points[12]={-1.0,+1.0,-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,-0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,-0.9739065285171717,0.9739065285171717};

	int i,j,falg=0;
	double theta,xi,m;

	for(i=1;i<Nlocal;i++)
	{
		m=0.0;
		for(j=0;j<12;j++)
		{
			xi=Points[j];
			if(int_newni(i,xi)<m)
				m=int_newni(i,xi);
		}

		theta=min(1.0,abs((0.0-new_ni_k[i][0])/(m-new_ni_k[i][0])));

		if(m<0)
			theta=theta-eps;

		for(j=1;j<dgaci;j++)
			new_ni_k[i][j]=theta*new_ni_k[i][j];

	}

	for(i=1;i<Nlocal;i++)
	{
		m=0.0;
		for(j=0;j<12;j++)
		{
			xi=Points[j];
			if(int_newne(i,xi)<m)
				m=int_newne(i,xi);
		}

		theta=min(1.0,abs((0.0-new_ne_k[i][0])/(m-new_ne_k[i][0])));

		if(m<0)
			theta=theta-eps;

		for(j=1;j<dgaci;j++)
			new_ne_k[i][j]=theta*new_ne_k[i][j];

	}


}

