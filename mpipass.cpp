#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <mpi.h>

#include "global_var.h"

using namespace std;



void globalvar::mpipassu()
{
	double dummy1[2],dummy2[2],dummy3[2],dummy4[2];

	MPI_Request stat;
	MPI_Status stat2;

	if(procid>0)
	{

		dummy1[0]=nirk[1];
		dummy1[1]=nerk[1];
		MPI_Send(&dummy1[0],2,MPI_DOUBLE,procid-1,100+iteration,MPI_COMM_WORLD);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&dummy2[0],2,MPI_DOUBLE,procid-1,200+iteration,MPI_COMM_WORLD,&stat2);
		nirk[0]=dummy2[0];
		nerk[0]=dummy2[1];
	}

	if(procid!=procsize-1)
	{
		dummy3[0]=nirk[2*Nlocal-2];
		dummy3[1]=nerk[2*Nlocal-2];
		MPI_Send(&dummy3[0],2,MPI_DOUBLE,procid+1,200+iteration,MPI_COMM_WORLD);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&dummy4[0],2,MPI_DOUBLE,procid+1,100+iteration,MPI_COMM_WORLD,&stat2);
		nirk[2*Nlocal-1]=dummy4[0];
		nerk[2*Nlocal-1]=dummy4[1];
	}
}

void globalvar::mpipassq()
{

	MPI_Request stat;
	MPI_Status stat2;
	double Dummy1[3], Dummy2[3], Dummy3[3], Dummy4[3];

	if(procid>0)
	{
		Dummy1[0]=qi[1];
		Dummy1[1]=qe[1];
		Dummy1[2]=Ef[1];

		MPI_Send(&Dummy1[0],3,MPI_DOUBLE,procid-1,300+iteration,MPI_COMM_WORLD);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&Dummy2[0],3,MPI_DOUBLE,procid-1,400+iteration,MPI_COMM_WORLD,&stat2);
		qi[0]=Dummy2[0];
		qe[0]=Dummy2[1];
		Ef[0]=Dummy2[2];
	}

	if(procid!=procsize-1)
	{
		Dummy3[0]=qi[2*Nlocal-2];
		Dummy3[1]=qe[2*Nlocal-2];
		Dummy3[2]=Ef[2*Nlocal-2];
		MPI_Send(&Dummy3[0],3,MPI_DOUBLE,procid+1,400+iteration,MPI_COMM_WORLD);
		//-------------------------------------------------------------------------------------------------------------------------------------------
		MPI_Recv(&Dummy4[0],3,MPI_DOUBLE,procid+1,300+iteration,MPI_COMM_WORLD,&stat2);
		qi[2*Nlocal-1]=Dummy4[0];
		qe[2*Nlocal-1]=Dummy4[1];
		Ef[2*Nlocal-1]=Dummy4[2];
	}



}
