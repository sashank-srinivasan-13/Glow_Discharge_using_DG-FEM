
#include "qreshape.h"

using namespace std;

double Quad_Points[7]={-1.0, 1.0 , 0.0 , -0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472};

void globalvar::reshapeq()
{
	int i,flag;

	for(i=1;i<Nlocal;i++)
	{
		flag=0;
		if(nirk[2*i-1]<eps && qi[2*i-1]!=0.0)
		{
			qi[2*i-1]=0.0;flag=1;
		}
		if(nirk[2*i]<eps && qi[2*i]!=0.0)
		{
			qi[2*i]=0.0;flag=1;
		}
		if(flag!=0)
			findniqk(i);

		//------------------------------------------------------------

		flag=0;
		if(nerk[2*i-1]<eps && qe[2*i-1]!=0.0)
		{
			qe[2*i-1]=0.0;flag=1;
		}
		if(nerk[2*i]<eps && qe[2*i]!=0.0)
		{
			qe[2*i]=0.0;flag=1;
		}
		if(flag!=0)
			findneqk(i);
	}

}

void globalvar::findniqk(int i)
{
	int dim,j;
	int *IPIV=new int[dgaci];
	int INFO;
	int NRHS;
	double nir[dgaci];
	double re[dgaci][dgaci],xi;

	for(int m=0;m<dgaci;m++)
	{
		nir[m]=reconstructqi(i,m);
		xi=Quad_Points[m];
		for(j=0;j<dgaci;j++)
			re[j][m]=int_P(xi,j);
	}

	nir[0]=qi[2*i-1];
	nir[1]=qi[2*i];

	dim=dgaci;
	NRHS=1;
	dgesv_(&dim,&NRHS,&re[0][0],&dim,&IPIV[0],&nir[0],&dim,&INFO);
	if(INFO!=0)
		cout<<"QiReshape Linear system solve not correct: INFO = "<<INFO<<endl;

	for(j=0;j<dgaci;j++)
		qi_k[i][j]=nir[j];

	delete [] IPIV;

}

double globalvar::reconstructqi(int i,int j)
{
	double ret,xi;
	xi=Quad_Points[j];
	ret=int_qi(i,xi);
	return ret;
}



void globalvar::findneqk(int i)
{
	int dim,j;
	int *IPIV=new int[dgace];
	int INFO;
	int NRHS;
	double ner[dgace];
	double re[dgace][dgace],xi;

	for(int m=0;m<dgace;m++)
	{
		ner[m]=reconstructqe(i,m);
		xi=Quad_Points[m];
		for(j=0;j<dgace;j++)
			re[j][m]=int_P(xi,j);
	}

	ner[0]=qe[2*i-1];
	ner[1]=qe[2*i];

	dim=dgace;
	NRHS=1;
	dgesv_(&dim,&NRHS,&re[0][0],&dim,&IPIV[0],&ner[0],&dim,&INFO);
	if(INFO!=0)
		cout<<"QeReshape Linear system solve not correct: INFO = "<<INFO<<endl;

	for(j=0;j<dgace;j++)
		qe_k[i][j]=ner[j];

	delete [] IPIV;

}


double globalvar::reconstructqe(int i,int j)
{
	double ret,xi;
	xi=Quad_Points[j];
	ret=int_qe(i,xi);
	return ret;
}
