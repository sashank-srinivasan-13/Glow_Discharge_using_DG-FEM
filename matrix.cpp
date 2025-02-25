#include "matrix.h"

using namespace std;



void globalvar::matrix_assemble()
{

	double phi_m_phi[25];
	double phi_p_phi[25];
	double q_m_q[25];
	double q_p_q[25];
	double P=pen;
	double Pd=pen_bound;
	double D=dx;

	phi_m_phi[0]=0.0;
	phi_m_phi[1]=0.0;
	phi_m_phi[2]=0.0;
	phi_m_phi[3]=0.0;
	phi_m_phi[4]=0.0;
	phi_m_phi[5]=0.0;
	phi_m_phi[6]=0.0;
	phi_m_phi[7]=0.0;
	phi_m_phi[8]=0.0;
	phi_m_phi[9]=0.0;

	phi_m_phi[10]=-1.0;
	phi_m_phi[11]=1.0;
	phi_m_phi[12]=-1.0;
	phi_m_phi[13]=1.0;
	phi_m_phi[14]=-1.0;
	phi_m_phi[15]=0.0;
	phi_m_phi[16]=0.0;
	phi_m_phi[17]=0.0;
	phi_m_phi[18]=0.0;
	phi_m_phi[19]=0.0;

	phi_m_phi[20]=1.0;
	phi_m_phi[21]=-1.0;
	phi_m_phi[22]=1.0;
	phi_m_phi[23]=-1.0;
	phi_m_phi[24]=1.0;


	//-----------------------------------------------------------

	phi_p_phi[0]=0.0;
	phi_p_phi[1]=0.0;
	phi_p_phi[2]=0.0;
	phi_p_phi[3]=0.0;
	phi_p_phi[4]=0.0;
	phi_p_phi[5]=0.0;
	phi_p_phi[6]=0.0;
	phi_p_phi[7]=0.0;
	phi_p_phi[8]=0.0;
	phi_p_phi[9]=0.0;

	phi_p_phi[10+0]=1.0;
	phi_p_phi[10+1]=-1.0;
	phi_p_phi[10+2]=1.0;
	phi_p_phi[10+3]=-1.0;
	phi_p_phi[10+4]=1.0;
	phi_p_phi[10+5]=0.0;
	phi_p_phi[10+6]=0.0;
	phi_p_phi[10+7]=0.0;
	phi_p_phi[10+8]=0.0;
	phi_p_phi[10+9]=0.0;

	phi_p_phi[10+10]=1.0;
	phi_p_phi[10+11]=-1.0;
	phi_p_phi[10+12]=1.0;
	phi_p_phi[10+13]=-1.0;
	phi_p_phi[10+14]=1.0;

	//----------------------------------------------------------------

	q_m_q[10-10]=P;
	q_m_q[10-9]=P;
	q_m_q[10-8]=P;
	q_m_q[10-7]=P;
	q_m_q[10-6]=P;
	q_m_q[10-5]=-1.0;
	q_m_q[10-4]=-1.0;
	q_m_q[10-3]=-1.0;
	q_m_q[10-2]=-1.0;
	q_m_q[10-1]=-1.0;

	q_m_q[10+0]=-2.0*P;
	q_m_q[10+1]=0.0;
	q_m_q[10+2]=-2.0*P;
	q_m_q[10+3]=0.0;
	q_m_q[10+4]=-2.0*P;
	q_m_q[10+5]=1.0;
	q_m_q[10+6]=1.0;
	q_m_q[10+7]=1.0;
	q_m_q[10+8]=1.0;
	q_m_q[10+9]=1.0;

	q_m_q[10+10]=P;
	q_m_q[10+11]=-P;
	q_m_q[10+12]=P;
	q_m_q[10+13]=-P;
	q_m_q[10+14]=P;



	//----------------------------------------------------------------


	q_p_q[10-10]=-P;
	q_p_q[10-9]=-P;
	q_p_q[10-8]=-P;
	q_p_q[10-7]=-P;
	q_p_q[10-6]=-P;
	q_p_q[10-5]=1.0;
	q_p_q[10-4]=1.0;
	q_p_q[10-3]=1.0;
	q_p_q[10-2]=1.0;
	q_p_q[10-1]=1.0;

	q_p_q[10+0]=0.0;
	q_p_q[10+1]=-2.0*P;
	q_p_q[10+2]=0.0;
	q_p_q[10+3]=-2.0*P;
	q_p_q[10+4]=0.0;
	q_p_q[10+5]=1.0;
	q_p_q[10+6]=1.0;
	q_p_q[10+7]=1.0;
	q_p_q[10+8]=1.0;
	q_p_q[10+9]=1.0;

	q_p_q[10+10]=P;
	q_p_q[10+11]=-P;
	q_p_q[10+12]=P;
	q_p_q[10+13]=-P;
	q_p_q[10+14]=P;


	//---------------------------------------------------------------------------------

	int i,j,k,start,factor;
	factor=10;


	for(i=0;i<10*(N-1);i++)
	{
		for(j=0;j<10*(N-1);j++)
			AA_prototype[i][j]=0.0;
	}

	i=1;
	start=factor*i-factor;

	for(k=0;k<5;k++)
		for(j=0;j<5;j++)
			AA_prototype[start+j+10][start+2*k]=pow(-1.0,j);

	AA_prototype[start+0][start+1]=-Pd-P;
	AA_prototype[start+1][start+1]=Pd-P;
	AA_prototype[start+2][start+1]=-Pd-P;
	AA_prototype[start+3][start+1]=Pd-P;
	AA_prototype[start+4][start+1]=-Pd-P;
	AA_prototype[start+6][start+1]=2.0;
	AA_prototype[start+8][start+1]=2.0;

	AA_prototype[start+0][start+5]=-Pd-P;
	AA_prototype[start+1][start+5]=Pd-P;
	AA_prototype[start+2][start+5]=-Pd-P;
	AA_prototype[start+3][start+5]=Pd-P;
	AA_prototype[start+4][start+5]=-Pd-P;
	AA_prototype[start+6][start+5]=2.0;
	AA_prototype[start+8][start+5]=2.0;

	AA_prototype[start+0][start+9]=-Pd-P;
	AA_prototype[start+1][start+9]=Pd-P;
	AA_prototype[start+2][start+9]=-Pd-P;
	AA_prototype[start+3][start+9]=Pd-P;
	AA_prototype[start+4][start+9]=-Pd-P;
	AA_prototype[start+6][start+9]=2.0;
	AA_prototype[start+8][start+9]=2.0;

	AA_prototype[start+0][start+3]=Pd-P;
	AA_prototype[start+1][start+3]=-Pd-P;
	AA_prototype[start+2][start+3]=Pd-P;
	AA_prototype[start+3][start+3]=-Pd-P;
	AA_prototype[start+4][start+3]=Pd-P;
	AA_prototype[start+5][start+3]=2.0;
	AA_prototype[start+7][start+3]=2.0;
	AA_prototype[start+9][start+3]=2.0;

	AA_prototype[start+0][start+7]=Pd-P;
	AA_prototype[start+1][start+7]=-Pd-P;
	AA_prototype[start+2][start+7]=Pd-P;
	AA_prototype[start+3][start+7]=-Pd-P;
	AA_prototype[start+4][start+7]=Pd-P;
	AA_prototype[start+5][start+7]=2.0;
	AA_prototype[start+7][start+7]=2.0;
	AA_prototype[start+9][start+7]=2.0;


	for(j=0;j<5;j++)
	{
		AA_prototype[start+10][start+2*j+1]=P;
		AA_prototype[start+11][start+2*j+1]=-P;
		AA_prototype[start+12][start+2*j+1]=P;
		AA_prototype[start+13][start+2*j+1]=-P;
		AA_prototype[start+14][start+2*j+1]=P;

	}


	AA_prototype[start+0][start+2]+=-2.0;
	AA_prototype[start+5][start+3]+=-2.0;
	AA_prototype[start+1][start+4]+=-2.0;
	AA_prototype[start+6][start+5]+=-2.0;
	AA_prototype[start+0][start+6]+=-2.0;
	AA_prototype[start+2][start+6]+=-2.0;
	AA_prototype[start+5][start+7]+=-2.0;
	AA_prototype[start+7][start+7]+=-2.0;
	AA_prototype[start+1][start+8]+=-2.0;
	AA_prototype[start+3][start+8]+=-2.0;
	AA_prototype[start+6][start+9]+=-2.0;
	AA_prototype[start+8][start+9]+=-2.0;



	AA_prototype[start+5][start]=AA_prototype[start+5][start]-D;
	AA_prototype[start+6][start+2]=AA_prototype[start+6][start+2]-D/3.0;
	AA_prototype[start+7][start+4]=AA_prototype[start+7][start+4]-D/5.0;
	AA_prototype[start+8][start+6]=AA_prototype[start+8][start+6]-D/7.0;
	AA_prototype[start+9][start+8]=AA_prototype[start+9][start+8]-D/9.0;

	for(i=2;i<N-1;i++)
	{
		start=factor*i-factor;

		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start]=phi_m_phi[j];
		}
		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start+1]=q_m_q[j];
		}
		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start+2]=phi_p_phi[j];
		}
		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start+3]=q_p_q[j];
		}
		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start+4]=phi_m_phi[j];
		}
		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start+5]=q_m_q[j];
		}
		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start+6]=phi_p_phi[j];
		}
		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start+7]=q_p_q[j];
		}
		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start+8]=phi_m_phi[j];
		}
		for(j=0;j<25;j++)
		{
			AA_prototype[start-10+j][start+9]=q_m_q[j];
		}


		AA_prototype[start+0][start+2]+=-2.0;

		AA_prototype[start+5][start+3]+=-2.0;

		AA_prototype[start+1][start+4]+=-2.0;

		AA_prototype[start+6][start+5]+=-2.0;

		AA_prototype[start+0][start+6]+=-2.0;
		AA_prototype[start+2][start+6]+=-2.0;

		AA_prototype[start+5][start+7]+=-2.0;
		AA_prototype[start+7][start+7]+=-2.0;

		AA_prototype[start+1][start+8]+=-2.0;
		AA_prototype[start+3][start+8]+=-2.0;

		AA_prototype[start+6][start+9]+=-2.0;
		AA_prototype[start+8][start+9]+=-2.0;




		AA_prototype[start+5][start]=AA_prototype[start+5][start]-D;
		AA_prototype[start+6][start+2]=AA_prototype[start+6][start+2]-D/3.0;
		AA_prototype[start+7][start+4]=AA_prototype[start+7][start+4]-D/5.0;
		AA_prototype[start+8][start+6]=AA_prototype[start+8][start+6]-D/7.0;
		AA_prototype[start+9][start+8]=AA_prototype[start+9][start+8]-D/9.0;
	}


	i=N-1;
	start=factor*i-factor;

	for(k=0;k<5;k++)
		for(j=0;j<5;j++)
			AA_prototype[start+j][start+2*k]=pow(-1.0,j+k+1);

	for(j=0;j<3;j++)
	{
		for(k=0;k<5;k++)
		{
			AA_prototype[start-10+k][start+4*j+1]=P;
			AA_prototype[start-5+k][start+4*j+1]=-1.0;

		}
	}

	for(j=0;j<2;j++)
	{
		for(k=0;k<5;k++)
		{
			AA_prototype[start-10+k][start+4*j+3]=-P;
			AA_prototype[start-5+k][start+4*j+3]=1.0;

		}
	}

	for(j=0;j<3;j++)
	{
		AA_prototype[start+0][start+4*j+1]=-Pd-P;
		AA_prototype[start+1][start+4*j+1]=-Pd+P;
		AA_prototype[start+2][start+4*j+1]=-Pd-P;
		AA_prototype[start+3][start+4*j+1]=-Pd+P;
		AA_prototype[start+4][start+4*j+1]=-Pd-P;
	}

	for(j=0;j<2;j++)
	{
		AA_prototype[start+0][start+4*j+3]=-Pd+P;
		AA_prototype[start+1][start+4*j+3]=-Pd-P;
		AA_prototype[start+2][start+4*j+3]=-Pd+P;
		AA_prototype[start+3][start+4*j+3]=-Pd-P;
		AA_prototype[start+4][start+4*j+3]=-Pd+P;
	}

	for(j=0;j<5;j++)
		for(k=0;k<5;k++)
			AA_prototype[start+5+k][start+2*j+1]=1.0;




	AA_prototype[start+0][start+2]+=-2.0;
	AA_prototype[start+5][start+3]+=-2.0;
	AA_prototype[start+1][start+4]+=-2.0;
	AA_prototype[start+6][start+5]+=-2.0;
	AA_prototype[start+0][start+6]+=-2.0;
	AA_prototype[start+2][start+6]+=-2.0;
	AA_prototype[start+5][start+7]+=-2.0;
	AA_prototype[start+7][start+7]+=-2.0;
	AA_prototype[start+1][start+8]+=-2.0;
	AA_prototype[start+3][start+8]+=-2.0;
	AA_prototype[start+6][start+9]+=-2.0;
	AA_prototype[start+8][start+9]+=-2.0;



	AA_prototype[start+5][start]=AA_prototype[start+5][start]-D;
	AA_prototype[start+6][start+2]=AA_prototype[start+6][start+2]-D/3.0;
	AA_prototype[start+7][start+4]=AA_prototype[start+7][start+4]-D/5.0;
	AA_prototype[start+8][start+6]=AA_prototype[start+8][start+6]-D/7.0;
	AA_prototype[start+9][start+8]=AA_prototype[start+9][start+8]-D/9.0;

	generate_matrix_inverse();
}




void globalvar::generate_matrix_inverse()
{
	int i,j;
	int order=dgacp;
	int R= 2*dgacp*(N-1);
	double *AA = new double[R*R];

	int count1=0,count2=0;
	int i7=0,j7=0;

	for(i=0;i<R;i++)
	{
		j7=0;
		count2=0;
		for(j=0;j<R;j++)
		{
			AA[j+R*i]=AA_prototype[i7][j7];
			count2++;
			if(count2%(2*order)==0)
			{
				j7=j7+11-2*order;
				count2=0;
			}
			else
			{
				j7=j7+1;
			}
		}
		count1++;
		if(count1%(order)==0)
		{
			i7=i7+6-order;
			count1=0;
		}
		else
		{
			i7=i7+1;
		}
	}


	int lwork=R*R;
	double *work;
	work=new double [lwork*sizeof(double)];
	int *IPIV=new int[R];
	int INFO=0;


	dgetrf_(&R, &R, AA, &R, IPIV, &INFO);
	dgetri_(&R, AA, &R, &IPIV[0], &work[0], &lwork, &INFO);


	int leap=2*dgacp*(N1);
	int r=2*dgacp*(Nlocal-1);
	for(i=0;i<r;i++)
		for(j=0;j<R;j++)
			AA_inverse[i][j]=AA[i+leap+j*R];

	for(i=0;i< 10*(N-1);i++)
		delete [] AA_prototype[i];

	delete [] AA_prototype;
	delete [] AA;
	delete [] work;
	delete [] IPIV;

}
