
#include "num_int.h"

#include <boost/lambda/lambda.hpp>
#include <boost/math/special_functions/expint.hpp>

using namespace std;
using namespace boost::lambda;
typedef std::istream_iterator<int> in;

double quad_points[10]={-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,-0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,-0.9739065285171717,0.9739065285171717};
double quad_weights[10]={0.2955242247147529,0.2955242247147529,0.2692667193099963,0.2692667193099963,0.2190863625159820,0.2190863625159820,0.1494513491505806,0.1494513491505806,0.0666713443086881,0.0666713443086881};
int order=10;



double newpoints[128]={
		0.0486909570091397,-0.0243502926634244,0.0486909570091397,0.0243502926634244,0.0485754674415034,-0.0729931217877990,0.0485754674415034,0.0729931217877990,0.0483447622348030,
		-0.1214628192961206,0.0483447622348030,0.1214628192961206,0.0479993885964583,-0.1696444204239928,0.0479993885964583,0.1696444204239928,0.0475401657148303,-0.2174236437400071,
		0.0475401657148303,0.2174236437400071,0.0469681828162100,-0.2646871622087674,0.0469681828162100,0.2646871622087674,0.0462847965813144,-0.3113228719902110,0.0462847965813144,
		0.3113228719902110,0.0454916279274181,-0.3572201583376681,0.0454916279274181,0.3572201583376681,0.0445905581637566,-0.4022701579639916,0.0445905581637566,0.4022701579639916,
		0.0435837245293235,-0.4463660172534641,0.0435837245293235,0.4463660172534641,0.0424735151236536,-0.4894031457070530,0.0424735151236536,0.4894031457070530,0.0412625632426235,
		-0.5312794640198946,0.0412625632426235,0.5312794640198946,0.0399537411327203,-0.5718956462026340,0.0399537411327203,0.5718956462026340,0.0385501531786156,-0.6111553551723933,
		0.0385501531786156,0.6111553551723933,0.0370551285402400,-0.6489654712546573,0.0370551285402400,0.6489654712546573,0.0354722132568824,-0.6852363130542333,0.0354722132568824,
		0.6852363130542333,0.0338051618371416,-0.7198818501716109,0.0338051618371416,0.7198818501716109,0.0320579283548516,-0.7528199072605319,0.0320579283548516,0.7528199072605319,
		0.0302346570724025,-0.7839723589433414,0.0302346570724025,0.7839723589433414,0.0283396726142595,-0.8132653151227975,0.0283396726142595,0.8132653151227975,0.0263774697150547,
		-0.8406292962525803,0.0263774697150547,0.8406292962525803,0.0243527025687109,-0.8659993981540928,0.0243527025687109,0.8659993981540928,0.0222701738083833,-0.8893154459951141,
		0.0222701738083833,0.8893154459951141,0.0201348231535302,-0.9105221370785028,0.0201348231535302,0.9105221370785028,0.0179517157756973,-0.9295691721319396,0.0179517157756973,
		0.9295691721319396,0.0157260304760247,-0.9464113748584028,0.0157260304760247,0.9464113748584028,0.0134630478967186,-0.9610087996520538,0.0134630478967186,0.9610087996520538,
		0.0111681394601311,-0.9733268277899110,0.0111681394601311,0.9733268277899110,0.0088467598263639,-0.9833362538846260,0.0088467598263639,0.9833362538846260,0.0065044579689784,
		-0.9910133714767443,0.0065044579689784,0.9910133714767443,0.0041470332605625,-0.9963401167719553,0.0041470332605625,0.9963401167719553,0.0017832807216964,-0.9993050417357722,
		0.0017832807216964,0.9993050417357722
};

int order2=64;

double globalvar::int_P(double xi, int k)
{
	double value=1.0;

	if(k==0)
		value=1.0;

	else if(k==1)
		value=xi;

	else if(k==2)
		value=(3.0*xi*xi-1.0)/2.0;

	else if(k==3)
		value=(5.0*pow(xi,3.0)-3.0*xi)/2.0;

	else if(k==4)
		value=(35.0*pow(xi,4.0)-30.0*pow(xi,2.0)+3.0)/8.0;

	else if(k==5)
		value=(63.0*pow(xi,5.0)-70.0*pow(xi,3.0)+15.0*xi)/8.0;

	return value;
}

double globalvar::int_gradP(double xi, int k)
{
	double value=0.0;

	if(k==0)
		value=0.0;

	else if(k==1)
		value=1.0;

	else if(k==2)
		value=3.0*xi;

	else if(k==3)
		value=(15.0*pow(xi,2.0)-3.0)/2.0;

	else if(k==4)
		value=(140*pow(xi,3.0)-60.0*pow(xi,1.0))/8.0;

	else if(k==5)
		value=(315.0*pow(xi,4.0)-210.0*pow(xi,2.0)+15.0)/8.0;

	return value;
}

double globalvar::int_ni(int i, double xi)
{
	int j;
	double sum=0.0;

	for(j=0;j<dgaci;j++)
		sum=sum+ni_k_RK[i][j]*int_P(xi,j);

	return sum;
}


double globalvar::int_ne(int i, double xi)
{
	int j;
	double sum=0.0;

	for(j=0;j<dgace;j++)
		sum=sum+ne_k_RK[i][j]*int_P(xi,j);

	return sum;
}

double globalvar::int_newni(int i, double xi)
{
	int j;
	double sum=0.0;

	for(j=0;j<dgaci;j++)
		sum=sum+new_ni_k[i][j]*int_P(xi,j);

	return sum;
}


double globalvar::int_newne(int i, double xi)
{
	int j;
	double sum=0.0;

	for(j=0;j<dgace;j++)
		sum=sum+new_ne_k[i][j]*int_P(xi,j);

	return sum;
}

double globalvar::int_qi(int i, double xi)
{
	int j;
	double sum=0.0;

	for(j=0;j<dgaci;j++)
		sum=sum+qi_k[i][j]*int_P(xi,j);

	return sum;
}


double globalvar::int_qe(int i, double xi)
{
	int j;
	double sum=0.0;

	for(j=0;j<dgace;j++)
		sum=sum+qe_k[i][j]*int_P(xi,j);

	return sum;
}

double globalvar::int_E(int i, double xi)
{
	double sum=0.0;

	int j;
	for(j=0;j<dgacp;j++)
		sum=sum+qphi_k[i][j]*int_P(xi,j);

	return -1.0*sum;
}

double globalvar::int_phi(int i, double xi)
{
	double sum=0.0;

	int j;
	for(j=0;j<dgacp;j++)
		sum=sum+phi_k[i][j]*int_P(xi,j);

	return sum;
}


double globalvar::int_qphi(int i, double xi)
{
	double sum=0.0;

	int j;
	for(j=0;j<dgacp;j++)
		sum=sum+qphi_k[i][j]*int_P(xi,j);

	return sum;
}


double globalvar::int_gradE(int i, double xi)
{
	int j;
	double sum=0.0;

	for(j=0;j<dgacp;j++)
		sum=sum+qphi_k[i][j]*int_gradP(xi,j);

	return -1.0*sum;
}

double globalvar::Alpha(int i)
{
	double ret=0.0;
	double xi;
	for(int j=0;j<order;j++)
	{
		xi=quad_points[j];
		if(int_E(i,xi)==0.0)
			ret=ret+0.0;
		else
			ret=ret+quad_weights[j]*p*A*exp(-B*p/abs(int_E(i,xi)));
	}


	return ret/2.0;     //for average of numerical integration 
}


double globalvar::getalpha(int i, double xi)
{
	double ret;

	if(t>pow(10,-6))
	{
		double arg=-B*p/(qphi_k[i][0]+qphi_k[i][1]*xi);
		ret=-(B*p-qphi_k[i][1]*xi)/qphi_k[i][1]*exp(-B*p/(qphi_k[i][0]+qphi_k[i][1]*xi));
		ret=ret+B*p/qphi_k[i][1]*boost::math::expint(arg);
	}
	else
	{
		if (t<pow(10,-6))
			if(Ef[2*i+1]==0.0)
				ret=0.0;
			else
				ret=p*A*exp(-B*p/abs(Ef[2*i+1]));
		else
			if(int_E(i,xi)==0.0)
				ret=0.0;
			else
				ret=p*A*exp(-B*p/abs(int_E(i,xi)));
	}

	return ret;
}


double globalvar::aG(int i, int k)
{
	int j;
	double sum=0.0;
	double xi;
	for(j=0;j<order;j++)
	{
		xi=quad_points[j];
		sum=sum+quad_weights[j]*abs(mu_e*int_ne(i,xi)*int_E(i,xi))*int_P(xi,k);;
	}

	return sum*(getalpha(i,+1.0)-getalpha(i,-1.0))/2.0;
}

double globalvar::bnn(int i, int k)
{
	double sum=0.0, xi;
	int j;
	for(j=0;j<order;j++)
	{
		xi=quad_points[j];
		sum=sum+quad_weights[j]*int_ni(i,xi)*int_ne(i,xi)*int_P(xi,k);
	}

	sum=sum*beta;
	return sum;
}

double globalvar::nmuE_i(int i, int k)
{
	double sum=0.0,xi;
	int j;
	for(j=0;j<order;j++)
	{
		xi=quad_points[j];
		sum=sum+quad_weights[j]*int_ni(i,xi)*int_E(i,xi)*int_gradP(xi,k);
	}

	sum=sum*mu_i;
	return sum;
}


double globalvar::nmuE_e(int i, int k)
{
	double sum=0.0,xi;
	int j;
	for(j=0;j<order;j++)
	{
		xi=quad_points[j];
		sum=sum+quad_weights[j]*int_ne(i,xi)*int_E(i,xi)*int_gradP(xi,k);
	}

	sum=sum*mu_e;
	return sum;
}





