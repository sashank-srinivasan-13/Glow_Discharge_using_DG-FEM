
#include "fluxes.h"
using namespace std;

double globalvar::Fni(int i)
{
	double f=nirk[2*i+1];
	return f;
}

double globalvar::Fne(int i)
{
	double f;
	f=(nerk[2*i]+nerk[2*i+1])/2.0;
	return f;
}


double globalvar::Qni(int i)
{
	double f;
	f=(2-positivity_on)*qi[2*i]+positivity_on*qi[2*i+1]+beta_i[i]*(nirk[2*i+1]-nirk[2*i]);
	f=f/2.0;
	return f;
}

double globalvar::Qne(int i)
{
	double f;
	f=(qe[2*i]+qe[2*i+1]+beta_e[i]*(nerk[2*i+1]-nerk[2*i]))/2.0;
	return f;
}

//--------------------------------------------------------------

double globalvar::CFni(int i, int j)
{
	double f;

	f=nirk[2*i]*Ef[2*i]+nirk[2*i+1]*Ef[2*i+1]-maxgradni[i]*(nirk[2*i+1]-nirk[2*i]);
	f=f/2.0;

	return f;
}

double globalvar::CFne(int i, int j)
{
	double f;

	f=nerk[2*i]*Ef[2*i]+nerk[2*i+1]*Ef[2*i+1]+maxgradne[i]*(nerk[2*i+1]-nerk[2*i]);
	f=f/2.0;

	return f;
}








