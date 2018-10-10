#include <cmath>
#include <vector>
using namespace std;

void vprod(double *a, double *b, double *c)   // cross product C = A x B 
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}
void vprod(double *a, vector<double> &b, double *c)   // cross product C = A x B 
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}


double  vMag(double *v, int n)
{
	double sum=0.0;
	for(int i=0;i<n;i++)
	    sum+=v[i]*v[i];
    return sqrt(sum);		
}
double  vMag(vector<double>v, int n)
{
	double sum=0.0;
	for(int i=0;i<n;i++)
	    sum+=v[i]*v[i];
    return sqrt(sum);		
}

void vproject(double *v, double *n, double *t, int itype)
{
	// calcs the project of vector V on the plane 
//whose normal vector is N
    double p[3],_n[3];
    vprod(n, v, p) ; // p = N x V
	
	for (int k= 0; k < 3; k++)  //_n is -n
	  { _n[k]=-n[k];}
    vprod(_n, p, t) ; // T =-N x P
    if( (itype == 0) || ( vMag(v,3) < 0.001)|| ( vMag(v,3) > (vMag(t,3) * 1000)) ) return;
	   for (int k= 0; k < 3; k++)
	   t[k] = t[k] / (vMag(t,3)/vMag(v,3)) ;
}

void vproject(vector<double> v, double *n, double *t, int itype)
{
	// calcs the project of vector V on the plane 
//whose normal vector is N
    double p[3],_n[3];
    vprod(n, v, p) ; // p = N x V
	
	for (int k= 0; k < 3; k++)  //_n is -n
	  { _n[k]=-n[k];}
    vprod(_n, p, t) ; // T =-N x P
    if( (itype == 0) || ( vMag(v,3) < 0.001) ||( vMag(v,3) > (vMag(t,3) * 1000)) ) return;
	   for (int k= 0; k < 3; k++)
	   t[k] = t[k] / (vMag(t,3)/vMag(v,3)) ;
}


