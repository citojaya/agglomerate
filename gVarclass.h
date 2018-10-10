#include <vector>
using namespace std;

#ifndef GVACLASS_H
#define GVACLASS_H


// Box parameters class
class CONTAINER
{   public:

	   double Rxy,th, bh;
       double gRxy, gth, gbh;
  	   double thv[3], bhv[3];
       int nWall;
       int nc;
	   void loadData(ifstream&);
};
//particle information (including wall)
class SIMULATION
{    
    public:	    
	    double cTime, chkTime, simTime,rlxtime;
        double dt, dtFactor;
        double rlunit, rmunit, rfunit, rsunit, rtunit, rvunit;
        double fdtime, tptime, anatime, dumptime;
        double fdpert, tppert, anapert, dumppert;
        double ctrF;
        double rmx;

		void loadData(ifstream&);
};
// --- liquid parameters ---
class LIQUIDTYPE
{    
    public:	    
	    double lqV, theta, brknmx;
        double sTension,cGapMn;
		double layer, cGapMx;
		int flag;

		void loadData(ifstream&);
};
// yield stress (plastic contact)
class YIELDTYPE
{    
    public:	    
	    double rstCefpw, rstCefpp;
        double stress,force,dsp;
        double stiffK;

		void loadData(ifstream&);
};
// parameters for compaction

class  COMPACTION
{    
public:
	    int stage;
        double strain, strain_mx;
		double l_org, dsp, dsp_mx;
		double stress, stress_mx;
		double zmn, zmx;
		void loadData(ifstream&);
};
class  MAT
{
public:
	    int np;
        double dia;
        double density;
        double emod, ymod, pois,yldp;
        double sfrc, rfrc;
        double dmpn,ha;
        void loadData(ifstream&);
};
class pInfo
{
	public:
        double upfc[3], upfv[3];     //forces from upward particles 
        double eg1, eg2,compress ;          //compress force
        int cn, upcn;            //total and up contact number
};
class IMPACT
{   public:
	double pos[3], vel[3];
    double vMag, angle;
}  ; 

const int mxnp = 10;
class ZONETYPE
	{
		public:
		double x, y, z;
        double velocity[3], force[3];
        double porosity, compress;
        int np;
        int pIndx[mxnp];
	};

#endif