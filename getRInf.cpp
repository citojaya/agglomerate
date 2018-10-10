#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include "gVarExt.h"
using namespace std;
double pVol(double r1, double r2,double d);

void getRInf()
{ 
    const double  dr = 1.00001;
    double rad, rad1, ipMag;
    double  v;
    int ibin, nbin;
    double  **rdc, **rdcl;

	nbin = int(dem.rmx / dr) + 1;
    makeArray(rdc,nbin,4);
    makeArray(rdcl,nbin,4);

    //layer information

    for(int ip=0; ip <np; ip++)
	{
	    double pPos_oPos[3];
		for(int k=0; k<3; k++)
		{

		pPos_oPos[k] = pPos[ip][k] - oPos[k];
		}

  		ipMag = vMag(pPos_oPos,3);
        ibin = int(ipMag / dr);
        rdcl[ibin][0] += 1.0 ;                //n particles
        rdcl[ibin][1] += pInf[ip].cn ;         //n contacts
        v = 0.0;
       for(int ib=ibin-1 ; ib<=ibin+1 ; ib++)
	   {
		   if ((ib >= 0) && (ib < nbin-1))               //particle volume
		   {
             // rdcl[ib][2] =  rdcl[ib][2] + pVol((ib+1)*dr, 0.5*pDia[ip], ipMag) - v ;
              v = pVol((ib+1)*dr, 0.5*pDia[ip], ipMag) ;
		   }
	   }
	}
	// sphere information


    for (ibin=0; ibin<nbin ; ibin++)
	{        
		rad = (ibin+1) * dr;
        rad1 = ibin * dr;
        rdcl[ibin][3] = 4.0 / 3.0 * (pow(rad,3.0) - pow(rad1,3.0));   //layer volume

        for(int k=0; k<=ibin; k++)
		{
		 rdc[ibin][0] += rdcl[k][0];
         rdc[ibin][1] += rdcl[k][1];
         rdc[ibin][2] += rdcl[k][2];
         rdc[ibin][3] += rdcl[k][3];
		 }
	}

   //output
    ofstream outfile;
    string  filename = genfile +"_rdc.dat";
	outfile.open(filename.c_str(),ios::out|ios::in|ios::ate);
   if(!outfile) 
     { 
		outfile.clear();
		outfile.open(filename.c_str(),ios::out);
		outfile<<"VARIABLES = RADI, NP, CN, DENSITY, nplayer, cnlayer, densitylayer"<<endl;
     } 
    outfile<<"ZONE T=   "<<"\""<<dem.cTime<<"\""<<endl; 
	outfile<<setprecision(5)<<setiosflags(ios::showpoint)<<setiosflags(ios::fixed);

    for (ibin=0; ibin<nbin ; ibin++)
	{
		if(rdc[ibin][0] > 0.0)  rdc[ibin][1] = rdc[ibin][1] / rdc[ibin][0];
        if(rdcl[ibin][0]> 0.0) rdcl[ibin][1] = rdcl[ibin][1] / rdcl[ibin][0];
        outfile<<setw(13)<<ibin*dr<<setw(13)<<rdc[ibin][0]<<setw(13)<<rdc[ibin][1]<<setw(13)<<rdc[ibin][2]/rdc[ibin][3]<<setw(13)<<
                                rdcl[ibin][0]<<setw(13)<<rdcl[ibin][1]<<setw(13)<<rdcl[ibin][2]/rdcl[ibin][3]<<endl;
	}
	outfile.close();

}

double pVol(double r1, double r2,double d)

{   double pVol_value=0.0;
	double h1=0.0, h2=0.0;
    double d1=0.0, d2=0.0;
    double rr=0.0;

    if (d <= (r1-r2)) 
        pVol_value = 4.0/3.0 * pow(r2 , 3.0);
    else if (d >= (r1+r2)) 
        pVol_value = 0.0;
    else
	{ d1 = 0.5 * (d * d + r1 * r1 - r2 * r2) / d;
        h1 = r1 - d1;
        d2 = d - d1;
        h2 = r2 - d2;
        rr = r1 * r1 - d1 * d1;
        pVol_value = (h1 * (h1 * h1 + 3.0 * rr) + h2 * (h2 * h2 + 3.0 * rr)) / 6.0;
	}
	    return pVol_value;
}
