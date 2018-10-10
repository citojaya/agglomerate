#include <iostream> 
#include <fstream>
#include <string>
#include<cmath>
#include<iomanip>
#include "gVarExt.h"
using namespace std;


void outfile(int flag)
{

     if (flag == 0) 
	 {
        ofstream outfile;
		string filename=genfile+".out";
		outfile.open(filename.c_str(),ios::out|ios::app);
       	if(!outfile)
		{cout<<"Unable to open the input file: "<<filename<<"\n"<<"Please try again."<<endl;}
        outfile<<setiosflags(ios::scientific);   
        outfile<<"-----------------------------------------"<<endl;
        outfile<< "REDUCED TIME UNIT = "<<dem.rtunit<< "  sec"<<endl;
        outfile<< "REDUCED FORCE UNIT = "<< dem.rfunit<< "  N"<<endl;
        outfile<< "REDUCED VELOCITY UNIT = "<< dem.rvunit<< "  m/s"<<endl;
        outfile<< "REDUCED EMODULUS = "<< pMat[0].emod<<endl;
        outfile<< "MAXIMUM FVW = "<< pMat[1].ha / pow(vGapMn,2)/pow(Min(pDia,npTotal),2)/4.0<< "  mg"<<endl;
        outfile<< "MAXIMUM VGAP = "<<vGapMx/Min(pDia,npTotal)<< "  d"<<endl;
        outfile<< "MAXIMUM BGAP = "<<bGapmx/Min(pDia,npTotal)<< "  d"<<endl;
		outfile<< "MAXIMUM lGAP = "<<lqd.brknmx/Min(pDia,npTotal)<< "  d"<<endl;
		outfile<< "MAXIMUM cutGAP = "<<cutGap<< "  d"<<endl;
        outfile<< "TIME STEP = "<< dem.dt * dem.rtunit<< "  sec"<<endl;
        outfile<< "P-P & P-W RESTITUTION CO (Plastic) = "<<  yld.rstCefpp<< "  "<< yld.rstCefpw<<endl;
        outfile<<"-----------------------------------------"<<endl;
		outfile.close();
	 }
        else if (flag == 1) 
      {
        ofstream outfile;
		string filename=genfile+".out";
		outfile.open(filename.c_str(),ios::out|ios::app);
        pclock();//cpu time
		outfile<<"STOP AT DEM TIME ="<<dem.cTime<<" CPU TIME = "<<etime/60<<" min"<<endl;
        outfile<<"-----------------------------------------"<<endl;
		outfile.close();

		/*string filenametxt=genfile+".txt";
		outfile.open(filenametxt.c_str(),ios::out|ios::app);
        outfile<<np<<endl;
		outfile<<box.hflx<<endl;
		for(int i=0;i<np;i++)
		{outfile<<pPos[i][0]<<"  "<<pPos[i][0]<<"  "<<pPos[i][0]<<"  "<<pDia[i]<<endl;}
		
		outfile.close();*/
		
	   }



}
        
        




