#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "gVarExt.h"
using namespace std;
void analyz()
{
	double  avgVel,vp1, vp2, czmn, czmx , vol, lz ,cn, strain;
	void outfile(int) ;

	avgVel = 0.0 ;vp1 = 0.0; vp2 = 0.0;
	lz = compact.zmx - compact.zmn ;
	czmn = compact.zmn + 0.25 * lz ;
	czmx = compact.zmn + 0.75 * lz;
	vol = lz * (PI*pow(box.Rxy,2));

	  for(int i=0; i<np ; i++)
        {
			avgVel += sqrt(pow(pVel[i][0],2)+pow(pVel[i][1],2)+pow(pVel[i][2],2));
			vp1 += PI * pow(pDia[i],3.0) / 6.0 ;
			if (pPos[i][2] > czmn && pPos[i][2] < czmx) 
				vp2 += PI * pow(pDia[i],3.0)/6.0 ;
		}
	 compact.stress = compact.stress / (PI*pow(box.Rxy,2));
	 tensile_stress = tensile_stress /vol /3.0  ; 
	 avgVel = avgVel / np ;
	
	/* cn =-1; //contact number
	for(int ip=0;ip<fList.size();ip++) 
	 {
	  if( cn < fList[ip][9])
		  cn = fList[ip][9]; 
	 }  
	 cn++;   //contact number from 0 start */
	
	 if (compact.l_org == 0.0)
		 strain = 0.0 ;
	 else
		 strain = compact.dsp/compact.l_org;
	 
	 // output the file
	 ofstream ofile;
     string  filename = genfile +"_analyz.txt";
     ofile.open(filename.c_str(),ios::out|ios::in|ios::ate);     
	 if(!ofile) 
     { 
		ofile.clear();
		ofile.open(filename.c_str(),ios::out);
		ofile<<setw(13)<<"Time"<<setw(13)<<"Solid_Fra"<<setw(13)<<"Solid_FraM"<<setw(13)<<"Strain"<<setw(13)<<"Cpress(Kpa)"<<setw(13)<<"Tstress(Kpa)"
			<<setw(13)<<"S_height"<<setw(13)<<"Avg_v(m/s)"<<setw(13)<<"ContactNum"<<setw(13)<<"overlap_max"<<endl;
     } 

	/* force of normal and tensile p-p
    double fn=0;
	double fc=0;
    double fList_max; //find the max value of fList[][9]
	fList_max=-1;
	for(int ip=0;ip<fList.size();ip++) 
	 {
	  if(fList_max < fList[ip][9])
		 fList_max=fList[ip][9]; 
	 }  
	int nfc = int(fList_max);

    for(int ifc=0;ifc<=nfc;ifc++) 
	 {
		fn+=pow(pow(fList[ifc][0],2)+pow(fList[ifc][1],2)+pow(fList[ifc][2],2),0.5);
		fc+=pow(pow(fList[ifc][3],2)+pow(fList[ifc][4],2)+pow(fList[ifc][5],2),0.5); 
	 }
	fn=fn/(nfc+1);
	fc=fc/(nfc+1);*/
	  ofile<<setprecision(5);
	  ofile<<setw(13)<<dem.cTime<<setw(13)<<vp1/vol<<setw(13)<<2.0*vp2/vol<<setw(13)<<strain<<setw(13)<<compact.stress * dem.rsunit/1.0e3 <<setw(13)<<
		  tensile_stress*dem.rsunit/1e3<<setw(13)<<lz<<setw(13)<<avgVel*dem.rvunit<<setw(13)<<2*double(tensile_cnv)/double(np)<<setw(13)<<mxOlp<<endl; 
	  ofile.close();
	  // if average velocity less than 1e-6 , then stop packing
	  if ( compact.stage == 0 && avgVel*dem.rvunit < 1e-7 )
		{
	     	dump();
			outfile(1);
			exit(1);
		}

}

