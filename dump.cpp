#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "gVarExt.h"
using namespace std;
void dump()
{
    ofstream outfile;
    string filename=genfile+".dump";
	outfile.open(filename.c_str(),ios::out);
   /****write out np,iter, position and velocity, displacement*****/
   
    outfile<<"DEM"<<endl;
    outfile<<dem.cTime<<endl;
    
    outfile<<"COMPACT"<<endl;
    outfile<<compact.stage<<endl;
	outfile<<compact.l_org<<endl;//
    
    outfile<<"CONTAINER"<<endl;
	outfile<<box.th<<" "<<box.bh<<" "<<endl;
	for(int j=0;j<3;j++)
		outfile<<box.thv[j]<<" "<<box.bhv[j]<<" ";
		outfile<<endl;
    outfile<<"BONDINITIAL"<<endl;
	    outfile<<b_flag<<endl;
	outfile<<"PARTICLE"<<endl;
    outfile<<np<<endl;

    outfile<<setiosflags(ios::scientific);   
    outfile<<setprecision(16);               

		for(int i=0;i< np;i++) 
		{
			for(int j=0;j<3;j++)
			{
		    outfile<<pPos[i][j]<<" "<<pVel[i][j]<<" "<<pAngp[i][j]<<" "<<pAngv[i][j]<<" ";

	           for(int k=0;k<box.nWall;k++)
			   {outfile<<dispt_pw[i][k][j]<<" ";}
			}

		     outfile<<pDia[i]<<" "<<pType[i]<<" "<<nPoint[i]<<endl;
			
			
		}

    outfile<<"NEIGHBOR"<<endl;
    outfile<<nne<<endl;

		for(int i=0;i<nne;i++) 
		{
			outfile<< nList[i]<<" "<<rbond[i]<<" "<<lqdbrg[i]<<" "<<dispAn_bond[i]<<" "<<ch_nDsp_pp[i]<<" "<<ch_nDsp_mx_pp[i]<<" "<<ch_plstDfm_pp[i]<<" "<<ch_plstRad_pp[i]<<" "<<initalGap_bond[i]<<" ";
            for(int j=0;j<3;j++)
				{outfile<<dispt_pp[i][j]<<" "<<dispAt_bond[i][j]<<" "<<dispt_bond[i][j]<<" ";}
			outfile<<endl;
		}
	outfile.close();


    cout<<"dumping finished"<<endl;




}
