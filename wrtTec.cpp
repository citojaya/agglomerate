#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "gVarExt.h"
using namespace std;

void wrtTec()
{
   void getZPor(); 
   void getRInf();
   void getfrq();
   void getGr();
   double ijPos[3],cPos[3];
   double X[2000],Y[2000],Z[2000];
   int bnum;

   //set the coordinates of the box  
			 bnum = 100 ; // bnum must less than 1000
            for (int i=0;i<2;i++)
			{
				X[i]= box.Rxy * cos ((double(i)/double(bnum))*2*PI);
				Y[i]= box.Rxy * sin ((double(i)/double(bnum))*2*PI);
				Z[i]= box.bh;
			}
				X[2]= box.Rxy * cos ((double(1)/double(bnum))*2*PI);
				Y[2]= box.Rxy * sin ((double(1)/double(bnum))*2*PI);
				Z[2]= box.th;
				X[3]= box.Rxy * cos ((double(0)/double(bnum))*2*PI);
				Y[3]= box.Rxy * sin ((double(0)/double(bnum))*2*PI);
				Z[3]= box.th;
		
			 for (int i=4;i<2*bnum-1;i=i+2)
			{
				X[i]= box.Rxy * cos ((double(i/2)/double(bnum))*2*PI);
				Y[i]= box.Rxy * sin ((double(i/2)/double(bnum))*2*PI);
				Z[i]= box.bh;
				X[i+1]= box.Rxy * cos ((double(i/2)/double(bnum))*2*PI);
				Y[i+1]= box.Rxy * sin ((double(i/2)/double(bnum))*2*PI);
				Z[i+1]= box.th;
			}
   /***********---- particle related information******/
    ofstream outfile;
    string  filename = genfile +"_ptec.dat";
    outfile.open(filename.c_str(),ios::out|ios::in|ios::ate);
       if(!outfile) 
     { 
		outfile.clear();
		outfile.open(filename.c_str(),ios::out);
		outfile<<"TITLE = \"PARTICLE FLOW\""<<endl;
		outfile<<"VARIABLES = X, Y, Z, VX, VY, VZ, FX, FY, FZ, WX, WY, WZ, DIA, IP"<<endl;
	    outfile<<"ZONE T=   "<<"\""<<"zhenbo"<<"\""<<",N="<<2*bnum<<",E="<<bnum<<", ZONETYPE=FEQuadrilateral"<<endl; 
	    outfile<<"DATAPACKING=POINT"<<endl; 
        for (int i=0;i<2*bnum;i++)
			{
				outfile<<X[i]<<" "<<Y[i]<<" "<<Z[i]<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "
				<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;
			}
				outfile<<4<<" "<<3<<" "<<2<<" "<<1<<endl;
				outfile<<3<<" "<<6<<" "<<5<<" "<<2<<endl;
		for (int i=6;i<2*bnum-1;i=i+2)
			{
				outfile<<i<<" "<<i+2<<" "<<i+1<<" "<<i-1<<endl;

			}
				outfile<<2*bnum<<" "<<4<<" "<<1<<" "<<2*bnum-1<<endl;
     } 
    outfile<<setprecision(5);

	outfile<<"ZONE T=   "<<"\""<<dem.cTime<<"\""<<",I="<<np<<",J="<<1<<",F=Point"<<endl; 

    for(int ip=0;ip<np;ip++) 
		outfile<<setw(13)<<pPos[ip][0]<<setw(13)<<pPos[ip][1]<<setw(13)<<pPos[ip][2]<<setw(13)
		       <<pVel[ip][0]<<setw(13)<<pVel[ip][1]<<setw(13)<<pVel[ip][2]<<setw(13)
			   <<pForce[ip][0]<<setw(13)<<pForce[ip][1]<<setw(13)<<pForce[ip][2]<<setw(13)
               <<pAngv[ip][0]<<setw(13)<<pAngv[ip][1]<<setw(13)<<pAngv[ip][2]<<setw(13)
			   <<pDia[ip]<<setw(13)<<ip<<endl;
	outfile.close();


  /****** force information at contact point****************/
/*	filename = genfile +"_cftec.dat";
	outfile.open(filename.c_str(),ios::out|ios::in|ios::ate);
   if(!outfile) 
     { 
		outfile.clear();
		outfile.open(filename.c_str(),ios::out);
		outfile<<"TILE = \"FORCE AT CONTACT POINT\""<<endl;
		outfile<<"VARIABLES = X, Y, Z, FCNX, FCNY, FCNZ,FTX, FTY, FTZ, FV, DIA, IP, JP"<<endl;
     } 
    outfile<<"ZONE T=   "<<"\""<<dem.cTime<<"\""<<endl; 
   
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
		int ip=int(fList[ifc][7]);
		int jp=int(fList[ifc][8]);
		for (int k= 0; k < 3; k++)
		ijPos[k]=pPos[ip][k]-pPos[jp][k];

		for (int k= 0; k < 3; k++)
        cPos[k] = pPos[ip][k] - pDia[ip]/(pDia[ip] + pDia[jp]) * ijPos[k];

        outfile<<setw(13)<<cPos[0]<<setw(13)<<cPos[1]<<setw(13)<<cPos[2]<<setw(13)
			   <<fList[ifc][0]<<setw(13)<<fList[ifc][1]<<setw(13)<<fList[ifc][2]<<setw(13)<<fList[ifc][3]<<setw(13)
			   <<fList[ifc][4]<<setw(13)<<fList[ifc][5]<<setw(13)<<fList[ifc][6]<<setw(13)
			   <<0.5*(pDia[ip]+pDia[jp])<<setw(13)<<ip<<setw(13)<<jp<<endl;
	 }
	outfile.close();*/
/************----- force information on wall******************************/
 /*  filename = genfile +"_f2wall.dat";
   outfile.open(filename.c_str(),ios::out|ios::in|ios::ate);
   if(!outfile) 
     { 
		outfile.clear();
		outfile.open(filename.c_str(),ios::out);
		outfile<<"TILE = \"FORCE AT CONTACT POINT\""<<endl;
		outfile<<"VARIABLES = X, Y, Z, FCNX, FCNY, FCNZ,FTX, FTY, FTZ, FV, FC+V"<<endl;
     } 
   if(box.nc>0)
   {
   outfile<<"ZONE T=   "<<"\""<<dem.cTime<<"\""<<endl; 

	for(int i=0; i<box.nc;i++)
	{
		for(int j=0; j<10; j++)
		outfile<<setw(13)<<wpfList[i][j];
        outfile<<setw(13)<<wpfList[i][5]+wpfList[i][9]<<endl;
	}
   }
   outfile.close();*/
	 
  /* filename = genfile +"_f2wall.txt";
   outfile.open(filename.c_str(),ios::out|ios::in|ios::ate);
   if(!outfile) 
     { 
		outfile.clear();
		outfile.open(filename.c_str(),ios::out);
		outfile<<"TILE = \"FORCE OF WALL\""<<endl;
		outfile<<"VARIABLES = I,J,FCNX, FCNY, FCNZ"<<endl;
     } 
    outfile<<"ZONE T=   "<<"\""<<dem.cTime<<"\""<<endl; 
	for(int i=0; i<box.nWall;i++)
		for(int j=0; j<2; j++)
     outfile<<i<<" "<<j<<" "<<setw(13)<<wForce[i][j][0]<<setw(13)<<wForce[i][j][1]<<setw(13)<<wForce[i][j][2]<<endl;
	 outfile.close();*/


	 // getRInf();
/************  ---- zone related information ----****************/    
/*if (dem.cTime <= 25*dem.tppert)  // calculate 20 times
 {
	 getZPor();
	 filename = genfile +"_ztec.dat";
	outfile.open(filename.c_str(),ios::out|ios::in|ios::ate);
   if(!outfile) 
     { 
		outfile.clear();
		outfile.open(filename.c_str(),ios::out);
		outfile<<"TILE = \"FIELD INFORMATION\""<<endl;
		outfile<<"VARIABLES = I, J, K, POROSITY"<<endl;
     } 
     outfile<<"ZONE T=   "<<"\""<<dem.cTime<<"\""<<",I="<<nxzone<<",J="<<nyzone<<",K="<<nzzone<<endl; 
        for(int k=0;k<nzzone;k++)
	    {
			for(int j=0; j<nyzone;j++)
			{    
				for (int i=0;i<nxzone;i++)
				{
					izone=zone[i][j][k];
					outfile<<izone.x<<" "<<izone.y<<" "<<izone.z<<" "<<izone.porosity<<endl;
				}
			}
		}
      deleteArray(zone,nxzone,nyzone);//delete the array of zone 
  }*/

   /*******---- cn distribution ---************/
  /*  getfrq();
    filename = genfile +"_cntec.dat";
	outfile.open(filename.c_str(),ios::out|ios::in|ios::ate);
    if(!outfile) 
     { 
		outfile.clear();
		outfile.open(filename.c_str(),ios::out);
		outfile<<"TILE = \"COORDINATION DISTRIBUTION\""<<endl;
		outfile<<"VARIABLES = CN, FREQUENCY"<<endl;
     } 
    outfile<<"ZONE T=   "<<"\""<<dem.cTime<<"\""<<endl; 
    const int  frqBin = 15;
    extern double freq[frqBin+1];
	for (int i=0; i<= frqBin; i++)
    outfile<<i<<" "<<freq[i]<<endl;
	outfile.close(); */

    /*---- rdf function ----*/
  /*  getGr();
    filename = genfile +"_grtec.dat";
	outfile.open(filename.c_str(),ios::out|ios::in|ios::ate);
    if(!outfile) 
     { 
		outfile.clear();
		outfile.open(filename.c_str(),ios::out);
		outfile<<"TILE = \"GR FUNCTION\""<<endl;
		outfile<<"VARIABLES = RADI, GR"<<endl;
     } 
    outfile<<"ZONE T=   "<<"\""<<dem.cTime<<"\""<<endl; 
    const int grBin = 250;
    extern double gr[grBin][2];

	 for(int i=0; i< grBin; i++) 
	 {
		 for(int j=0; j< 2; j++) 
		 {outfile<<gr[i][j]<<" ";}
		 outfile<<endl;
	 }

	outfile.close();*/

}
