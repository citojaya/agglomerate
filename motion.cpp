#include <cmath>
#include "gVarExt.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void bldNList();
void feed();
void force(); 
void wrtTec();
void outfile(int flag);
void analyz();

void motion()
{
	double **pMov;
    double wMov[3],qMovMx;
    bool rBldNList;
	makeArray(pMov,npTotal,3);

/**********---- set to zero ---******/
    rBldNList =true;
    qMovMx = 0.5 * cutGap;

     // dem.rmx = 15.0;
/*******---- starting simulation ----***********************************************************************/
while (dem.cTime <= dem.simTime)
 {
       //generate new particles ----
	  if ((np < npTotal) &&(dem.fdtime <= dem.cTime))   
		   { 
			   feed();
               rBldNList = true;
               dem.fdtime = dem.fdtime + dem.fdpert;
	       }

     //***********---- build neighbor list ----
        if (rBldNList) 
	   	  {
			bldNList();
            rBldNList = false;
               for(int   i=0;i<npTotal;i++)  
		      { 
			    for(int   j=0;j<3;j++)   
				 {    
                   pMov[i][j] = 0.0;
                   wMov[j] = 0.0;
				 }
			   }
		   }
   /************---- update force ---************/
		force();
   /***************---- update position, velocity ----*/
        //dem.rmx = 0.0;
        for(int ip=0;ip < np;ip++)
		{
			//---- translational velocity and position
		   for(int   k=0;k<3;k++)   
		    {
			   pVel[ip][k] = pVel[ip][k] + pForce[ip][k] * dem.dt / pMass[ip];
               pPos[ip][k] = pPos[ip][k] + pVel[ip][k] * dem.dt;
		    }
           //--- shift with oPos

            //pPos[ip][k] = pPos[ip][k] - oPos[k]
           
            //---- angular velocity and position ----
		   for(int   k=0;k<3;k++)   
		   {pAngv[ip][k] += pMom[ip][k] * dem.dt / pInert[ip];}
            pRMom[ip]   = pRMom[ip] * dem.dt / pInert[ip];      // rolling friction
           if (vMag(pAngv[ip],3) > 1e-7) 
		     { for(int   k=0;k<3;k++)   
			   pAngv[ip][k] =  pAngv[ip][k] * (1.0 - min(1.0, pRMom[ip]/vMag(pAngv[ip],3)));
		     }
		   for(int   k=0;k<3;k++)   
            pAngp[ip][k] +=  pAngv[ip][k] * dem.dt;

            //---- check particles out of box ----
            if ( sqrt(pow(pPos[ip][0],2)+pow(pPos[ip][1],2)) >= box.Rxy ||
                pPos[ip][2] >= box.th ||  pPos[ip][2] <= box.bh) 
			 {
			  cout<<ip<< "PARITCLE OUT OF BOX"<<endl;
			  cout<<pPos[ip][0]<<" "<<pPos[ip][1]<<" "<<pPos[ip][2]<<endl;
			  system("PAUSE");
			  wrtTec();
		     }
            /*--- assembly region ----***/
           // dem.rmx = max(dem.rmx, vMag(pPos[ip],3));
            	compact.zmx = -999999;
				compact.zmn = 999999;
				for(int i=0; i<np ; i++)
				  {
					if (compact.zmx < ( pPos[i][2]+0.5*pDia[i]) )
					compact.zmx = pPos[i][2] + 0.5*pDia[i];
					if (compact.zmn > ( pPos[i][2]-0.5*pDia[i]) )
					compact.zmn = pPos[i][2] - 0.5*pDia[i];
				  }

            //---- need to rebuild nList? ----
		    for(int   k=0;k<3;k++)   
            pMov[ip][k] += pVel[ip][k] * dem.dt;
            if (vMag(pMov[ip],3) > qMovMx) 
				rBldNList = true;
		 }
        /*******---- determine*********/ 
       // ---- --- update box position 
		if (compact.stage != 0) 
           {
			   box.bh  +=  box.bhv[2] *dem.dt;
			   box.th  +=  box.thv[2] * dem.dt;
			   wMov[2] +=  abs(box.thv[2] * dem.dt);
				if (vMag(wMov,3) > qMovMx) 
					rBldNList = true;
				compact.dsp = compact.l_org - (compact.zmx - compact.zmn);
		    }
        //--- update loading mode ---
        if (compact.stage == 1 && compact.stress >= compact.stress_mx && compact.dsp > compact.dsp_mx) //(0=packing  1=compaction 2=relax 3=unloading 4=stop) 
			{	
				compact.stage = 2 ;
				dem.rlxtime = 1.1* dem.cTime  ; 
				box.thv[2] = box.thv[2] * 0.001 ;  
				cout<<"holding"<<endl;
			}	
        else if (compact.stage == 2 && dem.cTime > dem.rlxtime) 
			{	
				compact.stage = 3 ;
				box.thv[2] = -box.thv[2]*500.0 ;    // unloading is half of the loading speed
				cout<<"unloading"<<endl;
			}	
        else if (compact.stage == 3 && compact.stress < 1.e-8 )
            {
				compact.stage = 4 ;
				dem.simTime = 1.1 * dem.cTime  ;   // relax time (0.1*dem%simTime)
				cout<<"relaxing"<<endl;
			}	
       /***********---- analyse ----**********/
        if (dem.cTime >= dem.anatime) 
		  { 
			 analyz();
			// cout<<"Time= "<<dem.cTime<<"  Max overlap= "<<mxOlp<<endl;
			 dem.anatime += dem.anapert;
		  }
        
       /***---- dump file and restart if over time ----****/
        if (dem.cTime >= dem.dumptime) 
		{
			dump();
            dem.dumptime += dem.dumppert;
		}
        
        /*****---- data for tecplot ----******/
        if (dem.cTime >= dem.tptime) 
		{
			wrtTec();
            dem.tptime += dem.tppert;
		}
		//---- stop and restart ----
	/*	pclock();
        if (etime > dem.chkTime * 3600) 
		{ 
			dump();
			outfile(1);
			exit(1);
		}*/
        /*******---- everyting done,  move to next ----********/
          dem.cTime +=  dem.dt * dem.rtunit;

 }


   /** analyze the final result  ********/
     dump();
     analyz();
     wrtTec();
     deleteArray(pMov,npTotal);
    //write the jobfinished
	ofstream outfile;
    string filename="jobfinished";
	outfile.open(filename.c_str(),ios::out);
	outfile.close();

}