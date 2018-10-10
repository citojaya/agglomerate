#include <iostream>
#include <ctime>
#include <cmath>
#include "gVarExt.h"

   void feed() 
{
    const int nLoop = 10000000;
    double p,q,pdiaMax,feedRad,feedAngle,feedzmin,feedzmax;
    double ipPos[3], ijPos[3];
 /**************---- feed new particles ----**********/
    cout<<"feeding new particles."<<endl;
    
   //call random_seed

  //first particle as nuclei
   /* if (np == 0) 
	  {
          for (int j = 0; j< 3; j++)
		    {
             pPos[np][j] = 0.0;
             pVel[np][j] = 0.0;
             pAngp[np][j] = 0.0;
             pAngv[np][j] = 0.0;
		    }
	   }*/
    srand(static_cast<unsigned int>(time(NULL)));
    //--- setting up feeding space ---
	pdiaMax = Max(pDia, npTotal);
	 
	 feedRad = box.gRxy - pdiaMax/2.0;
	 feedzmin = box.gbh + pdiaMax/2.0;
     feedzmax = box.gth - pdiaMax/2.0;

	// random feed
    
    int iLoop = 0;
    LOOP:
	while ((np < npTotal) & (iLoop < nLoop))
	{  
        iLoop = iLoop + 1;
	    p=static_cast<double>(rand())/RAND_MAX;
		       
	    q=static_cast<double>(rand())/RAND_MAX;
        
		ipPos[0] = (p*feedRad) * cos(2*q*PI);
		ipPos[1] = (p*feedRad) * sin(2*q*PI);
        
	    p=static_cast<double>(rand())/RAND_MAX;
        ipPos[2] = p * (feedzmax - feedzmin) + feedzmin;  
        
        //check overlap
		
        for(int jp = 0;jp<np;jp++)
		   {
			 for(int i=0;i<3;i++)
			 {ijPos[i] = ipPos[i] - pPos[jp][i];}
           
            if ( vMag(ijPos,3) < (0.5*(pDia[np]+pDia[jp])+vGapMx)) goto LOOP; 
           }

     /******** new particle fixed********************/
	      for(int i=0;i<3;i++)
	    	{pPos[np][i] = ipPos[i];}
    /*********---- assign velocity******************/

        pVel[np][2] = -initVel/dem.rvunit;
     
		/*double sum=0;
	     for(int i=0;i<3;i++)
			{ sum+=pVel[np][i]*pPos[np][i];}
         
         if ( sum > 0.0)    //go outward, change to inward
	        { 
	         for(int i=0;i<3;i++)
			   pVel[np][i] = -pVel[np][i];
		    }*/
	     for(int i=0;i<3;i++)
		    {
			 pAngp[np][i] = 0.0;
             pAngv[np][i] = 0.0;
		    }
            np = np + 1;
    }
    cout<<" particles were generated. "<<np<<endl;

 }
