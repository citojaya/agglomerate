#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "gVarExt.h"
using namespace std;
/*------------------------------------
! for all kinds of force calculation
!------------------------------------*/
    static int nfc;;
    static double gap, fv,fcap,b_fc[3], b_mom[3];
    static double ipPos[3], ipVel[3], ipAngv[3];
    static double jpPos[3], jpVel[3], jpAngv[3];
    static double ijPos[3], unitVec[3];
    static double ipRad, ipDia, ipMass, ipHa, ipEmod, ipPois, ipDmpn, ipDmpt, ipSfrc, ipRfrc,ipYldp;
    static double jpRad, jpDia, jpMass, jpHa, jpEmod, jpPois, jpDmpn, jpDmpt, jpSfrc, jpRfrc,jpYldp;
    static bool   contact_with_top;

void force()

{
    void fcontact(int ip, int jp, int ine);
    void fvdw(int ine);
	void fbond(int ine);
	void fcp(int ine);
   
	double ijMag;

/***************Initializing the value********************/
	
    for(int   i=0;i<npTotal;i++)   
	   {  
		   pRMom[i]=0;
      //     pCprss[i]=0;
          for (int j= 0; j < 3; j++)
		   {
			 pForce[i][j]=0.0;
		     pMom[i][j]=0.0;
		   }
	    }
    for(int i=0;i<nne+1;i++)   
	    { 
          for (int j= 0; j < 10; j++)
          fList[i][j]=0.0;
	    }
    mxOlp = 0.0;
    nfc = 0;
    cn = 0;
	box.nc = 0;
	compact.stress =0.0;
   /* for(int   i=0;i<npTotal;i++)   
	{
       for (int j= 0; j < 3; j++)
	   {
		 pInf[i].upfc[j] = 0.0;
         pInf[i].upfv[j] = 0.0;
	   }
	}
    for(int   i=0;i<npTotal;i++)   
	{	 
	  pInf[i].compress = 0.0;
	  pInf[i].upcn = 0;
	  pInf[i].cn = 0;
	}*/
	/*for(int i=0;i<box.nWall;i++)  
       for(int j=0;j<2;j++)   
        for(int k=0;k<3;k++)   
		   {wForce[i][j][k] =0.0;}
	
	for(int i=0;i<3*npTotal;i++)  
      for(int j=0;j<10;j++)   
	   {wpfList[i][j]=0.0;}*/
/********--- set the mass centre ---***/
 /*   double sum_pPos[3], sum_pMass=0.0;
	  for(int j=0;j<3;j++) 
	  {
        sum_pPos[j]=0.0;
	  }

	for(int i=0;i<npTotal;i++)   
	{
	  sum_pMass +=pMass[i];
	  for(int j=0;j<3;j++)   
      sum_pPos[j] += pPos[i][j] * pMass[i];
	}
	  double sum_pMass_np=np*sum_pMass;
	  for(int k=0;k<3;k++)   
	  { oPos[k] = sum_pPos[k] /sum_pMass_np;
	  }*/


/******--- force calcuation from here ****************/
    tensile_stress = 0.0;
    tensile_cnv = 0;
    tensile_cnc = 0 ;   
    
    /********---- ip alias ----*******/
 for(int  ip=0;ip < np; ip++)
 {
	  for (int k= 0; k < 3; k++)
	   {
		 ipPos[k]   = pPos[ip][k];
         ipVel[k]   = pVel[ip][k];
         ipAngv[k]  = pAngv[ip][k];
	   }
        ipDia   = pDia[ip];
        ipRad   = 0.5 * pDia[ip];
        ipMass  = pMass[ip];

        int itype   = pType[ip];
        ipHa    = pMat[itype].ha;
        ipEmod  = pMat[itype].emod;
        ipPois  = pMat[itype].pois;
        ipSfrc  = pMat[itype].sfrc;
        ipRfrc  = pMat[itype].rfrc;
        ipDmpn  = pMat[itype].dmpn ;    
       // ipDmpt  = pMat[itype].dmpt;
		ipYldp  = pMat[itype].yldp;
		
     /******---- p-p force ----*******/

    for(int ine = nPoint[ip];ine< nPoint[ip+1];ine++)
	  {
          /***---- jp alias -****/
            int jp      = nList[ine];
	      for (int k= 0; k < 3; k++)
	        {
	    	 jpPos[k]   = pPos[jp][k];
             jpVel[k]   = pVel[jp][k];
             jpAngv[k]  = pAngv[jp][k];
	        }

			jpDia   = pDia[jp];
            jpRad   = 0.5 * pDia[jp];
            jpMass  = pMass[jp];
            int jtype   = pType[jp];
            jpHa    = pMat[jtype].ha;
            jpEmod  = pMat[jtype].emod;
            jpPois  = pMat[jtype].pois;
            jpSfrc  = pMat[jtype].sfrc;
            jpRfrc  = pMat[jtype].rfrc;
            jpDmpn  = pMat[jtype].dmpn ;    
//            jpDmpt  = pMat[jtype].dmpt;
			jpYldp  = pMat[itype].yldp;
          /*******---- ij distance ----*******/
	        for (int k= 0; k < 3; k++)
			{ijPos[k]   = ipPos[k] - jpPos[k];}
				 
            ijMag   = vMag(ijPos,3);
			for (int k= 0; k < 3; k++)
			{unitVec[k] = ijPos[k] / ijMag;}
            gap   = ijMag - (ipRad + jpRad);

         /********---- start from the force with longest range***********/
            if (b_flag == 1 &&rbond[ine] )
				{
					fbond(ine);
					for (int k= 0; k < 3; k++)
					 {
					  pForce[ip][k] += b_fc[k];
					  pForce[jp][k] -= b_fc[k];
					  pMom[ip][k] += b_mom[k] ;
					  pMom[jp][k] -= b_mom[k] ;
	 				 }
				}
			if (lqd.flag)
				{		
					lqd.cGapMx = (1+0.5*lqd.theta) * pow(lqd.lqV * 1.33333333333333 * PI * (pow(ipRad,3)+pow(jpRad,3)),0.33333333333333);
						lqd.layer =  (ipRad + jpRad) *(pow(1.0+lqd.lqV,0.33333333333333) - 1.0) ;
						if ( -gap < lqd.layer)
							lqdbrg[ine]= 1; // to form liquid bridge
						else if (-gap > lqd.cGapMx)
							lqdbrg[ine]= 0; // to break liquid bridge

						if (lqdbrg[ine])
							{
								fcp(ine);
								 for (int k= 0; k < 3; k++)
								 { 
									 pForce[ip][k] = pForce[ip][k] + fcap * unitVec[k];
									 pForce[jp][k] = pForce[jp][k] - fcap * unitVec[k];
								 }

							}
				}
			if (gap < vGapMx) 
			  { 
				 fvdw(ine);   //p-p vdw force
	             for (int k= 0; k < 3; k++)
				 { pForce[ip][k] = pForce[ip][k] + fv * unitVec[k];
                   pForce[jp][k] = pForce[jp][k] - fv * unitVec[k];
				 }
               //tensile strength
			       tensile_stress += ipRad * abs(fv) + jpRad * abs(fv);
                   tensile_cnv ++;
	           } 

			if (gap < 0.0) 
		        {
				 fcontact(ip, jp, ine);
                // pInf[ip].cn++;
			    // pInf[jp].cn++;               
               
			    }
			else
	            {
					ch_nDsp_pp[ine]= 0.0 ;
					ch_nDsp_mx_pp[ine]= 0.0 ;
					ch_plstDfm_pp[ine]= 0.0 ;
					ch_plstRad_pp[ine]= 0.0;
					for (int k= 0; k < 3; k++)
						 dispt_pp[ine][k] = 0.0;
				}
	    }//end p-p force
     /***********---- p-w force ----********/
	  for (int iWall = 0; iWall < box.nWall; iWall++)   //0: Z, 1:XY,
	   {	 
            if (contact_with_wall[ip][iWall] == 0) continue;    
            if (iWall == 0)                            // z direction
			  { if (contact_with_wall[ip][iWall] == 1)     // top
			      {
				    jpPos[0] = ipPos[0];
			        jpPos[1] = ipPos[1];
				    jpPos[2] = box.th;
					for (int k= 0; k < 3; k++)
						jpVel[k]   = box.thv[k];
			      }
                else if (contact_with_wall[ip][iWall] == -1)   // bottom
			      {
				    jpPos[0] = ipPos[0];
			        jpPos[1] = ipPos[1];
				    jpPos[2] = box.bh;
					for (int k= 0; k < 3; k++)
						jpVel[k]   = box.bhv[k];
			      }
			   }                                  
            else if (iWall == 1)                   // xy direction
			   {
				  jpPos[0] = box.Rxy * ipPos[0]/(sqrt(pow(ipPos[1],2)+pow(ipPos[0],2)));
				  jpPos[1] = box.Rxy * ipPos[1]/(sqrt(pow(ipPos[1],2)+pow(ipPos[0],2)));
				  jpPos[2] = ipPos[2];
					 for (int k= 0; k < 3; k++)
					 jpVel[k]   = 0.0;
			   }
            
          /*******---- wall material property ----*********/
            jpRad   = 0.0 ;     //zero
            jpDia   = 0.0 ;
            jpMass  = 0.0 ;

	        for (int k= 0; k < 3; k++)
			{
			 jpAngv[k]  = 0.0;
			 ijPos[k]   = ipPos[k] - jpPos[k];
			}
            jpHa    = pMat[0].ha   ;
            jpEmod  = pMat[0].emod;
            jpPois  = pMat[0].pois;
            jpSfrc  = pMat[0].sfrc;
            jpRfrc  = pMat[0].rfrc;
            jpDmpn  = pMat[0].dmpn ;    
          //  jpDmpt  = pMat[0].dmpt;
			jpYldp  = pMat[0].yldp;
                                
            ijMag   = vMag(ijPos,3);
	        for (int k= 0; k < 3; k++)
			{unitVec[k] = ijPos[k] / ijMag;}
            gap     = ijMag - (ipRad+jpRad);
            
           if (gap < vGapMx) 
			 {  
				 fvdw(-1);
	            for (int k= 0; k < 3; k++)
			    pForce[ip][k] += fv * unitVec[k];
			  }
            if (gap < 0.0) 
			   fcontact(ip, iWall, -1);//iwall
            else
				{
					ch_nDsp_pw[ip][iWall]= 0.0 ;
					ch_nDsp_mx_pw[ip][iWall]= 0.0 ;
					ch_plstDfm_pw[ip][iWall]= 0.0 ;
					ch_plstRad_pw[ip][iWall]= 0.0;
					for (int k= 0; k < 3; k++)
						dispt_pw[ip][iWall][k] = 0.0;
				}
	    } 
        //---- gravity ----
        pForce[ip][2] = pForce[ip][2] - dem.ctrF*pMass[ip];
	    
		//---- seperate force ----
			//if ( ipPos[0] <= 0.0 )
	  	//		pForce[ip][1] = pForce[ip][1] - dem.ctrF *pMass[ip];
		//	else
	  	//		pForce[ip][1] = pForce[ip][1] + dem.ctrF *pMass[ip];
			
	   /* double ipMag = vMag(ipPos,3);
        if (ipMag >1e-10) 
		  {
			pForce[ip][0] = pForce[ip][0] - dem.ctrF * pMass[ip] * ipPos[0] / ipMag;
            pForce[ip][1] = pForce[ip][1] - dem.ctrF * pMass[ip] * ipPos[1] / ipMag;
            pForce[ip][2] = pForce[ip][2] - dem.ctrF * pMass[ip] * ipPos[2] / ipMag;
		  }*/
	  				/***write force **********/
			    //if (dem.cTime >=dem.anatime && ip == 0) 
			  //  if (dem.cTime <= 400*dem.dt * dem.rtunit && ip == 0) 
			   /* if ( ip == 0) 
					{ 
					   ofstream outfile;
					   string filename=genfile+"_force.txt";
					   outfile.open(filename.c_str(),ios::out|ios::app);
    				   outfile<<setprecision(10)<<setiosflags(ios::scientific);
    				   outfile<<dem.cTime<<" "<< pVel[ip][0]*dem.rvunit <<" "<<gap/ipDia <<" "<<fcap<<endl;
					   outfile.close();
					   dem.anatime += dem.anapert;
					   if (gap>0.0)
						   exit(1);
					} */

   }//end force calcuation  
}// end force fouction


void fcontact(int ip, int jp, int ine)
{
    void vprod(double *a, double *b, double*c);   // cross product C = A x B 
    void vproject(double *v, double *n, double *t, int itype);
    void vproject(vector<double> v, double *n, double *t, int itype);
	double ijDsmx, static_force_ratio, dispt_ratio;
    double nrmDsp;
    double pkt, ct;
    double ijRad, ijEmod, ijSfrc, ijPois, ijRfrc, ijDmpn, ijDmpt, ijMass;
    double rsum;
    double ipCPVel[3], jpCPVel[3], ijCPVel[3];
    double ijCPVn, ijCPVt[3];
    double disptTotal[3];
    double ipRVec[3], jpRVec[3], rot[3];
    double fn, fcn, fdn;
    double fc[3], ft[3], fct[3], fdt[3];
    double dMom[3];
	int ch_stage; // 0='nc',1='uld', 2='ld',3='rld'

	/*****--- reduced unit*****/
    nrmDsp  = -gap;
	ijEmod  = (ipEmod * jpEmod)/(ipEmod + jpEmod);
    ijSfrc = pow(ipSfrc * jpSfrc , 0.5);
    ijPois = pow(ipPois * jpPois , 0.5);
    ijDmpn = (ipDmpn * jpDmpn) / (ipDmpn + jpDmpn);
    
    if (ine == -1) 
	   { 
		 ijRad = ipRad;
         ijMass = ipMass;
	   }
    else 
	   {
		 ijRad = ipRad * jpRad / (ipRad + jpRad);
         ijMass = ipMass * jpMass / (ipMass + jpMass);
	   }
		rsum = ipRad + jpRad;
	
	/*****--- yielding point parameters*****/
		yld.stress  = 2.0* ipYldp * jpYldp / (ipYldp + jpYldp) ;
		yld.dsp     = pow(PI*yld.stress/(2.0*ijEmod),2.0) * ijRad;
		yld.force   = 4.0/3.0 * ijEmod * sqrt(ijRad * yld.dsp) * yld.dsp ;
		yld.stiffK  = 2.0*ijEmod * sqrt(ijRad*yld.dsp);

	//============================================
	//	load, unload and reload 0='no force although in contact 1='unloading', 2='loading',3='reloading'
	//============================================
		if (ine == -1) //wall
		   { 
				if (nrmDsp <= ch_plstDfm_pw[ip][jp])         // no force although in contact
					ch_stage = 0;
				else if (nrmDsp < ch_nDsp_pw[ip][jp])				  // unloading stage
					ch_stage = 1;
				else if (nrmDsp > ch_nDsp_mx_pw[ip][jp])	    // loading stage
					ch_stage = 2;
				else 
					ch_stage = 3;						// reloading stage
			/*****---- normal force ----*********/
				if (ch_stage == 0)          
					fcn = 0.0;
				else if (ch_stage == 2)       
				{
					if (nrmDsp < yld.dsp)     
					//-- hertzian force ---
						{
						fcn = 1.333333333333 * ijEmod * sqrt(ijRad * nrmDsp) * nrmDsp; 
						ch_plstRad_pw[ip][jp] = ijRad;  
						ch_plstDfm_pw[ip][jp] = 0.0;    
						}
					else    
						//---plastic force ---
						{
						fcn = yld.force + yld.stiffK * (nrmDsp - yld.dsp);
						//--- harden effect ---
						ch_plstRad_pw[ip][jp] = ijRad * (4.0/3.0 *ijEmod*sqrt(ijRad * nrmDsp)*nrmDsp)/fcn;
						ch_plstDfm_pw[ip][jp] = nrmDsp - pow(fcn/(1.333333333333 * ijEmod *sqrt(ch_plstRad_pw[ip][jp])),(2.0/3.0));
						}
					ch_nDsp_mx_pw[ip][jp] = nrmDsp;
				}
	    
				else 
				//--- unload and reload follows a same cuve
					{
						fcn =  1.333333333333 * ijEmod * sqrt(ch_plstRad_pw[ip][jp]) * pow(nrmDsp - ch_plstDfm_pw[ip][jp],1.5);
					}
				ch_nDsp_pw[ip][jp] = nrmDsp;

		   }
		else 
		   {

			   if (nrmDsp <= ch_plstDfm_pp[ine])         // no force although in contact
					ch_stage = 0;
				else if (nrmDsp < ch_nDsp_pp[ine])				  // unloading stage
					ch_stage = 1;
				else if (nrmDsp > ch_nDsp_mx_pp[ine])	    // loading stage
					ch_stage = 2;
				else 
					ch_stage = 3;						// reloading stage
			 /*****---- normal force ----*********/
				if (ch_stage == 0)          
					fcn = 0.0;
				else if (ch_stage == 2)       
				{
					if (nrmDsp < yld.dsp)     
					//-- hertzian force ---
						{
						fcn = 1.333333333333 * ijEmod * sqrt(ijRad * nrmDsp) * nrmDsp; 
						ch_plstRad_pp[ine] = ijRad;  
						ch_plstDfm_pp[ine] = 0.0;    
						}
					else    
						//---plastic force ---
						{
						fcn = yld.force + yld.stiffK * (nrmDsp - yld.dsp);
						//--- harden effect ---
						ch_plstRad_pp[ine] = ijRad * (1.333333333333 *ijEmod*sqrt(ijRad * nrmDsp)*nrmDsp)/fcn;
						ch_plstDfm_pp[ine] = nrmDsp - pow(fcn/(1.333333333333 * ijEmod *sqrt(ch_plstRad_pp[ine])),(2.0/3.0));
						//ch_plstRad_pp[ine] = 4.0/3.0 * ijEmod * pow((2.0 * fcn+ yld.force)/(2.0*PI*yld.stress),1.5)/fcn; 
						//ch_plstDfm_pp[ine] = nrmDsp - pow(fcn/(1.333333333333 * ijEmod *sqrt(ch_plstRad_pp[ine])),(2.0/3.0));
						}
					ch_nDsp_mx_pp[ine] = nrmDsp;
				}
	    
				else 
				//--- unload and reload follows a same cuve
					{
						fcn =  1.333333333333 * ijEmod * sqrt(ch_plstRad_pp[ine]) * pow(nrmDsp - ch_plstDfm_pp[ine],1.5);
					}
				ch_nDsp_pp[ine] = nrmDsp;
		   }
                
   /***----------- normal damping force --------------***/            
	 for (int k= 0; k < 3; k++)
	 {ipRVec[k] = -ipRad * (1 - nrmDsp /rsum) * unitVec[k];}
      vprod(ipAngv, ipRVec, rot);
	 for (int k= 0; k < 3; k++)
	  {
	   ipCPVel[k] = ipVel[k] + rot[k];
	   jpRVec[k] = jpRad * (1 - nrmDsp / rsum) * unitVec[k];
	  }
       vprod(jpAngv, jpRVec, rot);
	 
	 ijCPVn=0.0;
	 for (int k= 0; k < 3; k++)
	   {
		jpCPVel[k] = jpVel[k] + rot[k];
        ijCPVel[k] = ipCPVel[k] - jpCPVel[k];
        ijCPVn += ijCPVel[k] * unitVec[k];
	    }
    fdn    = -ijDmpn* sqrt(ijRad * nrmDsp) * ijCPVn; 
    /*------------- normal force ------------------***/
    fn = fcn + fdn;
    /**------------ tangential force ------------------***/
    
	   vproject(ijCPVel, unitVec, ijCPVt, 0);

    if (ine == -1)      //p-w force
       vproject(dispt_pw[ip][jp], unitVec, disptTotal, 1);
    else 
       vproject(dispt_pp[ine], unitVec, disptTotal, 1);
	  for (int k= 0; k < 3; k++)
	  {disptTotal[k] += ijCPVt[k] * dem.dt;}
    
    /*---- ratio with maximum displacement----******/

    ijDsmx  = ijSfrc * (2 - ijPois) / (2 - 2 * ijPois);
    
    dispt_ratio = vMag(disptTotal,3)/(ijDsmx * nrmDsp);
	   for (int k= 0; k < 3; k++)
       disptTotal[k] = disptTotal[k] / max(1.0, dispt_ratio);

    static_force_ratio = 1 - pow(1 - min(1.0, dispt_ratio),1.5);

    if (vMag(disptTotal,3) > 0.0) 
	   for (int k= 0; k < 3; k++)
		fct[k] = -ijSfrc * fn * static_force_ratio * disptTotal[k]/vMag(disptTotal,3);
    else
	   for (int k= 0; k < 3; k++)
        fct[k] = 0.0;
    
   /* --------- tangential damping force-------
     if (dispt_ratio > 1.0) 
	    for (int k= 0; k < 3; k++)
        fdt[k] = 0.0;
     else
	 {
		 pkt = 1.5 * ijSfrc * fn * pow( 1.0 - dispt_ratio ,0.5)/ (ijDsmx * nrmDsp);

        if (ine == -1) 
            ijMass = ipMass;
        else
            ijMass = ipMass * jpMass / (ipMass + jpMass);

        ijDmpt = ipDmpt * jpDmpt / (ipDmpt + jpDmpt);
		ct = 2.0 * ijDmpt * pow(ijMass * pkt , 0.5);
	    for (int k= 0; k < 3; k++)
		fdt[k] = -ct * ijCPVt[k];
	 }*/

    if (ine == -1) 
	   for (int k= 0; k < 3; k++)
        dispt_pw[ip][jp][k] = disptTotal[k];
    else
	   for (int k= 0; k < 3; k++)
        dispt_pp[ine][k] = disptTotal[k];
   /****------------- sum forces and moments-----********/
	 for (int k= 0; k < 3; k++)
	  {
        //ft[k] = fct[k] + fdt[k];
        ft[k] = fct[k];
        fc[k] = fn * unitVec[k] + ft[k];
	  }
    /****--- ip force and torque ---***********/
	for (int k= 0; k < 3; k++)
	 { pForce[ip][k] += fc[k];}
    vprod(ipRVec, fc, dMom);
	for (int k= 0; k < 3; k++)
	 {pMom[ip][k] = pMom[ip][k] + dMom[k];}
    ijRfrc = 0.5 * (ipRfrc + jpRfrc);
    pRMom[ip] += ijRfrc * fn;
  /*  pCprss[ip] += fcn;
 	pInf[ip].compress += fcn;
    
	if(ipPos[2] < jpPos[2]) 
	 {	
		 for (int k= 0; k < 3; k++)
		 {
			pInf[ip].upfc[k] +=  fc[k];
            pInf[ip].upfv[k] += fv * unitVec[k];
		 }

        pInf[ip].upcn ++;
	 }*/
   /********--- jp force and torque ---************/
    //int iw;
	if (ine == -1) 
	{  
       //*****  record the Wall force  ******
	/*	double sum=0.0;
	   for (int k= 0; k < 3; k++)
	    {
		  sum += unitVec[k];              // -1(top,right) or 1(bottom, left)
	    }
	      iw=int(sum);
          iw = (abs(iw) + iw)/2;   // 0 for iw =-1(top),  1 for iw = 1
	      
		  for (int k= 0; k < 3; k++)
		  {wForce[jp][iw][k] += - fc[k] ; }   
          
	      for (int k= 0; k < 3; k++)
		  {
		   wpfList[box.nc][k] = jpVel[k];
           wpfList[box.nc][k+3] = fn * unitVec[k];
           wpfList[box.nc][k+6] = ft[k];
		  }
           wpfList[box.nc][9]  = fv;
		   box.nc++ ; */   
		if (jp == 0 && contact_with_wall[ip][jp] == 1) //top wall
			compact.stress +=  fn ;
	}
    else 
	{
	    for (int k= 0; k < 3; k++)
		{
		  pForce[jp][k] = pForce[jp][k] - fc[k];
		}
		double _fc[3];   //_fc is -fc
		for (int k= 0; k < 3; k++)
		{		
          _fc[k]=-fc[k];
		}
		vprod(jpRVec, _fc, dMom);
		for (int k= 0; k < 3; k++)
        pMom[jp][k] +=  dMom[k];
        pRMom[jp] += ijRfrc * fn;
      //  pCprss[jp] += fcn ;
       
		/*pInf[jp].compress +=  fcn;
        if(ipPos[2] > jpPos[2]) 
		 {
		  for (int k= 0; k < 3; k++)
		   {  
			pInf[jp].upfc[k] +=  - fc[k];
            pInf[jp].upfv[k] +=  - fv * unitVec[k];
		   } 
			pInf[jp].upcn ++;
		 }*/
	
//**********   record the contact force in the array fList  **************
	/*	for (int k= 0; k < 3; k++)
	      fList[nfc][k] = fn * unitVec[k];

        if(ipPos[2] > jpPos[2]) 
		{
			for (int k= 3; k < 6; k++)
            fList[nfc][k] = ft[k-3];
		}
        else
		{
		  for (int k= 3; k < 6; k++)
			fList[nfc][k] = -ft[k-3];
        }
            fList[nfc][6] = fv;
            fList[nfc][7] = ip;
            fList[nfc][8] = jp;
            fList[nfc][9] = nfc;
			nfc++; */
	}
//*************    check the max overlap  *************		
    mxOlp = max(mxOlp, nrmDsp/(2*ipDia));
    if (mxOlp > 0.4) 
      {
	  cout<<mxOlp<< "overlap excedes 40%"<<endl;
      cout<<ip<<" "<<jp<<" "<<nrmDsp<<endl;
	  system("PAUSE");
	  mxOlp = 0.0;
      }
}

void fvdw(int ine)

{
    double ijRad, vGap, ijHa;
    double lamda, b, retarded_factor;
    
    vGap = max(gap, vGapMn);

    /* beware that Ha has already been divided by 6.0*/
    ijHa = pow(ipHa * jpHa,0.5) ;
    
    /*----------------------------------------------------
    ! pp = ha/6.0 * di^3 * dj^3 * (h +0.5di + 0.5 dj) /
    ! ((h^2 + dih + djh)^2 * (h^2+dih+djh+didj0^2)
    !
    ! pw = ha/6.0 * di^3 * 0.5 / (h * (h+di))^2
    !
    ! h << d: 
    ! pp = ha/24 * d / h^2
    ! pw = ha/12 * d / h^2
    !----------------------------------------------------
    if (ine == -1) 
        fv = -ijHa * pow(ipDia,3) * 0.5 / pow( vGap * (vGap+ipDia),2 );
    else
        fv = -ijHa * pow(ipDia * jpDia,3) * (vGap + 0.5 * (ipDia + jpDia)) 
        /pow( ( (pow(vGap,2) + vGap * ipDia + vGap * jpDia) 
        * (pow(vGap,2) + ipDia*vGap + jpDia*vGap + ipDia*jpDia) ),2);

    /*------------------------------------------------
    !consider retareded effect 
    !(J.Gregory, J. Col. Int. Sci., 83 1981)
    ! V=ha/6.0 * (r*) /h *(1-b*h/lamda*ln(1+lamda/b/h))
    ! f=ha/6.0 * (r*) /(h^2) * (1-1/(1+lamda/(b*h)))
    ! lamda = 100nm; b = 5.32
    !---------------------------------------------------*/
    if (ine == -1) 
        ijRad = 0.5 * ipDia; 
    else
        ijRad = 0.5 * (ipDia * jpDia)/(ipDia+jpDia);
    lamda = 1.0e-7;      //100nm
    b = 5.32;
    retarded_factor = 1.0 - 1.0 /(1.0 + lamda/(b*(vGap*dem.rlunit)));
    fv = -ijHa * ijRad /(vGap*vGap) * retarded_factor;
	}

void fbond(int ine)
{
	void vproject(double *v, double *n, double *t, int itype);
    void vproject(vector<double> v, double *n, double *t, int itype);
    void vprod(double *a, double *b, double*c);   // cross product C = A x B 
	double bond_rad,bond_A,bond_I,bond_J,bond_stiff_n,bond_stiff_s; // area, moment inertia and polar moment of inertia
    double nrmDsp, b_fcn,b_fct[3],b_mom_n, b_mom_t[3], b_tensile, b_shear, disptTotal[3],b_fdn,b_fn;
	double ipCPVel[3], jpCPVel[3], ijCPVel[3];
    double ipRVec[3], jpRVec[3], rot[3], ijCPVt[3],ijAngv[3],ijAngv_t[3], rsum,ijAdip_n,ijAdip_t[3],ijCPVn,ijDmpn;
		
		ijAdip_n = 0.0;
		rsum = ipRad + jpRad;
		bond_rad = rad_multi* min (ipRad, jpRad); // R = rad_multi * min (Rad_i, Rad_j)
		bond_A = PI * pow(bond_rad,2); // A= pi*r^2 
		bond_I = PI/4.0 * pow(bond_rad,4);// I= 1/4* pi* r^4
		bond_J = 2.0*bond_I ;// I= 1/2* pi* r^4 = I/2

		bond_stiff_n = b_emod_n /(rsum);   // kn = Emod / (Ri + Rj)
		bond_stiff_s = bond_stiff_n / ratio_stiff; //  ks = kn / ratio
		ijDmpn = (ipDmpn * jpDmpn) / (ipDmpn + jpDmpn);

        /*****----normal force from bond ----*********/
		 nrmDsp  = -gap - initalGap_bond[ine];
		 b_fcn   = bond_stiff_n * bond_A * nrmDsp ; 
        /*****--- tangential force from bond ----*********/
		for (int k= 0; k < 3; k++)
		 ipRVec[k] = -ipRad * (1 - nrmDsp /rsum) * unitVec[k];
		 vprod(ipAngv, ipRVec, rot);
		for (int k= 0; k < 3; k++)
		 {
			ipCPVel[k] = ipVel[k] + rot[k];
			jpRVec[k] = jpRad * (1 - nrmDsp / rsum) * unitVec[k];
		 }
         vprod(jpAngv, jpRVec, rot);
	 
		for (int k= 0; k < 3; k++)
		{
			jpCPVel[k] = jpVel[k] + rot[k];
			ijCPVel[k] = ipCPVel[k] - jpCPVel[k];
   	    }
		  vproject(dispt_bond[ine], unitVec, disptTotal, 1);
          vproject(ijCPVel, unitVec, ijCPVt, 0);
		for (int k= 0; k < 3; k++)
	      disptTotal[k] += ijCPVt[k] * dem.dt;
		for (int k= 0; k < 3; k++)
			b_fct[k] =-bond_stiff_s * bond_A *  disptTotal[k];
		for (int k= 0; k < 3; k++)
		  dispt_bond[ine][k] = disptTotal[k];
      /*****--- dumping force from bond ----*********/
			 ijCPVn=0.0;
			 for (int k= 0; k < 3; k++)
				ijCPVn += ijCPVel[k] * unitVec[k];
			b_fdn = -ijDmpn* sqrt( abs(nrmDsp) * bond_rad ) * ijCPVn; 
    /*------------- normal force ------------------***/
		b_fn = b_fcn + b_fdn;
      /*****--- sum force from bond ----*********/
		for (int k= 0; k < 3; k++)
		b_fc[k] = b_fn * unitVec[k] + b_fct[k];
      /*****----normal MOM from bond ----*********/
		for (int k= 0; k < 3; k++)
			ijAngv[k] = ipAngv[k]-jpAngv[k];
		
		for (int k= 0; k < 3; k++)
			ijAdip_n += ijAngv[k] * unitVec[k]* dem.dt; // theta_n = (wi-wj)*dt
			
		 dispAn_bond[ine] += ijAdip_n; // accumulate normall angular displacement 

		b_mom_n= -bond_stiff_s * bond_J * dispAn_bond[ine] ;
      /*****----tangential MOM from bond ----*********/
          vproject(ijAngv, unitVec, ijAngv_t, 0);
		for (int k= 0; k < 3; k++)
		  ijAdip_t [k] = ijAngv_t[k] * dem.dt; 
		for (int k= 0; k < 3; k++)
		  dispAt_bond[ine][k] += ijAdip_t[k]; // accumulate tangential angular displacement 
		for (int k= 0; k < 3; k++)
		  b_mom_t[k] = -bond_stiff_n * bond_I * dispAt_bond[ine][k];
      /*****--- sum MOM from bond ----*********/
		for (int k= 0; k < 3; k++)
			b_mom[k] = b_mom_n * unitVec[k] + b_mom_t[k];

      /*****--- judgement of bond break based on the maxmimum of tensile and shear stress  ----*********/
          
		b_tensile = (-b_fcn/bond_A) + (vMag(b_mom_t,3)*bond_rad/bond_I);
		b_shear = (vMag(b_fct,3)/bond_A) + (abs(b_mom_n)*bond_rad/bond_J);
		
        if ((b_tensile >= max_b_tensile) || (b_shear >= max_b_shear))
		{
		  rbond[ine] = 0 ; 
		}

}

void fcp(int ine)
{
	double ijDia,cGap,Vm; // effective radius, critical rupture separation, modified Liquid bridge volume

     cGap = max (gap,lqd.cGapMn);

	if (ine == -1) 
        {
			ijDia = 2* ipDia; 
			Vm = 0.5 * lqd.lqV /6.0 * PI * pow(ipDia,3); // liquid volume for p-w
		}
	else
        {
			ijDia =  2.0 * (ipDia * jpDia)/(ipDia+jpDia) ; 
			Vm = lqd.lqV /6.0 * PI * (pow(ipDia,3) + pow(jpDia,3)); // liquid volume for p-w
		}
			fcap= -PI*ijDia*lqd.sTension*cos(lqd.theta)/(1.0+1.0/(-1.0 +sqrt(1.0 + 4.0 * Vm/(PI*ijDia*pow(cGap,2))))); 



}