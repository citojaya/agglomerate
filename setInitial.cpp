#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cmath>
#include "gVarExt.h"
using namespace std;

void setInitial()
{
//declare function 
void setArray(int flag);
void setParticle();
void outfile(int);
void readDump(ifstream &infile);
void getvGapMx(double a, double hmn, double &hmx);
void setPV();


//define variables
 const double g = 9.81;
 int ip, itype;
 double massMn, sizeMn, stiffk,yld_vel,v_ratio,emod_star;
 string filename;

 /**** looking for file for restart ****/
   filename = genfile+".dump";
   ifstream  infile; 
   infile.open(filename.c_str(),ios::in); 
   if(!infile) 
     { 
        np  = 0;
        nne = 0; 
        dem.cTime = 0.0;

 	    setArray(1);
        setArray(2);
        setParticle();
     } 
     else 
     { 
     cout<<"found for restart for"<<" "<<filename<<endl; 

	    readDump(infile);
     } 
   infile.close();
  /******** particle mass*********************/ 

  for(ip=0;ip<npTotal;ip++)
     {
      itype = pType[ip];
      pMass[ip] = 1.0/6.0 * PI * pow(pSize*pDia[ip],3)* pMat[itype].density;

     }
  //---- time step ----*********
    massMn  = Min(pMass,npTotal);
    sizeMn  = Min(pDia,npTotal) * pSize;
	//langston
    //stiffk  = pMat[1].emod * pow(pSize , 0.5);
    //dem.dt  = dem.dtFactor * PI * pow(massMn / stiffk,0.5);  
    stiffk= 0.666666666666*pMat[1].emod*pow(sizeMn/4,0.5);
    dem.dt  =  dem.dtFactor*3.21*pow(massMn/stiffk ,0.4)*pow(10,-0.2); //  schafer
   //----------------------------------------------------------- 
  //estimate the restitution coefficient from the yield stress 
  //for sphere with plate at 1m/s  
  //--------------------------------------------------------------
	    //--- p-w resitution coefficient
    emod_star  = (pMat[0].emod * pMat[1].emod)/(pMat[0].emod + pMat[1].emod);
	yld_vel = 1.56 * pow(pow(pMat[0].yldp,5.0)/(pow(emod_star,4.0)*pMat[0].density),0.5);
    v_ratio = min(1.0, yld_vel);
    yld.rstCefpw = pow(1.2*pow(3.0,0.5),0.5)*pow(1-pow(v_ratio,2.0)/6.0,0.5)*pow(v_ratio/(v_ratio+2.0*pow(1.2-0.5*pow(v_ratio,2.0),0.5)),0.25);  
    //--- p-p restitution coefficient
    yld_vel = 3.194*pow(pow(pMat[1].yldp,5.0)*pow(sizeMn/4.0,3.0),0.5)/pow(pow(pMat[1].emod/2.0,4.0)*massMn/2.0,0.5); //
	v_ratio = min(1.0, yld_vel);
    yld.rstCefpp = pow(1.2*pow(3.0,0.5),0.5)*pow(1-pow(v_ratio,2.0)/6.0,0.5)*pow(v_ratio/(v_ratio+2.0*pow(1.2-0.5*pow(v_ratio,2.0),0.5)),0.25);
    //---- convert to the reduced unit ----
    dem.rlunit  = pSize;
    dem.rmunit  = 1.0/6.0 * PI * pMat[1].density * pow(pSize,3);
    dem.rfunit  = dem.rmunit * g;
    dem.rsunit  = dem.rfunit / pow(pSize , 2);
    dem.rtunit  = sqrt(pSize / g);
    dem.rvunit  = dem.rlunit / dem.rtunit;       //= sqrt(pSize * g)
    dem.dt      = dem.dt / dem.rtunit;
    
     for(int i=0; i<npTotal;i++)
	 {
	  pMass[i]       = pMass[i] /dem.rmunit;
      pInert[i]      = 0.1 * pMass[i] * pow(pDia[i],2);
	 }
     for(int i_ph=0; i_ph<=ntype;i_ph++)
	  {
	   pMat[i_ph].ha=	pMat[i_ph].ha/dem.rlunit / dem.rfunit / 6.0;
       pMat[i_ph].dmpn   = pMat[i_ph].dmpn * pMat[i_ph].emod * dem.rlunit * dem.rvunit / dem.rfunit;
	   pMat[i_ph].emod   = pMat[i_ph].emod  / dem.rsunit;
	   pMat[i_ph].yldp   = pMat[i_ph].yldp/ dem.rsunit;
      }
      //bond to reduced unit
      bGapmx = max_b_tensile/b_emod_n* Max(pDia,npTotal) ;
	  b_emod_n = b_emod_n / dem.rsunit;
	  max_b_tensile = max_b_tensile / dem.rsunit;
	  max_b_shear =	max_b_shear  / dem.rsunit;

	  // caplliary force reduced unit
	  lqd.cGapMn = lqd.cGapMn /	dem.rlunit;
	  lqd.sTension = lqd.sTension * dem.rlunit/dem.rfunit;
	  lqd.brknmx = (1+0.5*lqd.theta) * pow(lqd.lqV * PI / 3.0 * pow(Max(pDia,npTotal),3),0.33333333333333);
	  //vdw force reduced unit
      vGapMn = vGapMn / dem.rlunit;
      vGapMx = vGapMx / dem.rlunit;
	  getvGapMx(pMat[1].ha, vGapMn, vGapMx);
	  //compaction pressure
	  compact.stress_mx = compact.stress_mx / dem.rsunit * (PI*pow(box.Rxy,2));
    // set box velocity
	if (compact.stage == 0 )
	 {
		 for(int j=0;j<3;j++)
			{
				box.thv[j] = 0.0;
				box.bhv[j] = 0.0;
			}
	 }
	  // set the cutgap 
	  if (b_flag == 1 && lqd.flag ==1)
			{
			 cutGap = max( cutGap, bGapmx);
			 cutGap = max ( cutGap, lqd.brknmx);
			}	
	  else if (b_flag == 1)
		  {
			  cutGap = max( cutGap, bGapmx);
		  }
	      
	  else if (lqd.flag)
		 {
			 cutGap = max ( cutGap, lqd.brknmx);
		 }
    /*****---- intialise dem other time parameters ----****/
    dem.fdtime  = dem.cTime;
    dem.dumptime= dem.cTime;
    dem.anatime = dem.cTime;
    dem.tptime  = dem.cTime;

    /***---- output some essential information ----*/
    outfile(0);
}

//---- function of initialize global array ----
void setArray(int flag)
{
	if(flag==1)
	{   
        makeArray(ch_nDsp_pw,npTotal,box.nWall); 
        makeArray(ch_nDsp_mx_pw,npTotal,box.nWall);  
        makeArray(ch_plstDfm_pw,npTotal,box.nWall);   
        makeArray(ch_plstRad_pw,npTotal,box.nWall);   
		makeArray(pPos,npTotal,3);
        makeArray(pVel,npTotal,3);   
        makeArray(pAngp,npTotal,3);    
        makeArray(pAngv,npTotal,3);   
        makeArray(pForce,npTotal,3);   
		makeArray(pMom,npTotal,3); 
		makeArray(pType,npTotal);     
        makeArray(pDia,npTotal);   
        makeArray(pMass,npTotal);    
        makeArray(pInert,npTotal);   
        makeArray(pRMom,npTotal);       
        makeArray(pCprss,npTotal);   
        makeArray(nPoint,npTotal+1);
        makeArray(dispt_pw,npTotal,box.nWall,3); 
        makeArray(contact_with_wall,npTotal,box.nWall);    
        makeArray(wForce,box.nWall,2,3);  
        makeArray(wpfList,3*npTotal,10);  
		pInf=new   pInfo[npTotal]; 
	}
	else if(flag==2)
	{
	    //---- arrays relating to neighbor number ----
        nList.resize(nne);         
        ReSize(dispt_pp,3, nne);    
        ReSize(fList,fl_size, nne+1);         
		rbond.resize(nne);
		lqdbrg.resize(nne);
		dispAn_bond.resize(nne);
		initalGap_bond.resize(nne);
        ReSize(dispAt_bond,3, nne);
        ReSize(dispt_bond,3, nne);
	    ch_nDsp_pp.resize(nne);
		ch_nDsp_mx_pp.resize(nne);
		ch_plstDfm_pp.resize(nne);
		ch_plstRad_pp.resize(nne);
	}
	else
    cout<<"allocating array erro"<<endl;
    

}
//---- set all information for particles except position ----

void setParticle()
{

    int ip, itype;
   // double p, rnp, ratio;
   // bool fixed;
    
    srand(static_cast<unsigned int>(time(NULL)));  // random feed

   // vunit = sqrt(pSize * 9.81);
   // initVel = initVel / vunit;
      
        // wall properties
	ip = 0;
			
	for (itype = 1; itype <= ntype; itype++)
		{
		 for(int i = ip;i < ip + pMat[itype].np;i++)
			{
				pType[i] = itype;
                pDia[i] = pMat[itype].dia;
			}
		 ip += pMat[itype].np;
		}


    cout<<"setMat done."<<endl;
}

void readDump(ifstream &infile)
{   
	void fndrec(ifstream &infile,string str); 
	void setBond();
	const double g = 9.81;
/***********-- read information from dump file ----****************/
	int stage, b_flag_old;
    fndrec(infile,"DEM");
    infile>>dem.cTime;
    fndrec(infile,"COMPACT");
    infile>>stage;
    infile>>compact.l_org;
    fndrec(infile,"CONTAINER");
	infile>>box.th>>box.bh;
	if (compact.stage == 0 )
	 {
		 for(int j=0;j<3;j++)
			{
				box.thv[j] = 0.0;
				box.bhv[j] = 0.0;
			}
	 }
	
	else if ((stage != 0) || (compact.stage != 1) ) 
	 {
		 for(int j=0;j<3;j++)
			infile>>box.thv[j]>>box.bhv[j];
	 }
	fndrec(infile,"BONDINITIAL");
    infile>>b_flag_old;
	
	fndrec(infile,"PARTICLE");
    infile>>np;
    setArray(1);
		
	  for(int i=0;i<np;i++) 
		{
			for(int j=0;j<3;j++)
			{
		    infile>>pPos[i][j]>>pVel[i][j]>>pAngp[i][j]>>pAngv[i][j];

	           for(int k=0;k<box.nWall;k++)
			   {infile>>dispt_pw[i][k][j];}
			}
		     infile>>pDia[i]>>pType[i]>>nPoint[i];
		}
    
    fndrec(infile,"NEIGHBOR");
      infile>>nne;
     setArray(2);

	for(int i=0;i<nne;i++) 
		{
			infile>> nList[i] >> rbond[i] >>lqdbrg[i]>>dispAn_bond[i]>>ch_nDsp_pp[i]>>ch_nDsp_mx_pp[i]>>ch_plstDfm_pp[i]>>ch_plstRad_pp[i];//>>initalGap_bond[i];
            for(int j=0;j<3;j++)
			{infile>>dispt_pp[i][j]>>dispAt_bond[i][j]>>dispt_bond[i][j];}
		}

    cout<<"read dump file finish."<<endl;

	if ((stage == 0) & (compact.stage == 1) ) 
	   {
			//  reduce box velocity
	  		for(int i=0;i < 3;i++ )
			{
			  box.thv[i] = box.thv[i] / sqrt(pSize * g);
			  box.bhv[i] = box.bhv[i] / sqrt(pSize * g);
			}
			box.th = -9999;
		for(int i=0; i<npTotal ; i++)
		  {
            if (box.th < ( pPos[i][2]+0.5*pDia[i]) )
            box.th = pPos[i][2]+0.5*pDia[i];
		  }
			compact.l_org = box.th - box.bh ;
	   }
	compact.dsp_mx = compact.l_org * compact.strain_mx;
	if (b_flag == 1 && b_flag_old == 0 )
      { 
		  setBond();
	  }
}

void getvGapMx(double a, double hmn, double &hmx)

{
    double f, massMn,sizeMn,retarded_factor;
    hmx = 0.0;
    massMn  = Min(pMass,npTotal);
    sizeMn  = Min(pDia,npTotal)/4.0 ;
    if (a == 0.0)
	  hmx = 0.0;
	else
		{
			do
			{  hmx = hmx+ hmn;
			retarded_factor = 1.0 - 1.0 /(1.0 + 1.0e-7/(5.32*hmx*dem.rlunit));
			f = 6.0* a * sizeMn/(hmx*hmx) * retarded_factor/massMn;
			}while(f > 1.0);
			if (hmx > 0.5) 
				hmx = 0.5;
		}
}
void setPV()
//---- set agg position and velocity
    
{   
    double vunit;
    vunit = sqrt(pSize * 9.81);
    impact.vel[0] = -impact.vMag * cos(impact.angle * 3.14159 / 180.0) / vunit;
    impact.vel[1] = 0.0;
    impact.vel[2] = -impact.vMag * sin(impact.angle * 3.14159 / 180.0) / vunit;
       
		for(int ip=0;ip < np;ip++)
		{
			   for(int   k=0;k<3;k++)   
		    {
				pPos[ip][k] = pPos[ip][k] + impact.pos[k];
				pVel[ip][k] = impact.vel[k];
			}
		}

}
void setBond()
	{
    //----  if two particles contact together, then bonded them ------
		double gap,ipRad,jpRad,ijMag,ipPos[3],jpPos[3],ijPos[3];
		for(int  ip=0;ip < np; ip++)
			{  
				ipRad   = 0.5 * pDia[ip];
				 for (int k= 0; k < 3; k++)
					ipPos[k] = pPos[ip][k];
			 			
			 for(int ine = nPoint[ip];ine< nPoint[ip+1];ine++)
				{
					int jp = nList[ine];
					 jpRad   = 0.5 * pDia[jp];
					for (int k= 0; k < 3; k++)
						jpPos[k] = pPos[jp][k];
					for (int k= 0; k < 3; k++)
						ijPos[k]   = ipPos[k] - jpPos[k];
					ijMag   = vMag(ijPos,3);
					gap   = ijMag - (ipRad + jpRad);
					 if (gap <= 0.0) 
						 {
							rbond[ine] = 1;
							initalGap_bond[ine]= -gap;
						 }
				}
			 }	

	 }
           	
