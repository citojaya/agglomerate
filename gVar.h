#include <string>
#include "gVarclass.h"
using namespace std;
#ifndef GVAR_H
#define GVAR_H

////declare the global function
void loadData();
void setInitial();
void motion();
void outfile(int);
void delArray();

/********* global variables ********************************************************/
//container  
    CONTAINER box;
//Simulation parameters
	time_t start_time;
	double etime;
    int np, npTotal;
    int nne;
    int ntype;
    double pSize;
    double cutGap;
    double vGapMx, vGapMn;
    double initVel;
    double mxOlp;
    SIMULATION dem;
	YIELDTYPE yld;
	LIQUIDTYPE lqd;
	double tensile_stress, tensile_vol;
    int tensile_cnv, tensile_cnc;
//parameters for compaction
	COMPACTION compact;
// parameters for Mat
  MAT * pMat;
//parameters for impact
  IMPACT impact;
//parameters for bond
    double rad_multi,b_emod_n, ratio_stiff, max_b_tensile,max_b_shear,bGapmx;
	int b_flag; 
//Other parameters
	char sd;
    const int fl_size =10 ;
	const double PI=3.141592653589793238462643383279506;
	string genfile;
//zone
	int nxzone, nyzone, nzzone;
	ZONETYPE ***zone;

/********* global array ********************************************************/
	double **pPos, **pVel, **pForce;
    double **pAngp, **pAngv, **pMom;
    double *pRMom, *pDia, *pMass, *pInert, *pCprss;
	int *pType;
    double ***dispt_pw;
    double **wpfList;
    double ***wForce;
    vector<vector<double> > dispt_pp;
    vector<vector<double> >fList;
	vector<int> nList;

    vector<double> dispAn_bond; // bond normal angular displacement 
    vector<vector<double> > dispAt_bond; // bond tangential angular displacement 
    vector<vector<double> > dispt_bond; // bond tangentialdisplacement 
    vector<double> initalGap_bond; // bond initial gap 
	vector<int> rbond; //bond flag
	vector<int> lqdbrg; //liquid flag
	// varations for contact_history
	vector<double> ch_nDsp_pp,ch_nDsp_mx_pp,ch_plstDfm_pp,ch_plstRad_pp; //normal and max normal displacement,plastic deformation
	double **ch_nDsp_pw,**ch_nDsp_mx_pw,**ch_plstDfm_pw,**ch_plstRad_pw;  //normal and max normal displacement,plastic deformation

	int *nPoint;
    int ** contact_with_wall;
    int* cn;
    pInfo *pInf;
    double oPos[3];          // original point (mass centre, or predifined point)


#endif