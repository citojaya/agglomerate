#include <iostream>
#include <fstream>
#include <string>
#include "gVarExt.h"

void fndrec(ifstream &infile,string str); //declare fountion

/********** Load Data **********************************************/
/*******  prompt for the generic file name *************************/
void loadData()
{
ifstream infile;
  while(1)
   {
    cout<<"Please enter the generic file name >";
    string filename;
	cin>>genfile;
        filename=genfile+".in";
	infile.open(filename.c_str(),ios::in);
	if(!infile)
	  {
	   cout<<"Unable to open the input file: "
	   <<filename<<"\n"<<"Please try again."<<endl;
	   infile.clear();
          }
       else break;
	}
//container information
    box.loadData(infile);
//particle information (including wall)
	fndrec(infile,"REALSIZE");
    infile>> pSize;
    fndrec(infile,"SIZEDISTRIBUTION");
    infile>> sd;
    fndrec(infile,"MATERIALTYPE");
    infile>>ntype;

	pMat= new MAT[ntype+1];//dynamic initialize 
	npTotal=0;   //particle total number
	   for (int i=0;i<=ntype;i++)
	    {
		 pMat[i].loadData(infile);
		 npTotal+=pMat[i].np;
	    }
    fndrec(infile,"INITIALVELOCITY");
    infile>>initVel;
    fndrec(infile,"CUTGAP");
    infile>>cutGap;
    fndrec(infile,"VGAPMIN");
    infile>>vGapMn;	
	compact.loadData(infile);
//bond parameter
    fndrec(infile,"BOND");
	infile>>rad_multi;
	infile>>b_emod_n;
	infile>>ratio_stiff;
	infile>>max_b_tensile>>max_b_shear;
	infile>>b_flag;
// capillary force parameter
	lqd.loadData(infile);

//simulation paramters
	dem.loadData(infile);
    fndrec(infile,"IMPACT");
     for(int   k=0;k<3;k++) 
	 {infile>>impact.pos[k];}

     infile>>impact.vMag>>impact.angle;

	infile.close();

}

