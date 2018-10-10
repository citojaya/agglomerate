/************************************************************************   
*	Building the neighbor list:                                         *
*	Algorithm: neighborList is one dimension array holding neighbor     *
*		and neighborPoint is a array holding pointers                   *
*		the neighbor of particle ip in neighborList are from            *
*		neighborPoint(ip) -> neighborPoint(ip+1)-1                      *
*************************************************************************/
#include <iostream>
#include <cmath>
#include "gVarExt.h"
using namespace std;

void  bldNList()
{
    double ipRad, jpRad;
    double neGap;
    double  ipPos[3], jpPos[3], ijPos[3];

    int *nPointOld;
	vector<int> nListOld;
	vector<int> rbondOld; 
	vector<int> lqdbrgOld;
    vector<vector<double> > dispt_ppOld;
    vector<double> dispAn_bondOld;
	vector<double> initalGap_bondOld;
    vector<vector<double> > dispAt_bondOld, dispt_bondOld;
	vector<double> ch_nDsp_ppOld,ch_nDsp_mx_ppOld,ch_plstDfm_ppOld,ch_plstRad_ppOld; //normal and max normal displacement,plastic deformation

	//--- allocate array
	makeArray(nPointOld,np);
	nListOld.resize(nne);
	rbondOld.resize(nne);
	lqdbrgOld.resize(nne);
    ReSize(dispt_ppOld,3,nne);
	dispAn_bondOld.resize(nne);
	initalGap_bondOld.resize(nne);
    ReSize(dispAt_bondOld,3, nne);
    ReSize(dispt_bondOld,3, nne);
	ch_nDsp_ppOld.resize(nne);
	ch_nDsp_mx_ppOld.resize(nne);
	ch_plstDfm_ppOld.resize(nne);
	ch_plstRad_ppOld.resize(nne);

      for(int i=0;i<np;i++)   //np  is the particals' number
	  {
		nPointOld[i] = nPoint[i];
        nPoint[i] = 0;
	  }

      for(int i=0;i< nne;i++)      //nne is size of nList
	  {
        nListOld[i] = nList[i];
        nList[i] = 0;
		rbondOld[i] = rbond[i];
		rbond[i] = 0 ;
		lqdbrgOld[i] = lqdbrg[i];
		lqdbrg[i] = 0;
		dispAn_bondOld[i] = dispAn_bond[i];
		dispAn_bond[i] = 0.0;
		initalGap_bondOld[i]= initalGap_bond[i];
		initalGap_bond[i] = 0.0;
		ch_nDsp_ppOld[i] = ch_nDsp_pp[i];
		ch_nDsp_pp[i] = 0.0 ;
		ch_nDsp_mx_ppOld[i] = ch_nDsp_mx_pp[i];
		ch_nDsp_mx_pp[i] = 0.0 ;
		ch_plstDfm_ppOld[i] = ch_plstDfm_pp[i];
		ch_plstDfm_pp[i] = 0.0 ;
		ch_plstRad_ppOld[i] = ch_plstRad_pp[i];
		ch_plstRad_pp[i] = 0.0;
          for(int j=0;j< 3;j++)
		  {
		    dispt_ppOld[i][j] = dispt_pp[i][j];
            dispt_pp[i][j] = 0.0;
			dispAt_bondOld[i][j] = dispAt_bond[i][j];
			dispAt_bond[i][j] = 0.0;
			dispt_bondOld[i][j] = dispt_bond[i][j];
			dispt_bond[i][j] = 0.0;
		  }
	  }
    
 /******--- update neighbor list ---*****************************/
	  neGap = cutGap + vGapMx;  //r0 and Gap
      nne = 0;
   for(int i=0;i< np;i++)      //np is the number of particals
   {
	     nPoint[i] = nne;

		 for(int k=0;k< 3;k++)
		 {ipPos[k] = pPos[i][k];}

          ipRad = 0.5 * pDia[i];
        
		for( int j = i+1;j<np;j++)
		  {
            for(int m=0;m< 3;m++)
			{jpPos[m] = pPos[j][m];}

            jpRad = 0.5 * pDia[j];
            for(int n=0;n< 3;n++)
			{ijPos[n] = ipPos[n] - jpPos[n];}

            if ((vMag(ijPos,3)-(ipRad+jpRad)) > neGap) continue;    //sorry, not neighbor

			if (nne== nList.size()) 
			  {
				cout<<" inflate"<<endl;
				int new_size=nne+npTotal;
			    nList.resize(new_size);
			    ReSize(dispt_pp,3, new_size);
                ReSize(fList,fl_size,new_size+1); 
				rbond.resize(new_size);
				lqdbrg.resize(new_size);
				dispAn_bond.resize(new_size);
				initalGap_bond.resize(new_size);
				ReSize(dispAt_bond,3, new_size);
				ReSize(dispt_bond,3, new_size);
				ch_nDsp_pp.resize(new_size);
				ch_nDsp_mx_pp.resize(new_size);
				ch_plstDfm_pp.resize(new_size);
				ch_plstRad_pp.resize(new_size);
                cout<<"the nList new_size is "<<new_size<<endl;
			   }
			    nList[nne] = j;
				rbond[nne] = 0; 
				lqdbrg[nne] = 0;
				dispAn_bond[nne]= 0.0 ;
				initalGap_bond[nne]=0.0;
				ch_nDsp_pp[nne]= 0.0 ;
				ch_nDsp_mx_pp[nne]= 0.0 ;
				ch_plstDfm_pp[nne]= 0.0 ;
				ch_plstRad_pp[nne]= 0.0;
                for(int s=0;s< 3;s++)
			    {
					dispt_pp[nne][s] = 0.0;
					dispAt_bond[nne][s]= 0.0 ;
					dispt_bond[nne][s] = 0.0 ;
				}
            // check if j was also the old neighbor of ip
				for (int p=nPointOld[i];p<nPointOld[i+1];p++)
			    {
					//cout<<"   do it check old list!";
				  int jOld = nListOld[p];
                   if (jOld == j) 
				    {
						dispAn_bond[nne]= dispAn_bondOld[p];
						initalGap_bond[nne]=initalGap_bondOld[p];
						rbond[nne] = rbondOld[p];
						lqdbrg[nne] = lqdbrgOld[p];
						ch_nDsp_pp[nne] = ch_nDsp_ppOld[p];
						ch_nDsp_mx_pp[nne] = ch_nDsp_mx_ppOld[p];
						ch_plstDfm_pp[nne] = ch_plstDfm_ppOld[p];
						ch_plstRad_pp[nne] = ch_plstRad_ppOld[p];
						for(int q=0;q< 3;q++)
							{
							 dispt_pp[nne][q] = dispt_ppOld[p][q];
							 dispAt_bond[nne][q]= dispAt_bondOld[p][q] ;
							 dispt_bond[nne][q]= dispt_bondOld[p][q] ;
							}  
						break;
				    }
			     }
				nne++;
		 }

	 
     /*******--- contact with wall -----**********/
             for(int q=0;q< box.nWall;q++)
			 {contact_with_wall[i][q] = 0;}

			for(int iWall = 0;iWall < box.nWall;iWall++)
			 {
               if (iWall == 0)  //cotact with top or bottom
			    {   
				   if (abs(ipPos[2] - box.bh) < (ipRad+neGap)) 
                    contact_with_wall[i][iWall] = -1;
                   else if (abs(ipPos[2] - box.th) < (ipRad+neGap)) 
                    contact_with_wall[i][iWall] = 1;
			     }

               else if (iWall == 1) 
			     {
				    if ((box.Rxy - sqrt(pow(ipPos[1],2)+pow(ipPos[0],2))) < (ipRad+neGap))  
                     contact_with_wall[i][iWall] = 1;
		         }
              }
        
   }
			 delete[] nPointOld;
			 nPointOld=NULL; 
}