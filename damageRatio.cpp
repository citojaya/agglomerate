/**************************************************************************************   
*	Analyzing the number of Fragments:                                              *
*	Algorithm: app_nList is one dimension array holding neighbor                      *
*		and app_nPoint is a array holding pointers                                    *
*		the neighbor of particle ip in app_nList are from                             *
*		app_nPoin(ip) -> app_nPoin(ip+1)-1                                            *
*    Result is a array holding the number of every agglomerate with different number  *
***************************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "gVarExt.h"
using namespace std;

int *app_nPoint;
int app_sum;
vector<int> app_nList;
vector<int> app_nPointOld;

void  damageRatio()
{
    void fndpart(int ip);
    
    double ipRad, jpRad;
    double  ipPos[3], jpPos[3], ijPos[3];
    makeArray(app_nPoint,np+1);
	app_nPointOld.reserve(np);
	app_nList.resize(np);

	//--- allocate array

      for(int i=0;i<np;i++)   //np  is the particals' number
	  {
		app_nPoint[i] = 0;
		app_nList[i] = 0;
		result[i+1]=0;
	  }
    
 /******--- update neighbor list ---*****************************/
      int nnp = 0;
   for(int i=0;i< np;i++)      //np is the number of particals
   {
	     app_nPoint[i] = nnp;

		 for(int k=0;k< 3;k++)
		 {ipPos[k] = pPos[i][k];}

          ipRad = 0.5 * pDia[i];
        
		for( int j = 0;j<np;j++)
		  {
            if (i==j) continue;
            for(int m=0;m< 3;m++)
			{jpPos[m] = pPos[j][m];}

            jpRad = 0.5 * pDia[j];
            for(int n=0;n< 3;n++)
			{ijPos[n] = ipPos[n] - jpPos[n];}

			if ((vMag(ijPos,3)-(ipRad+jpRad)) > 0) continue;    //sorry, not neighbor

			if (nnp== app_nList.size()) 
			  {
				//cout<<" inflate"<<endl;
				
				int new_size=nnp+npTotal;
				app_nList.resize(new_size);
                //cout<<"the nList of agglomeration  new_size is "<<new_size<<endl;
			   }
                         
			    app_nList[nnp] = j;
				nnp++;
		 }
	     app_nPoint[np] = nnp;
   }
	// record the inter-particle bonds

   result[0]=nnp/2;
   
   for(int  ip=0;ip < np; ip++)
      {  
		  app_sum=0;
		  if (app_nPointOld.size()==0)
		  {
		  app_nPointOld.push_back(ip);
		  app_sum++;
		  fndpart(ip);
		  result[app_sum]++;
		  }
		  else
		  {
			  int Flag=0;
			  for(int i=0;i< app_nPointOld.size();i++)
		      {
				  if ( ip == app_nPointOld[i]) 
			      {
				      Flag=1;
					  break;
				  }
			  }
			  if(Flag==1) continue;
		      app_nPointOld.push_back(ip);
			  app_sum++;
			  fndpart(ip);
		      result[app_sum]++;
		 }
	  }


			delete[] app_nPoint;
			app_nPoint=NULL;
			app_nList.clear();
			app_nPointOld.clear();
}
   void fndpart(int ip)
   {
	 
	for(int ine = app_nPoint[ip];ine< app_nPoint[ip+1];ine++)
	       {
			   int EX=0;
			   for(int i=0;i< app_nPointOld.size();i++)
			  {
				if ( app_nList[ine] == app_nPointOld[i])
				{
					EX=1;
	                break;			
				}
			  }
			   if(EX==1) continue;
			  	app_nPointOld.push_back(app_nList[ine]);
                app_sum++;
			    fndpart(app_nList[ine]);
		   }
   }