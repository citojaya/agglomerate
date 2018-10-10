#include "gVarclass.h"
#include<ctime>
#ifndef GVAREXT_H
#define GVAREXT_H
using namespace std;

//function template
//make 2DArry function 

template<class   T>   
  void   makeArray(T   ***&x,int   rows,int   cols, int z)   
  {      
         x=new   T**[rows];   
	    for(int   i=0;i<rows;i++)  
	   {
         x[i] = new   T*[cols];   
         for(int   j=0;j<cols;j++)   
		  { x[i][j]=new   T[z];   
          for (int k = 0; k < z; k++)
              x[i][j][k]=0;
		  }
	   }         
  };


template<class   T>   
  void   makeArray(T   **&x,int   rows,int   cols)   
  {   
         x=new   T*[rows];   
         for(int   i=0;i<rows;i++)   
		 { x[i]=new   T[cols];   
         //initial 
         for (int k = 0; k < cols; k++)
            x[i][k]=0;
		  }
               
  };
   template<class   T>  
  void makeArray(T   *&x, int cols ) 
  {
           x=new   T[cols]; 
		   //initial
      for (int i = 0; i < cols; i++)
               
             x[i]= 0;
   };
  template<class   T>   
  void   deleteArray(T   ***&x,int   rows, int cols)   
  {   
         for(int   i=0;i<rows;i++)  
		 {
			 for(int   j=0;j<cols;j++) 
			 {
			  delete[] x[i][j];
			 }
		  delete [] x[i];
		 }
			 
		 delete[] x;
         //x= NULL;
  };
  template<class   T>   
  void   deleteArray(T   **&x,int   rows)   
  {   
         for(int   i=0;i<rows;i++)   
		 { 
			 delete[] x[i];
		 }
			 
		 delete[] x;
         //x= NULL;
  };
  template<class   T>   
  void   deleteArray(T   *&x)   
  {   
		 delete[] x;
        // x= NULL;
  };

 

template   <typename   T>     // two 
  void   ReSize(vector<vector<T>   >   &   vec,int   SizeX,int   SizeY)   
  {   
  vec.resize(SizeY);   
  for(int   i=0;i<SizeY;i++)   
  vec[i].resize(SizeX); 
  };


//return   the   index   of   an   array's   maximum   
  template<class   T>   

  T  Max  ( T   *&a , int n)   
  {   

        T temp=a[0];

        for   (int   i   =   1;   i   <   n;   ++i)   
        {   
                if   (temp   <  a[i])   
                {   
					temp  =   a[i]; 
                }   
        }     
          
        return   temp;   
  };
 
  //return   the   index   of   an   array's   minimum   

  template<class   T>   

  T  Min   ( T   *&  a, int n)   
  {   

        T temp=a[0];

        for   (int   i   =   1;   i   <   n;   i++)   
        {   
                if   (temp   >   a[i])   
                {   
					temp  =   a[i]; 
                }   
        }     
          
        return   temp;   
  };

  template<class   T>   
  void  Mutiple   ( vector<T>  &a, T Num)   
  {   
        int n;
        n=a.size();

        for   (int   i   =   0;   i   <   n;   i++)   
        {   
          a[i]=a[i]*Num;     
		   
        }     
          
           
  };

//declare globle function
void dump();
void pclock() ;
void oufile();
double  vMag(double *v, int n);

/********* global variables ********************************************************/

//container  
    extern CONTAINER box;
//Simulation parameters
    extern double etime;
    extern int np, npTotal;
    extern int nne;
    extern int ntype;
    extern double pSize;
    extern double cutGap;
    extern double vGapMx, vGapMn;
    extern double initVel;
    extern double mxOlp;
    extern SIMULATION dem;
	extern YIELDTYPE yld;
	extern LIQUIDTYPE lqd;
    extern double tensile_stress, tensile_vol;
    extern int tensile_cnv, tensile_cnc;
//parameters for compaction
	extern COMPACTION compact;
// parameters for Mat

    extern MAT *pMat;
// parameters for impact
	extern IMPACT impact;
//parameters for bond
    extern double rad_multi,b_emod_n, ratio_stiff, max_b_tensile,max_b_shear,bGapmx;
	extern int b_flag; 
//Other parameters
	extern char sd;
	const int fl_size =10 ;
    const double PI=3.141592653589793238462643383279506;
    extern string genfile;
//zone
    extern int nxzone, nyzone, nzzone;
    extern ZONETYPE ***zone;


/********* global array ********************************************************/
	extern double **pPos, **pVel, **pForce;
    extern double **pAngp, **pAngv, **pMom;
    extern double *pRMom, *pDia, *pMass, *pInert, *pCprss;
	extern int *pType;
    extern double ***dispt_pw;
    extern double **wpfList;
    extern double ***wForce;

    extern vector<vector<double> > dispt_pp;
    extern vector<vector<double> >fList;
	extern vector<int> nList; 

    extern vector<double> dispAn_bond; // bond normal angular displacement 
    extern vector<vector<double> > dispAt_bond; // bond tangential angular displacement 
    extern vector<vector<double> > dispt_bond; // bond tangential angular displacement 
    extern vector<double> initalGap_bond; // bond initial gap 
	extern vector<int> rbond; //bond flag
	extern vector<int> lqdbrg; //liquid flag

	// varations for contact_history
	extern vector<double> ch_nDsp_pp,ch_nDsp_mx_pp,ch_plstDfm_pp,ch_plstRad_pp; //normal and max normal displacement,plastic deformation
	extern double **ch_nDsp_pw,**ch_nDsp_mx_pw,**ch_plstDfm_pw,**ch_plstRad_pw;  //normal and max normal displacement,plastic deformation
	
    extern int* nPoint;

   extern  int ** contact_with_wall;


   extern  int* cn;
    extern pInfo *pInf;

   extern  double oPos[3];          // original point (mass centre, or predifined point)

	#endif

