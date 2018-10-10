#include "gVarExt.h"
#include <cmath>
#include <iostream>
using namespace std;

   static const double dzone = 1.3;     //must be larger than necutoff
   static vector<vector<double> > pZone;
   int neZone[27];

void setzone()
{ 
    //---- allocate zone number ----
    nxzone = int((2*box.Rxy / dzone)) + 2;
    nyzone = int((2*box.Rxy / dzone)) + 2;
    nzzone = int((box.th - box.bh) / dzone) + 2;
	double x_max=pPos[0][0];
	double x_min=pPos[0][0];
	double y_max=pPos[0][1];
	double y_min=pPos[0][1];
	double z_max=pPos[0][2];
	double z_min=pPos[0][2];
	for(int ip=0; ip<np; ip++)
	{
		if (x_max < pPos[ip][0]) x_max=pPos[ip][0];
		if (x_min > pPos[ip][0]) x_min=pPos[ip][0];
		if (y_max < pPos[ip][1]) y_max=pPos[ip][1];
		if (y_min > pPos[ip][1]) y_min=pPos[ip][1];
		if (z_max < pPos[ip][2]) z_max=pPos[ip][2];
		if (z_min > pPos[ip][2]) z_min=pPos[ip][2];
	}
    x_max+=4;x_min-=4;
    y_max+=4;y_min-=4;
    z_max+=4;z_min-=4;

	x_max=min(x_max,box.Rxy);
    x_min=max(x_min,-box.Rxy);
	y_max=min(y_max,box.Rxy);
	y_min=max(y_min,-box.Rxy);
	z_max=min(z_max,box.th);
	z_min=max(z_min,box.bh);

	nxzone = int((x_max-x_min)/ dzone) + 2;
    nyzone = int((y_max-y_min)/ dzone) + 2;
    nzzone = int((z_max-z_min)/ dzone) + 2;

        zone=new   ZONETYPE**[nxzone];   
	    for(int   i=0;i<nxzone;i++)  
	   {
         zone[i] = new   ZONETYPE*[nyzone];   
         for(int   j=0;j<nyzone;j++)   
		  { zone[i][j]=new   ZONETYPE[nzzone];}
	   }          

    //---- zone node (not center) position ----
    for (int ix = 0; ix<nxzone; ix++)
	{
		for (int iy = 0; iy<nyzone; iy++)
		{
			for (int iz = 0; iz<nzzone; iz++)
			{
				zone[ix][iy][iz].x = ix * dzone + x_min;
                zone[ix][iy][iz].y = iy * dzone + y_min;
                zone[ix][iy][iz].z = iz * dzone + z_min;
                
                zone[ix][iy][iz].x = min(x_max, zone[ix][iy][iz].x);
                zone[ix][iy][iz].y = min(y_max, zone[ix][iy][iz].y);
                zone[ix][iy][iz].z = min(z_max, zone[ix][iy][iz].z);
                
                zone[ix][iy][iz].porosity = 1.0;
                zone[ix][iy][iz].np = 0;
                for (int i=0; i<mxnp; i++)
				zone[ix][iy][iz].pIndx[i] = 0;
                for(int i=0; i<3; i++)
				{
				  zone[ix][iy][iz].force[i] = 0.0;
                  zone[ix][iy][iz].velocity[i] = 0.0;
				}
			}
		}
	}

}
void set_p2z()
{   
   pZone.resize(np);
   for(int i=0; i<np; i++)
   {pZone[i].resize(3);}
   int ix,iy,iz; 
   for(int ip=0; ip<np; ip++)
   {
	    ix = int((pPos[ip][0] - zone[0][0][0].x) / dzone) + 1;
        iy = int((pPos[ip][1] - zone[0][0][0].y) / dzone) + 1;
        iz = int((pPos[ip][2] - zone[0][0][0].z) / dzone) + 1;

        zone[ix][iy][iz].np = zone[ix][iy][iz].np + 1;
        if (zone[ix][iy][iz].np > mxnp)
		{
			cout<<"zone.np larger than zone.mxnp"<<endl;
		    exit(1);
		}
        zone[ix][iy][iz].pIndx[zone[ix][iy][iz].np-1] = ip;
		pZone[ip][0] = ix;
		pZone[ip][1] = iy;
		pZone[ip][2] = iz;
   }
}

//---- porosity distribution ----
void getZPor()
{  
    double  pVol(double r1,double r2, double d);
	double ipRad;
    double zRad, zVol, zPor;
    ZONETYPE izone, jzone;
    double zPos[3], ijPos[3];
	//---- establish the zone mesh ----
       setzone();
    //---- assign particle to zones ----
    set_p2z();
   //---- get zone related properties ----
    zRad = (0.5 * dzone) * sqrt(3.0);  //radius of sphere surrounding the cell
  for (int ix = 0; ix<nxzone; ix++)
   {
	for (int iy = 0; iy<nyzone; iy++)
	 {
	    for (int iz = 0; iz<nzzone; iz++)
	   {
		  zVol = 4.0/3.0 * pow(zRad,3);
		   if (ix == 0 || ix == nxzone-1) 
	 	  	zVol = 0.5 * zVol;
		   if (iy == 0 || iy == nyzone-1) 
			zVol = 0.5 * zVol;
		   if (iz == 0 || iz == nzzone-1) 
			zVol = 0.5 * zVol;
		  
        izone = zone[ix][iy][iz];
        zPos[0] = izone.x;
        zPos[1] = izone.y;
		zPos[2] = izone.z;
        zPor = 1.0;
        for( int jjx = ix-2; jjx<=ix+2 ; jjx++)
		{
			 if(jjx<0 || jjx>nxzone-1) continue; 
          for(int jjy =iy-2; jjy<=iy+2; jjy++)
		  {
			 if(jjy<0 || jjy>nyzone-1) continue; 
            for(int jjz = iz-2;jjz<=iz+2 ;jjz++)
			{
			 if(jjz<0 || jjz>nzzone-1) continue; 
			 jzone = zone[jjx][jjy][jjz];
              for(int idx=0; idx<jzone.np; idx++)
			  {
                int ip = jzone.pIndx[idx];
                for(int k=0; k<3;k++)
				ijPos[k] = pPos[ip][k] - zPos[k];

                //if (box.nWall < 3) ijPos[0] = ijPos[0] - int(ijPos[0] / box.hflx) * box.lx;
               // if (box.nWall < 2) ijPos[1] = ijPos[1] - int(ijPos[1] / box.hfly) * box.ly;
                ipRad = 0.5*pDia[ip];
                zPor = zPor - pVol(zRad, ipRad, vMag(ijPos,3))/zVol;
			  }
			}
		  }
		}
        zone[ix][iy][iz].porosity = zPor;
		}
	  }    
    }
}

