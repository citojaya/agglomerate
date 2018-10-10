#include "gVarExt.h"
   
const int grBin = 250;
const double dgr = 0.02;
double gr[grBin][2];

void getGr()
{
    double norm, radi;
    double sx, sy, szmx, szmn;
    double vol, den0, ratio;
    int iBin;
   
    double ijPos[3];
    
    sx = box.Rxy;
    sy = box.Rxy;

	for(int i=1; i< np; i++)    // minimum height
	{
		szmn = pPos[0][2];
		if(szmn >pPos[i][2])
         szmn =pPos[i][2];
	}
	for(int i=1; i< np; i++)    // maximum height 
	{
		szmx= pPos[0][2];
		if(szmx<pPos[i][2])
         szmx=pPos[i][2];
	}
	for(int i=0; i< grBin; i++)     
	 {for(int j=0; j< 2; j++)    
	   gr[i][j] = 0.0;}
  
	for(int ip=0; ip< np; ip++)   
	  {  
        for(int jp=ip+1; jp< np; jp++)
	    { 
			for(int k=0; k< 3; k++)
            ijPos[k] = pPos[ip][k] - pPos[jp][k];

            
            iBin = int(vMag(ijPos,3)/dgr)-1;
            
            if (iBin < 0) continue;
            if (iBin >= grBin) continue;
            gr[iBin][1]++;
		 } 
	   }
    vol = sx * sy * (szmx - szmn);
    den0 = np/vol;

	for(iBin=0; iBin< grBin; iBin++)     
	{
		radi = (iBin+1.5) * dgr;
        norm = 4 * PI * radi * radi * dgr;
        gr[iBin][0] = radi;
        gr[iBin][1] = 2*gr[iBin][1]/(np*norm)/den0;
	}
    //--- to unity ---
	 ratio=0.0;
    for(int i=grBin-49; i<grBin;i++ )
	{ratio += gr[i][1];}
     ratio=ratio/50;
	for(int i=0; i< grBin; i++)     
    gr[i][1] = gr[i][1] / ratio;
}
