#include "gVarExt.h"
    const int  frqBin = 15;
    double freq[frqBin+1];
void getfrq()
{
	for (int i=0; i<= frqBin; i++)
        freq[i] = 0.0;
	for (int ip=0; ip< np; ip++)
        freq[pInf[ip].cn] ++ ;
	for (int i=0; i<= frqBin; i++)
        freq[i] = freq[i] / np;
 
}
