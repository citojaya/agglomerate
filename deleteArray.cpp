#include "gVarExt.h"
void delArray()
//delete globle array
{       deleteArray(pPos,npTotal);
        deleteArray(pVel,npTotal);
        deleteArray(pAngp,npTotal);    
        deleteArray(pAngv,npTotal);   
        deleteArray(pForce,npTotal);   
 	    deleteArray(pMom,npTotal); 
	 	deleteArray(pType);     
        deleteArray(pDia);   
        deleteArray(pMass);    
        deleteArray(pInert);   
        deleteArray(pRMom);       
        deleteArray(pCprss);   
        deleteArray(nPoint);
        deleteArray(contact_with_wall,npTotal);
        deleteArray(pInf);
        deleteArray(dispt_pw,npTotal,box.nWall);

}
