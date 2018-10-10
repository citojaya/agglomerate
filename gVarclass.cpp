#include <iostream>
#include <fstream>
#include <string>
#include<cmath>
using namespace std;
#include "gVarclass.h"

void fndrec(ifstream &infile,string str);

// class function

void CONTAINER::loadData(ifstream &infile)
{   
	fndrec(infile,"BOX");
    infile>>nWall;
	infile>>Rxy>>th>>bh;
	infile>>gRxy>>gth>>gbh;
	infile>>thv[0]>>thv[1]>>thv[2];
	infile>>bhv[0]>>bhv[1]>>bhv[2];
	nc=0;
}
void SIMULATION::loadData(ifstream &infile)
{   
	fndrec(infile,"SIMULATION");  
    infile>>simTime>>cTime>>chkTime;
    infile>>fdpert>>tppert>>anapert>>dumppert;
    infile>>dtFactor;
    infile>>ctrF;

}
void COMPACTION::loadData(ifstream &infile)
{   
	fndrec(infile,"COMPACTION");
	infile>>stage;
    infile>>strain_mx;
	infile>>stress_mx;
}
void MAT::loadData(std::ifstream &infile)
{
	infile>>np>>dia>>density>>ha>>ymod>>pois>>sfrc>>rfrc>>dmpn>>yldp;
	emod = ymod / (1 - pow(pois,2));
    
}
void LIQUIDTYPE::loadData(ifstream &infile)
{   
	fndrec(infile,"LIQUIDPROPERTY");
	infile>>sTension>>lqV>>cGapMn>>theta;
	infile>>flag;

}


