#include <ctime>
#include "gVar.h"
using namespace std;

int  main()
{
start_time=time(NULL);
loadData();        //   read in data
setInitial()  ;     // check restart file 
motion();          // carry out the simulation
outfile(1) ;         // write out results
delArray();        //free mem
return 0;
exit(1);
}


