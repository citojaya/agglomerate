#include <ctime>
//  returns the CPU time elasped
extern time_t start_time;
extern double etime;
void pclock() 

 {  
	 time_t curent_time =time(NULL);
     etime=double(curent_time-start_time);
 }

