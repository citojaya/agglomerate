#include <iostream>
#include <fstream>
#include <string>
using namespace std;
//find the position of a record
void fndrec(ifstream &infile,string str)
{
	 char ch_temp[700];
	 string str_temp;

	 infile.seekg(0,ios::beg);  //move the file's pointer to begain
     while(infile.getline(ch_temp,700))
	 {
   		str_temp=ch_temp;
		if (str_temp==str) break;
	 }
     if (str_temp!=str) cout<<"Unable to find record"<<" "<<str<<endl;
}



