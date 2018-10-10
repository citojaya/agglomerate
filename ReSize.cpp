#include <vector>
  
using namespace std;

template   <typename   T>     
  void   ReSize(vector<vector<T>   >   &   vec,size_t   SizeX,size_t   SizeY)   
  {   
  vec.resize(SizeX);   
  for(int   i=0;i<SizeX;i++)   
  vec[i].resize(SizeY); 
  };

template   <typename   T>     
  void   ReSize(vector<vector<T>   >   &   vec,size_t   SizeX,size_t   SizeY, T initial)   
  {   
vec.resize(SizeX);   
  for(int   i=0;i<SizeX;i++)   
  vec[i].resize(SizeY); 
  for(int i=0;i<SizeX;i++)  
  for(int j=0;j<SizeY;j++)  
   vec[i][j]=Initial;  
  };
 






