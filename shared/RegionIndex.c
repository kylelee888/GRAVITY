#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"


int main()
{

int BoxLevel
{
  int i;
    
for(i=0,i<=pow(8,N),i++)
{
  lvl = i/(pow(8,N));
}
return lvl;
    
printf("we are in level %d", lvl);
}




vector BoxLocation()
{
  double L = 4; //length
  int N = 2; //level
  int i;
    
    vector pos;
    for (i=0,i<=pow(8,N),i++)
    {

      pos.x = i*(L/pow(2,N));
      pos.y = (i/4)*(L/pow(2,N));
      pos.z = (i/16)*(L)/pow(2,N));
    }
    return pos;
}
    
    
    
    
}