#include <math.h>
#include <stdio.h>
#include "vector.h"
#include "Region4.h"


//Given i,j,k of our box and the level, we return its spot in the giant array

//(pow(8,t)) counts all the spots from previous levels



/*int arrayIndex_iV(intVec patsmom)
{
  int prevL=0;

  for(t=0; t<patsmom.Lvl; t++)
  {
    prevL += pow(8,t);
  }

  arrayIndex_iV = prevL + patsmom.i + pow(2,patsmom.Lvl)*patsmom.j + pow(pow(2,patsmom.Lvl),2)*patsmom.k;

  return arrayIndex_iV;
}

*/


//-------------------------------------------------------------------------------


//void Children(int Levl)
int main()
{

Region4 child;
int Levl = 5;
int s;
int t;
int q;
int g;
int m;
int n;
int k;
int var;
int shit;
    
for (s=0; s<Levl-1; s++)                //for loop looping through each level up to Levl-1, our highest level
{
    var = 0;
    shit = 8;
    for (t=0; t<pow(8,s); t++)          //loops through the regions in each lvl, starting at 0 for each
  {
    for (q=0; q<8; q++)               //loops through the eight children of each region, finding the childs i,j,k
    {
        
      n = pow(2,s);
      m = (t/8)%n;
      g = (t/64)%n;
        
      child.i = q%2 + (t%n)*2; //(var/k)%k; //+ 2*t - g*(4*s);
      child.j = (q/2)%2 + ((var/16)%2)*2;
        //child.j = ((var/16)%2)*2;
        child.k = (q/4)%2;
      child.Lvl = s+1;

      printf("the childs position is: (%d, %d, %d), level: %d\n", child.i,child.j,child.k,child.Lvl);
        
      var++;
        shit = shit*2;
        
        //var is a sequential positive integer (including 0) that begins at zero again for each s loop.
        
      //Region[arrayIndex(s,t)].kids[q] = arrayIndex_iV(child);
    }
  }
}
    
}
        
    


















