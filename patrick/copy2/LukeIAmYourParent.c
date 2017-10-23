#include <math.h>
#include "vector.h"
#include "planet.h"
#include "pln.h"
#include "region.h"
#include "Region4.h"
#include "getArrayIndex.h"
#include <stdlib.h>
#include <stdio.h>

//Given i,j,k of our box and the level, we return its spot in the giant array

//(pow(8,t)) counts all the spots from previous levels



int arrayIndex_R4(Region4 patsmom)
{
  int prevL=0;
  int t;

  for(t=0; t<patsmom.Lvl; t++)
  {
    prevL += pow(8,t);
  }

  return prevL + patsmom.i + pow(2,patsmom.Lvl)*patsmom.j + pow(pow(2,patsmom.Lvl),2)*patsmom.k;
}




//-------------------------------------------------------------------------------


void Children(int Levl)
{

Region4 child;
int s, t, q, g, m, n, var;

for (s=0; s<=Levl-1; s++)            //Clear lists of planets FIXME Do we need to do this here? 
    for (t=0; t<pow(8,s); t++)      
      regions[getArrayIndex(s,t)].planets=NULL;
    
for (s=0; s<Levl-1; s++)                //for loop looping through each level up to Levl-1, our highest level
{
    var = 0;
    for (t=0; t<pow(8,s); t++)          //loops through the regions in each lvl, starting at 0 for each
    {
      printf("!Babyfinder: level %d region %d out of %.0f\n",s,t,pow(8,s));
      for (q=0; q<8; q++)               //loops through the eight children of each region, finding the childs i,j,k
      {
        
        n = pow(2,s+1);
        m = pow(4,s+1);
        g = (t/(pow(2,s)));
        
        child.i = ((q)%2) + 2*t - g*(4*s);
        child.j = (var/n)%n;
        child.k = (var/m)%m;


        child.Lvl = s+1;

        printf("!  Child %d lives at: %d %d %d, level: %d\n", q,child.i,child.j,child.k,child.Lvl);
        
        var++;
        
        //var is a sequential positive integer (including 0) that begins at zero again for each s loop.
        
        regions[getArrayIndex(s,t)].kids[q] = arrayIndex_R4(child);
      }
    }
}
    
}
        
    


















