/**********************************************
  The force is what gives a jedi his power.  It's an energy field created by all living things.  It surrounds us and penetrates us.  It binds the galaxy together.  Except that Meghan woman, she's a butt.
 **********************************************/

// Aaron is a poopy-head. --WF
//He's referring to you, Aaron Buttchins
//Shut up Megbutt

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "planet.h"
#include "pln.h"
#include "region.h"
#include "intVec.h"
#include "wellSepCoM.h"
#include "displist.h"
#define DEBUG false
extern double G;


/**********************************************
  Finds net force on planet "myplanet" from region "regionToTest"
 **********************************************/

void space(int n)
{
  for (int i=0;i<=n;i++) if (DEBUG)printf("%d%d%d ",i,i,i);
}

void dash(int n)
{
  for (int i=0;i<n;i++) if (DEBUG)printf("--- ");
}

void forceMagic(region regionToTest, planet& myplanet, planet *BD){
  //assuming implementation in recover.c, N is number of planets
  int i,j,k,a;
  double rad;
  vector F;
  intVec my, other;
  pln *currPln;
  planet current;
  int lev=regionToTest.level;

  if (DEBUG)
  {
  space(lev);
  printf("BEGIN computing force on planet from region %d:%d,%d,%d containing %d planets : ",regionToTest.level,regionToTest.i,regionToTest.j,regionToTest.k,regionToTest.numPln); displist(regionToTest.planets);
  }
  if(regionToTest.level == 4 || (regionToTest.numPln==1))
  {
    currPln=regionToTest.planets;
    space(lev); if (DEBUG)printf("exactly, by iterating over planets.\n");
    while (currPln != NULL)
    {
      space(lev); if (DEBUG){printf("Remaining list of planets in this region: "); displist(currPln);}
      current=BD[currPln->plnNum];
      rad = radius(myplanet.pos,current.pos);
      if (rad > 0) // very crude check to make sure it's not the same planet pulling on itself
      {
        F = (current.pos-myplanet.pos) * G * myplanet.m * BD[currPln->plnNum].m / (rad*rad*rad) ;//Direct force calculation
        myplanet.acc = myplanet.acc + F/myplanet.m;
        space(lev); if (DEBUG)printf("Exact force from planet %d at r=(%.2e, %.2e, %.2e) on me at r=(%.2e, %.2e, %.2e) : acc is now = (%.2e, %.2e, %.2e)\n", currPln->plnNum, 
            current.pos.x, current.pos.y, current.pos.z,
          myplanet.pos.x, myplanet.pos.y, myplanet.pos.z,
          myplanet.acc.x, myplanet.acc.y, myplanet.acc.z);
      }
      else
      {
        space(lev); if (DEBUG)printf("This is myself, skipping\n");
      }
      currPln=currPln->nextPln;
   }
  }

  else
  {
    //		intVec my = getIJKVec(n, myplanet.pos);  //intVec to store coordinates of region containing myplanet

    k = wellSepCoM(myplanet,regionToTest);
    k=0;
    if (k==1)
    {
      space(lev); if (DEBUG)printf("region is well-separated, computing force via multipole\n");
      rad = radius(myplanet.pos,regionToTest.com);
      vector Fhere;
      Fhere = G * myplanet.m * regionToTest.mass * (regionToTest.com-myplanet.pos) / (rad*rad*rad);
      myplanet.acc = myplanet.acc + Fhere/myplanet.m;
      space(lev); if (DEBUG)printf("   --> Force from multipole is (%.2e, %.2e, %.2e), mass is %.2e, acceleration is now (%.2e, %.2e, %.2e)\n",
                    Fhere.x, Fhere.y, Fhere.z, myplanet.m,
                    myplanet.acc.x, myplanet.acc.y, myplanet.acc.z);
    }

    if (k==0) 
    {
      space(lev); if (DEBUG)printf("region is not well-separated, descending to children\n");
      for(int l=0;l<8;l++) 
      {
	if(regions[regionToTest.kids[l]].numPln == 0)
	{
//          space(lev); printf("child %d is region %d, containing %d planets: empty\n",l,regionToTest.kids[l],regions[regionToTest.kids[l]].numPln);
	  continue;
	}
	else
        {
          space(lev); if (DEBUG)printf("child %d is region %d, containing %d planets: recursing\n",l,regionToTest.kids[l],regions[regionToTest.kids[l]].numPln);
	  forceMagic(regions[regionToTest.kids[l]],myplanet,BD);
        }
      }
    }
  }
  if (DEBUG) {
  space(lev);printf("END   force on planet from region %d:%d,%d,%d containing %d planets : ",regionToTest.level,regionToTest.i,regionToTest.j,regionToTest.k,regionToTest.numPln); displist(regionToTest.planets);
  }
}	
