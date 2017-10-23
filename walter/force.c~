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
#include "timing.h"
extern double G;
extern int exact;
extern int approx;
int DEBUG;
extern int LVL;
extern long int approx_time, exact_time;
extern int coll1[], coll2[];
extern int collnum;
extern bool collide_check;

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

void forceMagic(region& regionToTest, planet& myplanet, planet *BD, double alpha2){
  //assuming implementation in recover.c, N is number of planets
   int i,j,a;
   double k;
   double rad;
   vector F;
   pln *currPln;
   planet current;
   int lev;
  lev=regionToTest.level;

  if(regionToTest.level == LVL-1 || (regionToTest.numPln==1))
  {
    starttimer(3);
    currPln=regionToTest.planets;
    while (currPln != NULL)
    {
      current=BD[currPln->plnNum];
      rad = radius(myplanet.pos,current.pos);
      if (rad > 0) // very crude check to make sure it's not the same planet pulling on itself
      {
        if (collide_check) 
        {
//          printf("!Collide check on: looking at %d/%d, radius %f, sum of sizes %f\n",myplanet.num,current.num,rad,current.r + myplanet.r);
          if ((current.r + myplanet.r > rad*2) && (myplanet.num > current.num))
          {
            coll1[collnum] = myplanet.num;
            coll2[collnum] = current.num;
            collnum++;
          }
        }
        F = (current.pos-myplanet.pos) * G * myplanet.m * BD[currPln->plnNum].m / (rad*rad*rad) ;//Direct force calculation
        myplanet.acc = myplanet.acc + F/myplanet.m;
        exact++;
      }
      currPln=currPln->nextPln;
   }
   exact_time+=stoptimer(3);   
  }

  else
  {
    k = wellSepCoM(myplanet,regionToTest);
    if (k > alpha2)
    {
//      printf("!Planet at %.2e,%.2e,%.2e well separated from region at %.2e,%.2e,%.2e: k = %e, alpha2 = %e\n",
//            myplanet.pos.x,myplanet.pos.y,myplanet.pos.z,
 //           regionToTest.com.x,regionToTest.com.y,regionToTest.com.z,
  //          k,alpha2);
      rad = radius(myplanet.pos,regionToTest.com);
      vector Fhere;
      myplanet.acc = myplanet.acc + regionToTest.mass * (regionToTest.com-myplanet.pos) / (rad*rad*rad);
      approx++;
    }

    else
     {
      for(int l=0;l<8;l++) 
      {
	if(regions[regionToTest.kids[l]].numPln == 0)
	{
	  continue;
	}
	else
        {
	  forceMagic(regions[regionToTest.kids[l]],myplanet,BD,alpha2);
        }
      }
    }
  }
}	
