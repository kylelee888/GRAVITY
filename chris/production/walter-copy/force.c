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
extern double G;

/**********************************************
  Finds net force on planet "myplanet" from region "regionToTest"
 **********************************************/

//void space(int n)
//{
//  for (int i=0;i<=n;i++) printf("! %d%d%d ",i,i,i);
//}

/*void dash(int n)
{
  for (int i=0;i<n;i++) printf("! --- ");
}*/

void forceMagic(region regionToTest, planet &myplanet, planet *BD){
  //assuming implementation in recover.c, N is number of planets
  int i,j,k,a;
  double rad;
  vector F;
  intVec my, other;
  pln *currPln;
  planet current;
  int lev=regionToTest.level;

  //space(lev); //printf("! BEGIN computing force on planet from region %d:%d,%d,%d containing %d planets : ",regionToTest.level,regionToTest.i,regionToTest.j,regionToTest.k,regionToTest.numPln); displist(regionToTest.planets);
  if(regionToTest.level == 4 || (regionToTest.numPln==1))
  {
    currPln=regionToTest.planets;
    //space(lev); //printf("! exactly, by iterating over planets.\n");
    while (currPln != NULL)
    {
      //space(lev); //printf("! Remaining list of planets in this region: "); displist(currPln);
      current=BD[currPln->plnNum];
      rad = radius(myplanet.pos,current.pos);
      if (rad > 0) // very crude check to make sure it's not the same planet pulling on itself
      {
        F = (current.pos - myplanet.pos) * G * myplanet.m * BD[regionToTest.planets->plnNum].m / (rad*rad*rad) ;//Direct force calculation
        myplanet.acc = myplanet.acc + F/myplanet.m;
      }
      currPln=currPln->nextPln;
   }
  }

  else
  {
    //		intVec my = getIJKVec(n, myplanet.pos);  //intVec to store coordinates of region containing myplanet

    k = wellSepCoM(myplanet,regionToTest);

    if (k==1)
    {
      //space(lev); //printf("! region is well-separated, computing force via multipole\n");
      rad = radius(myplanet.pos,regionToTest.com);
      F =  G * myplanet.m * regionToTest.mass * (regionToTest.com-myplanet.pos) / (rad*rad*rad);
      myplanet.acc = myplanet.acc + F/myplanet.m;
    }

    if (k==0) 
    {
      //space(lev); //printf("! region is not well-separated, descending to children\n");
      for(int l=0;l<8;l++) 
      {
	if(regions[regionToTest.kids[l]].numPln == 0)
	{
//          space(lev); printf("child %d is region %d, containing %d planets: empty\n",l,regionToTest.kids[l],regions[regionToTest.kids[l]].numPln);
	  continue;
	}
	else
        {
          //space(lev); //printf("! child %d is region %d, containing %d planets: recursing\n",l,regionToTest.kids[l],regions[regionToTest.kids[l]].numPln);
	  forceMagic(regions[regionToTest.kids[l]],myplanet,BD);
        }
      }
    }
  }
  //space(lev); //printf("! END   force on planet from region %d:%d,%d,%d containing %d planets : ",regionToTest.level,regionToTest.i,regionToTest.j,regionToTest.k,regionToTest.numPln); displist(regionToTest.planets);
}	
