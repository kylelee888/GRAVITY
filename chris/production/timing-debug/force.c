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
extern int collision_pair_1[];
extern int collision_pair_2[];
extern int collision_number;
extern int collision_check;
extern int LVL;
extern int N;
extern int exact, approx;

/**********************************************
  Finds net force on planet "myplanet" from region "regionToTest"
 **********************************************/

void force_naive(planet &myplanet, planet BD[])
{
  double rad;
  vector F;
  planet current;
  for (int i=0; i<N; i++)
  {
    current=BD[i];
    rad = radius(myplanet.pos,current.pos);
    if (rad > 0)
    {
      exact++;
      F = (current.pos - myplanet.pos) * G * myplanet.m * current.m / (rad*rad*rad) ;//Direct force calculation
      myplanet.acc = myplanet.acc + F/myplanet.m;
      
      if( collision_check ){
	if( 2*rad < (myplanet.r + current.r) && (myplanet.num > current.num )){//this actually decides if things hit.
	  collision_pair_1[collision_number] = myplanet.num;
	  collision_pair_2[collision_number] = current.num;
	  collision_number++;

	}
      }
    }
  }
}

void forceMagic(region regionToTest, planet &myplanet, planet *BD, double alpha){
  if (alpha == 0)
  {
    force_naive(myplanet, BD);
    return;
  }
  //assuming implementation in recover.c, N is number of planets
  int i,j,k,a;
  double rad;
  vector F;
  intVec my, other;
  pln *currPln;
  planet current;
  int lev=regionToTest.level;
  

  
  if(regionToTest.level == (LVL-1) || (regionToTest.numPln==1))
  {
    currPln=regionToTest.planets;
    while (currPln != NULL)
    {
      //space(lev); //printf("! Remaining list of planets in this region: "); displist(currPln);
      current=BD[currPln->plnNum];
      rad = radius(myplanet.pos,current.pos);
      
	
      if (rad > 0){ // very crude check to make sure it's not the same planet pulling on itself

	if( collision_check ){
	  if( 2*rad < (myplanet.r + current.r) && (myplanet.num > current.num )){//this actually decides if things hit.
	    collision_pair_1[collision_number] = myplanet.num;
	    collision_pair_2[collision_number] = current.num;
	    collision_number++;

	  }
	}
	    
        F = (current.pos - myplanet.pos) * G * myplanet.m * BD[currPln->plnNum].m / (rad*rad*rad) ;//Direct force calculation
        exact++;
        myplanet.acc = myplanet.acc + F/myplanet.m;
      }
      currPln=currPln->nextPln;
   }
  }

  else
  {
    //		intVec my = getIJKVec(n, myplanet.pos);  //intVec to store coordinates of region containing myplanet
    //k = 0;
    k = wellSepCoM(myplanet,regionToTest);
    if(k > (alpha*alpha) ){
      //space(lev); //printf("! region is well-separated, computing force via multipole\n");
      rad = radius(myplanet.pos,regionToTest.com);
      // F =  G * myplanet.m * regionToTest.mass * (regionToTest.com-myplanet.pos) / (rad*rad*rad);
      myplanet.acc = myplanet.acc + regionToTest.mass * (regionToTest.com-myplanet.pos) / (rad*rad*rad);
      approx++;
      //myplanet.acc = myplanet.acc + F/myplanet.m;
    }

    else
    {
      //space(lev); //printf("! region is not well-separated, descending to children\n");
      for(int l=0;l<8;l++) 
      {
	if(regions[regionToTest.kids[l]].numPln == 0)
	{
	  continue;
	}
	else
        {
	  forceMagic(regions[regionToTest.kids[l]],myplanet,BD,alpha);
        }
      }
    }
  }
}	
