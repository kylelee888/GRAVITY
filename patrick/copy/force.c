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
/**********************************************
  Finds net force on planet "myplanet" from region "regionToTest"
 **********************************************/

void forceMagic(region regionToTest, planet &myplanet, planet *BD, double alpha){
	//assuming implementation in recover.c, N is number of planets
	int i,j,k,a;
	double rad;
	vector F;
	intVec my, other;
	pln *currPln;
	planet current;
	int lev=regionToTest.level;


	currPln=regionToTest.planets;
	current=BD[currPln->plnNum];
	//space(lev); //printf("! BEGIN computing force on planet from region %d:%d,%d,%d containing %d planets : ",regionToTest.level,regionToTest.i,regionToTest.j,regionToTest.k,regionToTest.numPln); displist(regionToTest.planets);
	if(regionToTest.level == (LVL-1) || (regionToTest.numPln==1))
	{
		//space(lev); //printf("! exactly, by iterating over planets.\n");
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

				F = (current.pos - myplanet.pos) * G * myplanet.m * BD[regionToTest.planets->plnNum].m / (rad*rad*rad) ;//Direct force calculation
				myplanet.acc = myplanet.acc + F/myplanet.m;
				if(F.x!=F.x || F.y!=F.y || F.z!=F.z)
				{
					printf("1--ERROR: Planet %d receiving force %e %e %e from planet %d located at %e %e %e\n",myplanet.num, F.x, F.y, F.z, current.num, current.pos.x, current.pos.y, current.pos.z);
				}
			}
			currPln=currPln->nextPln;
		}
	}

	else
	{
		//		intVec my = getIJKVec(n, myplanet.pos);  //intVec to store coordinates of region containing myplanet
		//k = 0;
		k = wellSepCoM(myplanet,regionToTest, alpha);

		if (k==1)
		{
			//space(lev); //printf("! region is well-separated, computing force via multipole\n");
			rad = radius(myplanet.pos,regionToTest.com);
			F =  G * myplanet.m * regionToTest.mass * (regionToTest.com-myplanet.pos) / (rad*rad*rad);
			if(F.x!=F.x || F.y!=F.y || F.z!=F.z)
			{
				printf("2--ERROR: Planet %d receiving force %e %e %e from planet %d located at %e %e %e\n",myplanet.num, F.x, F.y, F.z, current.num, current.pos.x, current.pos.y, current.pos.z);
			}
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
					forceMagic(regions[regionToTest.kids[l]],myplanet,BD,alpha);
				}
			}
		}
	}
	//space(lev); //printf("! END   force on planet from region %d:%d,%d,%d containing %d planets : ",regionToTest.level,regionToTest.i,regionToTest.j,regionToTest.k,regionToTest.numPln); displist(regionToTest.planets);
}	
