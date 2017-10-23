/**********************************************
  The force is what gives a jedi his power.  It's an energy field created by all living things.  It surrounds us and penetrates us.  It binds the galaxy together.  Except that Meghan woman, she's a butt.
 **********************************************/

// Aaron is a poopy-head. --WF
//He's referring to you, Aaron Buttchins
//Shut up Megbutt
#include <omp.h>
#include <unistd.h>
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
extern double M_Sat;
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

double forceMagic(const region regionToTest, planet &myplanet, planet *BD, int fLoop_i){
	//assuming implementation in recover.c, N is number of planets
	int i,j,k,a,thread;
	double rad;
	vector F;
	pln *currPln;
	planet current;
	int lev=regionToTest.level;

	thread = omp_get_thread_num();
//	printf("!...... BEGIN FM lvl %d ..  i: %d .. thread %d\n", lev, myplanet.num, thread); 

	//space(lev); //printf("! BEGIN computing force on planet from region %d:%d,%d,%d containing %d planets : ",regionToTest.level,regionToTest.i,regionToTest.j,regionToTest.k,regionToTest.numPln); displist(regionToTest.planets);
	if(regionToTest.level == 4 || (regionToTest.numPln==1))
	{
		currPln=regionToTest.planets;
		//space(lev); //printf("! exactly, by iterating over planets.\n");
		while (currPln != NULL)
		{
			//space(lev); //printf("! Remaining list of planets in this region: "); displist(currPln);
//			printf("!...... FM WHILE lvl %d ..  i: %d .. thread %d .. plnNum %d\n", lev, myplanet.num, thread, currPln->plnNum);
			current=BD[currPln->plnNum];
			rad = radNew(myplanet.pos,current.pos);
			if(myplanet.num==0)
				printf("!..If..rad: %.16e\n", rad);
			if (rad > 0) // very crude check to make sure it's not the same planet pulling on itself
			{
				F = (current.pos - myplanet.pos) * G * myplanet.m * BD[regionToTest.planets->plnNum].m / (rad*rad*rad) ;//Direct force calculation
				myplanet.acc = myplanet.acc + F/myplanet.m;
				if(myplanet.acc.x > 10000 || myplanet.acc.y > 10000 || myplanet.acc.z > 10000 )
					printf("!...... ERROR FM lvl %d .. myplanet %d .acc = %f %f %f .. i: %d  current.num: %d\n", lev, myplanet.num, myplanet.acc.x, myplanet.acc.y, myplanet.acc.z, fLoop_i, current.num);
			
				if(myplanet.pos.x > 10000 || myplanet.pos.y > 10000 || myplanet.pos.z > 10000)
					printf("!...... ERROR FM lvl %d .. myplanet %d .pos = %f %f %f .. i: %d  current.num: %d\n", lev, myplanet.num, myplanet.pos.x, myplanet.pos.y, myplanet.pos.z, fLoop_i, current.num);
			}
			currPln=currPln->nextPln;
//			printf("!...... END FM Expl lvl %d .. i=%d .. thread %d\n", lev, myplanet.num, thread);
		}
	}

	else
	{
		//		intVec my = getIJKVec(n, myplanet.pos);  //intVec to store coordinates of region containing myplanet

		k = wellSepCoM(myplanet,regionToTest);

		if (k==1)
		{
			//space(lev); //printf("! region is well-separated, computing force via multipole\n");
			rad = radNew(myplanet.pos,regionToTest.com);
			if(myplanet.num==0)
				printf("!..Else..rad: %.16e\n", rad);
//			if(rad != rad)
//				printf("!...... ERROR rad = %f .... myplanet.pos = (%.3f, %.3f, %.3f) .. regionCoM = (%.3f, %.3f, %.3f)\n", rad, myplanet.pos.x, myplanet.pos.y, myplanet.pos.z, regionToTest.com.x, regionToTest.com.y, regionToTest.com.z);
			F =  G * myplanet.m * regionToTest.mass * (regionToTest.com-myplanet.pos) / (rad*rad*rad);
			myplanet.acc = myplanet.acc + F/myplanet.m;
			if(myplanet.acc.x > 10000 || myplanet.acc.y > 10000 || myplanet.acc.z > 10000 )
				printf("!...... ERROR FM lvl %d .. myplanet %d .acc = %f %f %f .. fLoop_i: %d  reg.com (%f %f %f)\n", lev, myplanet.num, myplanet.acc.x, myplanet.acc.y, myplanet.acc.z, fLoop_i, regionToTest.com.x, regionToTest.com.y, regionToTest.com.z);
			if(myplanet.pos.x > 10000 || myplanet.pos.y > 10000 || myplanet.pos.z > 10000)
				printf("!...... ERROR FM lvl %d .. myplanet %d .pos = %f %f %f .. fLoop_i: %d  reg.com (%f %f %f)\n", lev, myplanet.num, myplanet.pos.x, myplanet.pos.y, myplanet.pos.z, fLoop_i, regionToTest.com.x, regionToTest.com.y, regionToTest.com.z);

//			if(myplanet.acc.x != myplanet.acc.x || myplanet.acc.y != myplanet.acc.y || myplanet.acc.z != myplanet.acc.z)
//				printf("!...... ERROR FM lev %d ..  myplanet %d .acc = %f %f %f\n", lev, myplanet.num, myplanet.acc.x, myplanet.acc.y, myplanet.acc.z);
//			printf("!...... END FM Approx lvl %d .. i=%d .. thread %d\n", lev, myplanet.num, thread);
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
					forceMagic(regions[regionToTest.kids[l]],myplanet,BD,myplanet.num);
				}
			}
		}
	}
	//space(lev); //printf("! END   force on planet from region %d:%d,%d,%d containing %d planets : ",regionToTest.level,regionToTest.i,regionToTest.j,regionToTest.k,regionToTest.numPln); displist(regionToTest.planets);

	//Finally, do Saturn's contribution
//	vector zero; zero.x=0; zero.y=0; zero.z=0;
//	rad = radNew(myplanet.pos, zero);
//	myplanet.acc = myplanet.acc - G*M_Sat*myplanet.pos/(rad*rad*rad);
}	
