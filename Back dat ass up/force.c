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
extern double G;

/**********************************************
Finds net force on planet "myplanet" from region "regionToTest"
**********************************************/

void forceMagic(region regionToTest, planet myplanet){
	//assuming implementation in recover.c, N is number of planets
	int i,j,k,a;
	double rad;
	vector F;
	intVec my, other;
	
	if(regionToTest.level==5 || (regionToTest.numPln==1))
	{
		rad = radius(myplanet.pos,regionToTest.com);
		F = F + G * myplanet.m * regionToTest.mass * (regionToTest.com-myplanet.pos) / (rad*rad*rad);
	 	myplanet.acc = myplanet.acc + F/myplanet.m;
	}
		
	else
	{
//		intVec my = getIJKVec(n, myplanet.pos);  //intVec to store coordinates of region containing myplanet
		
		for (a=0;a<8;a++)
		{ //cyles through each region in regionsToTest
			
			k = wellSepCoM(myplanet,regionToTest);
				
			if (k==1)
			{
				rad = radius(myplanet.pos,regionToTest.com);
				F = F + G * myplanet.m * regionToTest.mass * (regionToTest.com-myplanet.pos) / (rad*rad*rad);
	 			myplanet.acc = myplanet.acc + F/myplanet.m;
			}
			if (k==0) 
			{
				for(int l=0;l<8;l++) 
				{
					forceMagic(regions[regionToTest.kids[l]],myplanet);
				}
			}
		}	
	}		
}	
