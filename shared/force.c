/**********************************************
force.c, used to find the force on each planet
by Aaron Hutchins
**********************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "list.h"
#include "FindRe[li]gion.h"

/**********************************************
Finds net force on planet "myplanet"
**********************************************/

/**********************************************
First use in code:  -> need to populate an array for regionsToTest, so start on level 1 to simplify

for(i=0,i>7,i++){
kids[i]=i+1;
}
forceMagic(BD[],kids[i],myplanet,1);
**********************************************/
void forceMagic(planet BD[], int regionsToTest[], int myplanet, int L){
	//assuming implementation in recover.c, N is number of planets
	int i,j,k,n;
	double rad;
	vector F;
	intVec my, other;
	int kids[8];
	
	int n = pow(2,L);
	
	if(L==5 || (/*there is only one planet in region*/)
	{
		//FIXME: can this be done with plnt_num on line 134 of recover.c?
		/* do force calculation */
	}
		
	else
	{
		intVec my = getIJKVec(n, BD[myplanet].level[L]);  //intVec to store coordinates of region containing myplanet
		
		for (a=0;a<8;a++)
		{ //cyles through each region in regionsToTest
			
			k = wellSepCoM(myplanet,regionsToTest[a]);
				
			if (k==1)
			{
				rad = radius(BD[myplanet].pos,//vector to regionsToTest[a]);	
				/* FIXME: Need to implement CoM code etc in a master struct before this can be effective */
			}
			if (k==0) 
			{
				for(l=0;l<8;l++) 
				{
					forceMagic(BD[],kids[l],myplanet,L-1);
					//kids[l] is from struct
				}
			}
		}	
	}		
}	
