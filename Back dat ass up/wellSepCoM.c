// code by Patrick, blame him for bugs -WF
#include <math.h>
#include "vector.h"
#include "planet.h"
#include "pln.h"
#include "region.h"

/*
Notes: could scale (aka multiply) alpha by some factor of how concentrated the center of mass is?
	-> Do this via comparing separations of CoM and the CoM's of each child region?
*/

int wellSepCoM(planet currPlanet, region tempReg)
{
        double beta,beta2,alpha2,length;
	int l = tempReg.level;
	vector tempCoM=tempReg.com;
	beta = length/pow(2,l);//length assumed to be global variable, length of one side of full space
		//beta is the width of a region for a given level
	beta2 = beta*beta;//Use beta^2 for cheaper math
	//alpha = 2; alpha2 = alpha*alpha;
	alpha2 = 4;//This factor determines "how many beta's" a region's center of mass needs to be from the planet we're concerned with in order to be considered well-separated. In this case, > 2 betas is well separated.

	vector sep = currPlanet.pos - tempCoM;
	if( sep*sep/beta2 > alpha2)
	{
		return(1);
	}
	else
	{
		return(0);
	}
}
