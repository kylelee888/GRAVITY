#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "planet.h"
#include "pln.h"
#include "region.h"

extern planet BD[];

double * quadrupole(region jim){

	int i,j,k;
	double ru;
	double rv;
	double rcu;
	double rcv;
	vector rc = jim.com; 
	double earth.m;
	planet earth;
	double Q[9];
	double r2;
											
	for(i=0; i<9; i++)
	{
		Q[i] = 0;
	}
	
	r2 = earth.pos.x*earth.pos.x + earth.pos.y*earth.pos.y + earth.pos.z*earth.pos.z;
//you have to take the sqrt you dumbass

	for(i=0;i<jim.numPln;i++)
	{
		//printf("%d\n",jim.planets[i].plnNum);
		earth = BD[jim.planets[i].plnNum];

    		for(j=0;j<3;j++)
		{
			if(j==0)
			{
				ru = earth.pos.x;
				rcu = rc.x;
			}
			if(j==1)
			{
				ru = earth.pos.y;
				rcu = rc.y;
			}
			if(j==2)
			{
				ru = earth.pos.z;
				rcu = rc.z;
			}
			for(k=0;k<3;k++)
			{
				if(k==0)
				{
					rv = earth.pos.x;
					rcv = rc.x;
				}
				if(k==1)
				{
					rv = earth.pos.y;
					rcv = rc.y;
				}
				if(k==2)
				{
					rv = earth.pos.z;
					rcv = rc.z;
				}

				if(j==k)
				{
					Q[3*j+k] += (3*(ru-rcu)*(rv-rcv) - r2)*earth.m;
				}
				else
				{
					Q[3*j+k] += 3*(ru-rcu)*(rv-rcv)*earth.m;
				}

			}
		}
	}


}
