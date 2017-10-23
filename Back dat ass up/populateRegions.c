#include <stdio.h>
#include <math.h>
#include "getArrayIndex.h"
#include "vector.h"
#include "pln.h"
#include "region.h"
#include "planet.h"

extern int LVL,N;
extern double L;

double mass( planet BD[] ){
  
  int i;
  double mass=0.0;
  for( i=0; i<N; i++ ){
    mass += BD[i].m;
  }
  return mass;
}

//===============================
void populateRegions(planet BD[])//FIXME should this be void since it isn't returning anything?
{
	//and put regions.MASS, regions.COM.x, regions.N and so forth

        int i,j,totalIndeces;//here i will represent the level and j is planet number
	int a,b,c;//a b c are how many rows, columns, or slabs something is
	int n,k;
	int z,sum;
	int plnt_num;
	int index;
	double plnt_mass;
	vector plnt_COM;
 
	//sum now declared globally, called totalIndeces
	int NMBR[sum];//FIXME Not sure this is necessary b/c of how linked lists and addPlanet works
 
	for( i=0; i<totalIndeces; i++ )//why was this only a problem for last level and not all??
	{
		regions[i].mass = 0;//need to assign 0 as it will be +=
		regions[i].com.x = 0;
		regions[i].com.y = 0;
		regions[i].com.z = 0;
		NMBR[i]	= 0;
	}
	 
	for( i=0; i<LVL; i++ )
	{
		n = pow(2,i);//n is number of boxes per side		
		for( j=0; j<N; j++ )
		{
			//EDIT 3/10: Only need to call arrayIndex (now getArrayIndex) once instead of 5 times. Call it and store in the int 'index'
 
			a = BD[j].pos.x*n/L;
			b = BD[j].pos.y*n/L;
			c = BD[j].pos.z*n/L;
			// FIXME: MASS, COM, etc., are pieces of a larger array
			//split this into two functions, one that returns a,b,c and one that calls that to find the box number
			//also find MASS,COM,planet NUMBER in here too, way too expsneive the way i did it
			z = a + b*n + c*n*n;
			BD[j].level[i] = z;
			index = getArrayIndex(i,z);//Where i is the current level and z is the region in that level
			regions[index].mass = regions[index].mass + BD[j].m;
			regions[index].com = regions[index].com + BD[j].m*BD[j].pos;//not divided by total M yet
			regions[index].numPln++;//might be a simpler way
 
			//printf("Examining planet %d in level %d at location %.2e, %.2e, %.2e (L=%.2e): planet is in region (%d,%d,%d) = %d\n\n",j,i,BD[j].pos.x,BD[j].pos.y,BD[j].pos.z,L,a,b,c,BD[j].level[i]);
		}
	}
 
	//Divide COM by total mass to get actual COM
	for( i=0; i<totalIndeces; i++ )
	{
		if(regions[i].mass == 0)
		{
			//skip cus you cant divide by zero and COM will be zero anyway
		}
		else
		{
			regions[i].com = regions[i].com / regions[i].mass;
			//printf("The com of box %d = (%f,%f,%f)\n",i,COM[i].x,COM[i].y,COM[i].z);
		}
	}
 
	//Check total number of planets in each region
	for( i=0; i<LVL; i++ )
	{
		n = pow(8,i);
		plnt_num = 0;
		for(j=0; j<n; j++)
		{
			plnt_num += regions[getArrayIndex(i,j)].numPln;
		}
		printf("The total planets inside region	%d of level %d is %d with N = %d\n",j,i,plnt_num,N);
	}
 
	//Check total mass in each region
	for( i=0; i<LVL; i++ ){
		n = pow(8,i);
		plnt_mass = 0;
		for(j=0; j<n; j++)
		{
			plnt_mass += regions[getArrayIndex(i,j)].mass;
		}
		printf("The total mass inside boxes in level %d is %f with M = %f\n",i,plnt_mass,mass(BD));//FIXME Could probably just calculate total mass once at beginning of progra
	}
 
	//Check COM of each region	FIXME Shouldn't we just check the location of the center of mass?
	for( i=0; i<LVL; i++ ){
		n = pow(8,i);
		plnt_COM.x = 0;//reset the COM before each level
		plnt_COM.y = 0;
		plnt_COM.z = 0;
		for( j=0; j<n; j++ ){
		  plnt_COM = plnt_COM + regions[getArrayIndex(i,j)].com*regions[getArrayIndex(i,j)].mass;//COM*MASS is the pos*mass
		}
		plnt_COM = plnt_COM / mass(BD); //can divide by mass(BD) here as it was proved the mass is good
		printf("The COM of boxes in level %d is (%f,%f,%f) with COM = (%f,%f,%f)\n",i,plnt_COM.x,plnt_COM.y,plnt_COM.z,regions[getArrayIndex(i,j)].com.x,regions[getArrayIndex(i,j)].com.y,regions[getArrayIndex(i,j)].com.z);
	}
}
