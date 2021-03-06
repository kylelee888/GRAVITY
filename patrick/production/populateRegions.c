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

void clamp(int &a, int &b)
{
   if (a<0) a=0;
   if (a>=b) a=b-1;
}

//===============================
void populateRegions(planet BD[], int totalIndeces)
{
	//and put regions.MASS, regions.COM.x, regions.N and so forth

        int i,j;//here i will represent the level and j is planet number
	int a,b,c;//a b c are how many rows, columns, or slabs something is
	int n,k;
	int z,sum;
	int plnt_num;
	int index;
	double plnt_mass;
	vector plnt_COM;
 
	//sum now declared globally, called totalIndeces
	
 
	for( i=0; i<totalIndeces; i++ )//why was this only a problem for last level and not all??
	{
		regions[i].mass = 0;//need to assign 0 as it will be +=
		regions[i].com.x = 0;
		regions[i].com.y = 0;
		regions[i].com.z = 0;
		regions[i].numPln = 0;
	}
//	printf("!POPREG: finished initializing shit\n"); 
	for( i=0; i<LVL; i++ )
	{
	  //printf("!POPREG: working on level %d\n",i);
		n = pow(2,i);//n is number of boxes per side		
		for( j=0; j<N; j++ )
		{
			//EDIT 3/10: Only need to call arrayIndex (now getArrayIndex) once instead of 5 times. Call it and store in the int 'index'
 
			a = BD[j].pos.x*n/L;
			b = BD[j].pos.y*n/L;
			c = BD[j].pos.z*n/L;
			clamp(a,n);
			clamp(b,n);
			clamp(c,n);
			//split this into two functions, one that returns a,b,c and one that calls that to find the box number
			//also find MASS,COM,planet NUMBER in here too, way too expsneive the way i did it
			z = a + b*n + c*n*n;
			index = getArrayIndex(i,z);//Where i is the current level and z is the region in that level
			BD[j].level[i] = index;
			regions[index].mass = regions[index].mass + BD[j].m;
			regions[index].com = regions[index].com + BD[j].m*BD[j].pos;//not divided by total M yet
			regions[index].numPln++;//might be a simpler way
			//printf("!POPREG:   Planet %d (%.2f,%.2f,%.2f) is in region %d:%d,%d,%d/%d, which now has %d planets of mass %f\n",j,BD[j].pos.x,BD[j].pos.y,BD[j].pos.z,i,a,b,c,index,regions[index].numPln,regions[index].mass);
 
			//printf("Examining planet %d in level %d at location %.2e, %.2e, %.2e (L=%.2e): planet is in region (%d,%d,%d) = %d\n\n",j,i,BD[j].pos.x,BD[j].pos.y,BD[j].pos.z,L,a,b,c,BD[j].level[i]);
		}
	}
//	printf("!POPREG: Finished the second thing\n");

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
//	printf("!POPREG: Finished diving COM by total mass\n");

	//Check total number of planets in each region
	for( i=0; i<LVL; i++ )
	{
		n = pow(8,i);
		plnt_num = 0;
		for(j=0; j<n; j++)
		{
			plnt_num += regions[getArrayIndex(i,j)].numPln;
		}
//		printf("The total planets in level %d is %d with N = %d\n",i,plnt_num,N);
	}
//	printf("!POPREG: Finished czeching total # of planets in each region\n");
 
	//Check total mass in each region
	for( i=0; i<LVL; i++ ){
		n = pow(8,i);
		plnt_mass = 0;
		for(j=0; j<n; j++)
		{
			plnt_mass += regions[getArrayIndex(i,j)].mass;
		}
//		printf("The total mass in level %d is %f with M = %f\n",i,plnt_mass,mass(BD));//FIXME Could probably just calculate total mass once at beginning of progra
	} 	
//	printf("!POPREG: Finished chekking total mass in each level\n");
	
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
//		printf("The COM of boxes in level %d is (%f,%f,%f) with COM = (%f,%f,%f)\n",i,plnt_COM.x,plnt_COM.y,plnt_COM.z,regions[getArrayIndex(i,j)].com.x,regions[getArrayIndex(i,j)].com.y,regions[getArrayIndex(i,j)].com.z);
	}
//	printf("!POPREG: Finished checking COM of each level\n");
}
