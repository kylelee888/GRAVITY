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

int getregion(vector p, int lev)
{
  int n=pow(2,lev);
  int a = p.x*n/L;
  int b = p.y*n/L;
  int c = p.z*n/L;
  clamp(a,n);
  clamp(b,n);
  clamp(c,n);
  int z = a + b*n + c*n*n;
  return getArrayIndex(lev, z);
}

void recurse_divide_by_mass(region &r)
{
  int i;
  if (r.mass == 0) return;
  r.com = r.com / r.mass;
  if (r.level != LVL-1)
  {
    for (i=0; i<8; i++)
    {
      recurse_divide_by_mass(regions[r.kids[i]]);
    }
  }
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
	
 
//	for( i=0; i<totalIndeces; i++ )//why was this only a problem for last level and not all??
//	{
//		regions[i].mass = 0;//need to assign 0 as it will be +=
//		regions[i].com.x = 0;
//		regions[i].com.y = 0;
//		regions[i].com.z = 0;
//		regions[i].numPln = 0;
//	}
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
                        if (regions[index].numPln == 0) // this is the first time we've seen it, so we need to clear
			  regions[index].mass = regions[index].com.x = regions[index].com.y = regions[index].com.z = 0;
			
                        regions[index].mass = regions[index].mass + BD[j].m;
			regions[index].com = regions[index].com + BD[j].m*BD[j].pos;
			regions[index].numPln++;//might be a simpler way
//			printf("!POPREG:   Planet %d (%.2f,%.2f,%.2f) is in region %d:%d,%d,%d/%d, which now has %d planets of mass %f\n",j,BD[j].pos.x,BD[j].pos.y,BD[j].pos.z,i,a,b,c,index,regions[index].numPln,regions[index].mass);
 
			//printf("Examining planet %d in level %d at location %.2e, %.2e, %.2e (L=%.2e): planet is in region (%d,%d,%d) = %d\n\n",j,i,BD[j].pos.x,BD[j].pos.y,BD[j].pos.z,L,a,b,c,BD[j].level[i]);
		}
	}
        recurse_divide_by_mass(regions[0]);
}
