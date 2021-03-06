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

void recurse_divide_by_mass(region &r){
  int i;
  
  if( r.mass == 0 ) return;
  
  r.com = r.com / r.mass;
  if( r.level != LVL-1 ){
    for(i=0; i<8; i++){
      recurse_divide_by_mass(regions[r.kids[i]]);
    }
  }
}
 
//===============================
void populateRegions(planet BD[], int totalIndeces){
  
  int i,j;//here i will represent the level and j is planet number
  int a,b,c;//a b c are how many rows, columns, or slabs something is
  int n,k;
  int z,sum;
  int plnt_num;
  int index;
  double plnt_mass;
  vector plnt_COM;
  
  //sum now declared globally, called totalIndeces
  
  
  for( i=0; i<LVL; i++ ){
    //printf("!POPREG: working on level %d\n",i);
    n = pow(2,i);//n is number of boxes per side		
    for( j=0; j<N; j++ ){
      
      a = BD[j].pos.x*n/L;
      b = BD[j].pos.y*n/L;
      c = BD[j].pos.z*n/L;
      clamp(a,n);
      clamp(b,n);
      clamp(c,n);
      z = a + b*n + c*n*n;
      index = getArrayIndex(i,z);
      BD[j].level[i] = index;

      
      if(regions[index].numPln == 0){	
	regions[index].mass = 0;
	regions[index].com.x = 0;
	regions[index].com.y = 0;
	regions[index].com.z = 0;
      }
      
      regions[index].mass = regions[index].mass + BD[j].m;
      regions[index].com = regions[index].com + BD[j].m*BD[j].pos;
      regions[index].numPln++;
      
    }
  }
  recurse_divide_by_mass(regions[0]);
}
