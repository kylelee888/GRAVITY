#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "planet.h"
#include "pln.h"
#include "region.h"
#include "force.h"
#include "populateRegions.h"
#include "planetListFunctions.h"
#include "LukeIAmYourParent.h"

#define THETA 1.35120719195966
//UNITS
//Mass   : Earth Mass
//Length : AU (astronomical unit)
//Time   : YEAR
//fix planet sizes, way too big

// link for Forest-Ruth Algorithm 4th order
//http://young.physics.ucsc.edu/115/leapfrog.pdf

int LVL = 5; 

region *regions;

//GLOBAL variable declaration
int N; //number of planets    
double dt; //time step
double G  = 4*M_PI*M_PI; //newtons constant 
int a; //number of collisions
double L; //inital grid size
vector com; //center of mass to test my functions
//double THETA = 1/(2-pow(2,1./3));


//=======================
int arrayIndex( int level, int box ){

  //this will take the planet box at a level
  //and return where to put it in the region =array
  
  int i;
  int TOT = 0;
  for( i=0; i<level; i++ ){//because only want the boxes in previous levels
    TOT += pow(8,i);
  }
  return TOT + box;
}

//=======================
vector posHalfStep( planet BD[] ){
  
  int i;
  
  for( i=0; i<N; i++ ){
    BD[i].pos = BD[i].pos + BD[i].vel*dt/2;
  }
}

//=======================
vector velFullStep( planet BD[] ){
  
  int i;
  
  for( i=0; i<N; i++ ){
    BD[i].vel = BD[i].vel + BD[i].acc*dt;
  }
}

//=======================
vector pos4_FWD( planet BD[] ){
  
  int i;

  for( i=0; i<N; i++ ){
    BD[i].pos = BD[i].pos + THETA*BD[i].vel*dt/2;
  }
}

//=======================
vector pos4_BCK( planet BD[] ){

  int i;

  for( i=0; i<N; i++ ){
    BD[i].pos = BD[i].pos + (1-THETA)*BD[i].vel*dt/2;
  }
}

//=======================
vector vel4_FWD( planet BD[] ){

  int i;

  for( i=0; i<N; i++ ){
    BD[i].vel = BD[i].vel + THETA*BD[i].acc*dt;
  }
}

//=======================
vector vel4_BCK( planet BD[] ){

  int i;

  for( i=0; i<N; i++ ){
    BD[i].vel = BD[i].vel + (1-2*THETA)*BD[i].acc*dt;
  }
}

//=======================
vector Force( planet BD[] ){

  int i,j;
  double rad;
  vector F;
  
  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++){//the i+1 is what makes it faster instead of j=0 and going through to N

      rad = radius(BD[i].pos,BD[j].pos);
      F   = G * BD[j].m * BD[i].m * (BD[j].pos-BD[i].pos) / (rad*rad*rad);
      BD[i].acc = BD[i].acc + F/BD[i].m;
      BD[j].acc = BD[j].acc - F/BD[j].m;
    }
  }
}
//=======================

double rElem( planet BD[] ){

  int i,j;
  double rad;
  double rho = 20;

  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++ ){
      rad = radius(BD[i].pos,BD[j].pos);
            
      if( rad < 0.8*(BD[i].r+BD[j].r) ){
      
	BD[i].pos = (BD[i].m*BD[i].pos + BD[j].m*BD[j].pos)/(BD[i].m + BD[j].m);//find new properties 
	BD[i].vel = (BD[i].m*BD[i].vel + BD[j].m*BD[j].vel)/(BD[i].m + BD[j].m);
	BD[i].m   = (BD[i].m + BD[j].m);
	BD[i].r   = pow( BD[i].m/rho, 1.0/3 );
	
	BD[j].pos = BD[N-1].pos;//move information
	BD[j].vel = BD[N-1].vel;
	BD[j].m   = BD[N-1].m;
	BD[j].r   = BD[N-1].r;
	  
	N--;
	a++;
      }
    }
  }
}
  
//=======================

double leap( planet BD[] ){
  printf("LEAP: Start\n"); 
  vector a[N];
  int i;

  for( i=0; i<N; i++ ){
    BD[i].acc.x = 0;
    BD[i].acc.y = 0;
    BD[i].acc.z = 0;
  }

  posHalfStep(BD);

  //Force(BD,a); old N^2 calc
  for( i=0; i<N; i++ ){
    //just need to pass planet and region 0
    forceMagic( regions[0], BD[i] );
  }
  //FIXME: change planet struct to have .acc
  velFullStep(BD);
    
  posHalfStep(BD);
}  


//===================

double energy( planet BD[] ){
  double KE = 0,PE = 0;
  int i,j;
  
  for( i=0; i<N; i++ ){
    KE += 0.5*BD[i].m*(BD[i].vel*BD[i].vel);
  }

  for( i=0; i<N; i++){
    for( j=i+1; j<N; j++ ){
      PE += BD[i].m*BD[j].m / radius(BD[i].pos,BD[j].pos);
    }
  }
  return KE - G*PE;
}

//==================
double mom( planet BD[] ){

  int i;
  vector p;
  p.x = 0;
  p.y = 0;
  p.z = 0;
  
  for( i=0; i<N; i++ ){
    p = p + BD[i].m*BD[i].vel;
  }

  return sqrt(p*p);
}

//==================

int main(int nParam, char **paramList){

  if( nParam < 6 ) {printf("\n!Usage: ./string <N> <dt> <frameskip> <initial velocity> <volume> \n\n"); exit(0);}

  N          = atoi(paramList[1]);
  dt         = atof(paramList[2]);
  int fSkip  = atoi(paramList[3]);
  double v   = atof(paramList[4]);
  L          = atof(paramList[5]);

  planet BD[N];
  vector vcm;
  int i;
  int totalIndeces;
  int frame = 0;
  double MASS = 0;

  for(i=0; i<LVL; i++)
  {
   totalIndeces += pow(8,i);
  }
  regions=(region*)malloc(sizeof(region)*totalIndeces);//FIXME Is this correct?
  printf("!Total # of indeces in region array: %d\n", totalIndeces);
    
  Children(LVL);

  vcm.x = 0;
  vcm.y = 0;
  vcm.z = 0;
  
  double t,theta,rho=20;
  com.x = 0;//to check the COM finding of my function
  com.y = 0;
  com.z = 0;
  
  

  for( i=0; i<N; i++){

    BD[i].pos.x = L*drand48();//not neg so fits planetRegion
    BD[i].pos.y = L*drand48();
    BD[i].pos.z = L*drand48();

    //printf("x = %.2e, y = %.2e, z = %.2e\n",BD[i].pos.x, BD[i].pos.y, BD[i].pos.z);

    BD[i].vel.x = v*(drand48()-.5);
    BD[i].vel.y = v*(drand48()-.5);
    BD[i].vel.z = v*(drand48()-.5);

    BD[i].r   = 0.1;
    BD[i].m  = rho*BD[i].r*BD[i].r*BD[i].r;

    com = com + BD[i].pos * BD[i].m;
  }

  for( i=0; i<N; i++ ){
    MASS += BD[i].m;
  }
  
  printf("!Finished initializing planet locations\n");

  for( i=0; i<N; i++ ){
    vcm = vcm + BD[i].m*BD[i].vel;
  }
  
 com = com/MASS;
 printf("!The center of mass is (%f,%f,%f)\n",com.x,com.y,com.z);

 vcm = vcm/MASS;
  
  for( i=0; i<N; i++ ){
    BD[i].vel = BD[i].vel - vcm;
  }

  //where the magic happens
  
  for( t=0; 1; t+=dt ){
    printf("!Time: t = %f\n", t);
    populateRegions(BD, totalIndeces);//find region properties
    printf("!called popRegions \n");
    for( i=0; i<totalIndeces; i++ ){//clear lists
      clearList(&(regions[i].planets));//if we call the pointer .list
    }
    printf("!Called clearList\n");
    popLists(BD);//populate lists
    printf("!Called popLists\n");
    leap(BD);
    printf("!Called leap\n");
    rElem(BD);
    printf("!Called rElem\n");

    if( frame % fSkip == 0 ){
      printf("T -0.8 0.82\nE = %f : P = %f : # Collisions = %d\n",energy(BD),mom(BD),a);
      for( i=0; i<N; i++){
	printf("c3 %e %e %e %e\n",BD[i].pos.x,BD[i].pos.y,BD[i].pos.z,BD[i].r);
      }
      printf("F\n");
    }
    frame++;
  }
}

