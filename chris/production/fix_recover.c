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
#include "resetregion.h"
#include "parentChild.h"
#include "omelyan.h"
#include "timing.h"

//UNITS
//Mass   : Earth Mass
//Length : AU (astronomical unit)
//Time   : YEAR
//fix planet sizes, way too big

int LVL = 6; 
int numDel;
region *regions;
int totalIndeces;
double alpha;
int collision = 0;
int collision_check;
int collision_number;
int collision_pair_1[30000];
int collision_pair_2[30000];

//GLOBAL variable declaration
int N; //number of planets    
double dt; //time step
double G  = 4*M_PI*M_PI; //Newton's constant 
int a; //number of collisions
double L; //inital grid size
vector com; //center of mass to test my functions


//=======================
void reset_acc( planet BD[] ){
  int i;
  for(i=0;i<N;i++){
    BD[i].acc.x = 0;
    BD[i].acc.y = 0;
    BD[i].acc.z = 0;
  }
}

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
vector Force( planet BD[] ){

  int i,j;
  double rad;
  vector F;

  for(i=0; i<N; i++){
    BD[i].acc = 0;
  }
  
  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++){

      rad = radius(BD[i].pos,BD[j].pos);
      F   = G * BD[j].m * BD[i].m * (BD[j].pos-BD[i].pos) / (rad*rad*rad);
      BD[i].acc = BD[i].acc + F/BD[i].m;
      BD[j].acc = BD[j].acc - F/BD[j].m;
    }
  }
}
//=================================

void compactify(planet BD[]){
  printf("!I am trying to compactify but why?\n");
  for(int i=0; i<N; i++){
    if (BD[i].m==0){
      N--;
      BD[i].pos = BD[N].pos;
      BD[i].vel = BD[N].vel;
      BD[i].m   = BD[N].m;
      BD[i].r   = BD[N].r;
    }
  }
}

//===================================
double mag_sq(vector vec){
  return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}
//===================================
double mag(vector vec){
  return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}
//===================================

int rElemNew( planet BD[] ){

  int i,j;
  int ndel=0;
  double rho = 20;
  double rad;

  double bf = 0;
  double mi,mj;
  vector vcm, v1, v2;
  double v1_sq, v2_sq;
  double theta, phi;
  //printf("!%d\n",collision_number);
  for(int k=0; k<collision_number; k++){
    i = collision_pair_1[k];
    j = collision_pair_2[k];
    if(i != j){
      rad = radius(BD[i].pos, BD[j].pos);
      if(rad < (BD[i].r + BD[j].r) ){
	printf("!N=%d\n",N);
	if( rad >= radius(BD[i].pos+BD[i].vel*dt, BD[j].pos+BD[j].vel*dt) ){//meaning they are moving towards each other
	  mi = BD[i].m;
	  mj = BD[j].m;
	  
	  vcm = (BD[i].vel*mi + BD[j].vel*mj)/(mj+mi);
	  
	  //v1 = BD[i].vel - vcm;
	  //v2 = BD[j].vel - vcm;
	  // printf("!%f %f\n",v1.x,v2.x);
	  //v1_sq = mag_sq(v1);
	  // v2_sq = mag_sq(v2);
	  //printf("!%f %f\n",v1_sq, v2_sq);
	  //now we are in reference frame with zero momentum
	  /*
	  if( 0 ){
	    //printf("!They stick\n");
	    BD[i].pos = (BD[i].m*BD[i].pos + BD[j].m*BD[j].pos)/(BD[i].m + BD[j].m); 
	    BD[i].vel = (BD[i].m*BD[i].vel + BD[j].m*BD[j].vel)/(BD[i].m + BD[j].m);
	    BD[i].m   = (BD[i].m + BD[j].m);
	    BD[i].r   = pow( BD[i].m/rho, 1.0/3 );
	    //don't shift planets yet, since if planet N-1 also collides this breaks things. 
	    //just mark planet j as "to be deleted" by setting its mass to zero.
	    //we can then compactify the planet list later.
	    BD[j].m=0;
	    ndel++;
	  }
	  */
	  //else{
	  printf("!Things are colliding\n");
	  theta = acos((BD[i].vel*vcm)/(mag(BD[i].vel)+mag(vcm)));
	    printf("!%f\n",theta);

	    BD[i].vel.x = BD[i].vel.x*(cos(theta)*cos(theta) - sin(theta)*sin(theta)) - sin(2*theta)*BD[i].vel.y;
	    BD[i].vel.y = -sin(2*theta)*BD[i].vel.x + BD[i].vel.y*(sin(theta)*sin(theta) - cos(theta)*cos(theta));

	    theta = acos((BD[j].vel*vcm)/(mag(BD[j].vel)+mag(vcm)));
	    printf("%f\n",phi);
	    BD[j].vel.x = BD[j].vel.x*(cos(theta)*cos(theta) - sin(theta)*sin(theta)) - sin(2*theta)*BD[j].vel.y;
	    BD[j].vel.y = -sin(2*theta)*BD[j].vel.x + BD[j].vel.y*(sin(theta)*sin(theta) - cos(theta)*cos(theta));
	    
	}
      }		
      //else{/
      //	  //printf("!They just bounced yo, give them a chance to get away\n");
      //	}
    }
  }
  //printf("!ndel = %d\n",ndel);
  if(ndel) compactify(BD);
  return ndel;
}

//==============================
void omelyan( planet BD[], int method ){

  int i;

  //r_step1(BD);

  //reset_acc(BD);
  
  Force(BD);
  for(i=0; i<N; i++){
    BD[i].vel = BD[i].vel + BD[i].acc*dt;
    BD[i].pos = BD[i].pos + BD[i].vel*dt;
  }
  
  //v_step1(BD);
  //r_step2(BD);

  //reset_acc(BD);
  //Force(BD);
  
  //v_step2(BD);
  //r_step3(BD);

  //reset_acc(BD);
  //Force(BD);
  
  //v_step3(BD);
  //r_step4(BD);

  //reset_acc(BD);
  //Force(BD);
  
  //v_stepf(BD);
  //r_stepf(BD);
  
}
//========================================
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

  if( nParam != 8 ) {printf("\n!Usage: ./<this> <N> <dt> <frameskip> <initial velocity> <box width> <initial conditions width> <alpha>\n\n"); exit(0);}

  N          = atoi(paramList[1]);
  dt         = atof(paramList[2]);
  int fSkip  = atoi(paramList[3]);
  double v   = atof(paramList[4]);
  L          = atof(paramList[5]);
  double scatter  = atof(paramList[6]);
  alpha      = atof(paramList[7]);

  planet BD[N];
  vector vcm;
  int i;
  totalIndeces = 0;
  int frame = 0;
  double MASS = 0;
  double Ei,Pi;
  
  double force_t;
  vector F_exact[N];

  for(i=0; i<LVL; i++)
  {
   totalIndeces += pow(8,i);
  }
  regions=(region*)malloc(sizeof(region)*totalIndeces);
      
  findKids(LVL);

  vcm.x = 0;
  vcm.y = 0;
  vcm.z = 0;
  
  double t,theta,rho=20;
  com.x = 0;
  com.y = 0;
  com.z = 0;
  int q;
  int count = 0;
      
  for( i=0; i<N; i++){
    
    BD[i].num = i;
    
    BD[i].pos.x = L * 0.5 + (drand48()-0.5)*scatter; //not neg so fits planetRegion
    BD[i].pos.y = L * 0.5 + (drand48()-0.5)*scatter; //not neg so fits planetRegion
    //BD[i].pos.z = 0;
    BD[i].pos.z = L * 0.5 + (drand48()-0.5)*scatter; //not neg so fits planetRegion
    
    BD[i].vel.x = v*(drand48()-.5);
    BD[i].vel.y = v*(drand48()-.5);
    //BD[i].vel.z = 0;
    BD[i].vel.z = v*(drand48()-.5);
    
    BD[i].r   = 0.3;
    BD[i].m  = rho*BD[i].r*BD[i].r*BD[i].r;
    
    com = com + BD[i].pos * BD[i].m;
  }
  
  for( i=0; i<N; i++ ){
    MASS += BD[i].m;
  }
    
  for( i=0; i<N; i++ ){
    vcm = vcm + BD[i].m*BD[i].vel;
  }
  
  com = com/MASS;
  vcm = vcm/MASS;
  
  for( i=0; i<N; i++ ){
    BD[i].vel = BD[i].vel - vcm;
  }
  
  Ei = energy(BD);
  Pi = mom(BD);
  
  for( t=0; 1; t+=dt ){
    
    //reset_region(regions[0]);
        
    //populateRegions(BD, totalIndeces);//find region properties
    
    //popLists(BD);//populate lists

    //euler(BD);
    omelyan(BD,0);
  
    //a+=rElemNew(BD);
    
    if( frame % fSkip == 0 ){
    
      printf("T -0.8 0.82\nE = %f : P = %f : # Collisions = %d : t = %f\n",energy(BD),mom(BD),a,t);
      for( i=0; i<N; i++){
	printf("c3 %e %e %e %e\n",BD[i].pos.x,BD[i].pos.y,BD[i].pos.z,BD[i].r);
      }
      printf("F\n"); 
    }
    frame++;
  }
}




