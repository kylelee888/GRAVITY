#include <stdio.h>
#include <math.h>
#include "planet.h"
#include "vector.h"

double bf = 5*5; //bounce factor

double mag_sq(vector vec){

  return vec.x*vex.x + vec.y*vec.y + vec.z*vec.z;
  
}
void collision(planet BD[], int i, int j){

  double mi,mj;
  vector vcm, v1, v2, v1_sq, v2_sq;
//getting passed the ideces of planets colliding
  mi = BD[i].m;
  mj = BD[j].m;
  
  vcm = (BD[i].vel*mi + BD[j].vel*mj)/(mj+mi);

  v1 = BD[i].vel - vcm;
  v2 = BD[j].vel - vcm;

  v1_sq = mag_sq(v1);
  v2_sq = mag_sq(v2);
  //now we are in reference frame with zero momentum
  if( v1_sq/v2_sq < bf || v1_sq/v2_sq < 1/bf ){

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

  else{
    v1 = -v1;
    v2 = -v2;
  }

  BD[i].vel = v1 + vcm;
  BD[j].vel = v2 + vcm;
}
