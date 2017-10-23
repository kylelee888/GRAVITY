#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "planet.h"

#define THETA 1.35120719195966

//link for Forest-Ruth Algorithm 4th order
//http://young.physics.ucsc.edu/115/leapfrog.pdf

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
