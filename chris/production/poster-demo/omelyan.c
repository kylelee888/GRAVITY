#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "planet.h"

double eps = 0.1786178958448;
double lam = -0.2123418310626;
double chi = -0.06626458266982;

//Link to the paper I got this from
//https://arxiv.org/pdf/cond-mat/0110585.pdf

extern int N;
extern double dt;

double r_step1(planet BD[]){

  int i;
  for(i=0;i<N;i++){
    BD[i].pos = BD[i].pos + BD[i].vel*eps*dt;
  }
}

double v_step1(planet BD[]){
  int i;
  for(i=0;i<N;i++){
    BD[i].vel = BD[i].vel + BD[i].acc*(1-2*lam)*dt/2;
  }
}

double r_step2(planet BD[]){
  int i;
  for(i=0;i<N;i++){
    BD[i].pos = BD[i].pos + BD[i].vel*chi*dt;
  }
}

double v_step2(planet BD[]){
  int i;
  for(i=0;i<N;i++){
    BD[i].vel = BD[i].vel + BD[i].acc*lam*dt;
  }
}

double r_step3(planet BD[]){
  int i;
  for(i=0;i<N;i++){
    BD[i].pos = BD[i].pos + BD[i].vel*(1-2*(chi+eps))*dt;
  }
}

double v_step3(planet BD[]){
  int i;
  for(i=0;i<N;i++){
    BD[i].vel = BD[i].vel + BD[i].acc*lam*dt;
  }
}

double r_step4(planet BD[]){
  int i;
  for(i=0;i<N;i++){
    BD[i].pos = BD[i].pos + BD[i].vel*chi*dt;
  }
}

double v_stepf(planet BD[]){
  int i;
  for(i=0;i<N;i++){
    BD[i].vel = BD[i].vel + BD[i].acc*(1-2*lam)*dt/2;
  }
}

double r_stepf(planet BD[]){
  int i;
  for(i=0;i<N;i++){
    BD[i].pos = BD[i].pos + BD[i].vel*eps*dt;
  }
}
