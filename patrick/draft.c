#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "planet.h"


//so the problem occurs after the first level loops through!!
//UNITS
//Mass   : Earth Mass
//Length : AU (astronomical unit)
//Time   : YEAR

int LVL = 5; //number of levels?
int N;      //do a planet struct
double dt;
double G  = 4*M_PI*M_PI;
int a;
double L;

//=======================
planet planetRegion( planet BD[] ){//believe it works

  int i,j;//here i will represent the level and j is planet number
  int a,b,c;//a b c are how many rows, columns, or slabs something is
  int n; 

  printf("\n\nYou are now looking at the function planetRegion\n\n");
 
  for( i=0; i<LVL; i++ ){
    
    n = pow(2,i);//n is number of boxes per side    
    
    for( j=0; j<N; j++ ){

      a = BD[j].pos.x*n/L; 
      b = BD[j].pos.y*n/L; 
      c = BD[j].pos.z*n/L; 

      BD[j].level[i] = a + n*b + n*n*c;
      
      //printf("Examining planet %d in level %d at location %.2e, %.2e, %.2e (L=%.2e): planet is in region (%d,%d,%d) = %d\n\n",j,i,BD[j].pos.x,BD[j].pos.y,BD[j].pos.z,L,a,b,c,BD[j].level[i]);
    }
  }
  //for( j=0; j<N; j++ ){
  //  for( i=0; i<LVL; i++ ){
  //    printf("Planet %d: Level %d = %d ",j,i,BD[j].level[i]);
  //  }
  //  printf("\n");
  //}
}

//=======================
planet regionCOM( planet BD[] ){

  int i,j;

  for( i=0; i<LVL; i++ ){

    for( j=0; j<N; j++) {
      //some kind of if statement
      

  }

}
//=======================
double radius( vector BD1, vector BD2 ){//vector r1, vector r2){
  return sqrt( (BD1 - BD2) * (BD1 - BD2) );
  //return sqrt( (r1-r2) * (r1-r2) );
}

//=======================
vector posHalfStep( planet BD[] ){//vector pos[], vector vel[] ){

  int i;
  
  for( i=0; i<N; i++ ){
    BD[i].pos = BD[i].pos + BD[i].vel*dt/2;
    //pos[i] = pos[i] + vel[i]*(dt/2);
  }
}

//=======================
vector velFullStep( vector a[], planet BD[] ){//vector vel[] ){

  int i;
  for( i=0; i<N; i++ ){
    BD[i].vel = BD[i].vel + a[i]*dt;
    //vel[i] = vel[i] + a[i]*dt;
  }
}

//=======================
vector Force( planet BD[], vector a[] ){//vector pos[], vector F, vector a[], double m[] ){

  int i,j;
  double rad;
  vector F;
  
  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++){//the i+1 is what makes it faster instead of j=0 and going through to N

      rad = radius(BD[i].pos,BD[j].pos);
      F   = G * BD[j].m * BD[i].m * (BD[j].pos-BD[i].pos) / (rad*rad*rad);
      a[i] = a[i] + F/BD[i].m;
      a[j] = a[j] - F/BD[j].m;
      //rad  = radius(pos[i],pos[j]);
      //F    = G * m[j] * m[i] * (pos[j]-pos[i]) / (rad*rad*rad );
      //a[i] = a[i] + F/m[i];
      //a[j] = a[j] - F/m[j];
    }
  }
}
//=======================

double rElem( planet BD[] ){//vector pos[], vector vel[], double m[],  double r[] ){//int i, int j){

  int i,j;
  double rad;
  double rho = 20;

  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++ ){
      rad = radius(BD[i].pos,BD[j].pos);
      //rad = radius(pos[i],pos[j]);
      
      if( rad < 0.8*(BD[i].r+BD[j].r) ){
      //if( rad < 0.8*(r[i]+r[j]) ){

	  BD[i].pos = (BD[i].m*BD[i].pos + BD[j].m*BD[j].pos)/(BD[i].m + BD[j].m);
	  BD[i].vel = (BD[i].m*BD[i].vel + BD[j].m*BD[j].vel)/(BD[i].m + BD[j].m);
	  BD[i].m   = (BD[i].m + BD[j].m);
	  BD[i].r   = pow( BD[i].m/rho, 1.0/3 );
	  
	  //pos[i] = ( m[i]*pos[i] + m[j]*pos[j] ) / ( m[i] + m[j] );
	  //vel[i] = ( m[i]*vel[i] + m[j]*vel[j] ) / ( m[i] + m[j] );
	  //m[i]   = m[i] + m[j];
	  //r[i]   = pow( m[i]/rho, 1.0/3 );

	  BD[j].pos = BD[N-1].pos;
	  BD[j].vel = BD[N-1].vel;
	  BD[j].m   = BD[N-1].m;
	  BD[j].r   = BD[N-1].r;
	  
	  // pos[j] = pos[N-1]; //as N is total planets, index is N-1
	  //vel[j] = vel[N-1];
	  //m[j]   = m[N-1];
	  //r[j]   = r[N-1];
	  N--;
	  a++;
      }
    }
  }
}
  
//=======================

double leap( planet BD[] ){//vector pos[], vector vel[], double m[], double r[]){

  vector a[N];
  int i;

  for( i=0; i<N; i++ ){//why do i have to do this?? cuz u suck
    a[i].x = 0;
    a[i].y = 0;
    a[i].z = 0;
  }

  planetRegion(BD);//finds the regions a planet is in
  
  regionCOM(BD);
  
  posHalfStep(BD);
  //posHalfStep(pos,vel);

  Force(BD,a);
  //Force(pos,F,a,m);
  
  velFullStep(a,BD);
  //velFullStep(a,vel);
  
  posHalfStep(BD);
  //posHalfStep(pos,vel);

  rElem(BD);
  //rElem(pos,vel,m,r);
}  


//===================

double energy( planet BD[] ){//vector pos[], vector vel[], double m[]){
  double KE = 0,PE = 0;
  int i,j;
  
  for( i=0; i<N; i++ ){
    KE += 0.5*BD[i].m*(BD[i].vel*BD[i].vel);
    // KE += 0.5*m[i]*(vel[i]*vel[i]);
  }

  for( i=0; i<N; i++){
    for( j=i+1; j<N; j++ ){
      PE += BD[i].m*BD[j].m / radius(BD[i].pos,BD[j].pos);
      //PE += m[i] * m[j] / radius(pos[i],pos[j]);
    }
  }
  return KE - G*PE;
}

//==================
double mass( planet BD[] ){//double m[]){
  
  int i;
  double mass=0.0;
  for( i=0; i<N; i++ ){
    mass += BD[i].m;
    //mass += m[i];
  }
  return mass;
}

//==================

double mom( planet BD[] ){//vector vel[], double m[]){
//more like urmom
  int i;
  vector p;
  p.x = 0;
  p.y = 0;
  p.z = 0;
  
  for( i=0; i<N; i++ ){
    p = p + BD[i].m*BD[i].vel;
    //p = p + m[i]*vel[i];
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

  //vector pos[N],vel[N],vcm;
  planet BD[N];
  vector vcm;
  int i;
  int frame = 0;

  // for( i=0; i<N; i++ ){//how to initialize it??
  //   BD[i].level[LVL];
  // }
  
  vcm.x = 0;
  vcm.y = 0;
  vcm.z = 0;
  
  //double m[N],r[N];
  double t,theta,rho=20;

  
  

  for( i=0; i<N; i++){

    BD[i].pos.x = L*drand48();//(drand48()-.5);not neg so fits planetRegion
    BD[i].pos.y = L*drand48();//(drand48()-.5);
    BD[i].pos.z = L*drand48();//(drand48()-.5);

    printf("x = %.2e, y = %.2e, z = %.2e\n",BD[i].pos.x, BD[i].pos.y, BD[i].pos.z);

    BD[i].vel.x = v*(drand48()-.5);
    BD[i].vel.y = v*(drand48()-.5);
    BD[i].vel.z = v*(drand48()-.5);

    BD[i].r   = 0.1;
    BD[i].m  = rho*BD[i].r*BD[i].r*BD[i].r;
    
    //pos[i].x = L*(drand48()-.5);
    //pos[i].y = L*(drand48()-.5);
    //pos[i].z = L*(drand48()-.5);
    
    //vel[i].x = v*(drand48()-.5);
    //vel[i].y = v*(drand48()-.5);
    //vel[i].z = v*(drand48()-.5);

    //r[i] = 0.1;
    //m[i] = rho*r[i]*r[i]*r[i];
  }

  for( i=0; i<N; i++ ){
    vcm = vcm + BD[i].m*BD[i].vel;
    //vcm = vcm + m[i]*vel[i];
  }

  vcm = vcm/mass( BD );
  //vcm = vcm/mass( m );
  
  for( i=0; i<N; i++ ){
    BD[i].vel = BD[i].vel - vcm;
    //vel[i] = vel[i] - vcm;
  }
      
  for( t=0; 1; t+=dt ){
    leap(BD);
    if( frame % fSkip == 0 ){
      printf("T -0.8 0.82\n E = %f : M = %f : P = %f : # Collisions = %d\n",energy(BD),mass(BD),mom(BD),a);
      //printf("T -0.8 0.82\n E = %f : M = %f : P = %f : # of Collisions = %d\n",energy(pos,vel,m),mass(m),mom(vel,m),a);
      
      for( i=0; i<N; i++){
	printf("c3 %e %e %e %e\n",BD[i].pos.x,BD[i].pos.y,BD[i].pos.z,BD[i].r);
	       //printf("c3 %e %e %e %e\n",pos[i].x,pos[i].y,pos[i].z,r[i]);
      }
      printf("F\n");
    }
    frame++;
  }
}
