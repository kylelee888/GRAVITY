//On level 7, it seems like force is no longer taking the longest. Check to make
//sure that everything is recursing to eliminate looping over all regions
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
//#include "LukeIAmYourParent.h"

#include "parentChild.h"
#include "omelyan.h"
#include "timing.h"

//UNITS
//Mass   : Earth Mass
//Length : AU (astronomical unit)
//Time   : YEAR
//fix planet sizes, way too big


//./FUNCTION 100 1e-3 100 2.5 20 15 1, this is the initial condition for the graphs
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
  
  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++){//the i+1 is what makes it faster instead of j=0 and going through to N

      rad = radius(BD[i].pos,BD[j].pos);
      F   = G * BD[j].m * BD[i].m * (BD[j].pos-BD[i].pos) / (rad*rad*rad);
      BD[i].acc = BD[i].acc + F/BD[i].m;
      BD[j].acc = BD[j].acc - F/BD[j].m;
    }
  }
}
//=================================

void compactify(planet BD[]){
  
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

int rElemNew( planet BD[] ){

  int i,j;
  int ndel=0;
  double rho = 20;
  double rad;

  for(int k=0; k<collision_number; k++){
    i = collision_pair_1[k];
    j = collision_pair_2[k];
    if(i != j){
      rad = radius(BD[i].pos, BD[j].pos);
      if(rad < (BD[i].r + BD[j].r)){
	//new properties of collided planets
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
    }
  }
  if(ndel) compactify(BD);
  return ndel;
}

double rElem( planet BD[] ){

  int i,j;
  double rad;
  double rho = 20;
  int other;
 
  int regionIndex;
 
 
 
  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++ ){
      rad = radius(BD[i].pos,BD[j].pos);
      
      if( rad < 0.8*(BD[i].r+BD[j].r) ){
	
	BD[i].pos = (BD[i].m*BD[i].pos + BD[j].m*BD[j].pos)/(BD[i].m + BD[j].m);//find new properties 
	BD[i].vel = (BD[i].m*BD[i].vel + BD[j].m*BD[j].vel)/(BD[i].m + BD[j].m);
	BD[i].m   = (BD[i].m + BD[j].m);
	BD[i].r   = pow( BD[i].m/rho, 1.0/3 );

	BD[j].num = BD[N-1].num;
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

double omelyan( planet BD[], int method ){

  int i;
  long int force = 0;
  //need to calc f_apprx and f_exact for all planets
  
  
  r_step1(BD);

  reset_acc(BD);
  starttimer(2);
  if( method == 1 ){
    for( i=0; i<N; i++ ){
      forceMagic( regions[0], BD[i], BD, alpha );
    }
  }
  else{
    Force(BD);
  }
  force+=stoptimer(2);
  v_step1(BD);

  r_step2(BD);

  reset_acc(BD);
  starttimer(2);
  if( method == 1 ){
    for( i=0; i<N; i++ ){
      forceMagic( regions[0], BD[i], BD, alpha );
    }
  }
  else{
    Force(BD);
  }
  force+=stoptimer(2);
  v_step2(BD);

  r_step3(BD);

  reset_acc(BD);
  starttimer(2);
  if( method == 1 ){
    for( i=0; i<N; i++ ){
      forceMagic( regions[0], BD[i], BD, alpha );
    }
  }
  else{
    Force(BD);
  }
  force+=stoptimer(2);
  v_step3(BD);

  r_step4(BD);

  reset_acc(BD);
  collision_check  = 1;
  collision_number = 0;
  starttimer(2);
  if( method == 1 ){
    for( i=0; i<N; i++ ){
      forceMagic( regions[0], BD[i], BD, alpha );
    }
  }
  else{
    Force(BD);
  }
  force+=stoptimer(2);
  collision_check = 0;
  v_stepf(BD);

  r_stepf(BD);
  return (double) force;
}


double leap( planet BD[] ){
//  printf("!LEAP: Start\n");  
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
    //printf("!  Computing force on planet %d\n",i);
    forceMagic( regions[0], BD[i], BD, alpha );
  }
  
  velFullStep(BD);
    
  posHalfStep(BD);
}  

double calc_dF( planet BD[], vector F_exact[], double alpha ){
  
  int i;
  double dF = 0;
  double Ftot = 0;
  //vector F_apprx[N];
  //vector F_exact[N];
  // vector F_exact[N];
  vector F_apprx[N];

  reset_acc(BD);
  
  for( i = 0; i<N; i++ ){
    forceMagic(regions[0], BD[i], BD, alpha);
  }
  for(i=0;i<N;i++){
    //F_apprx += (BD[i].acc*BD[i].acc)*(BD[i].m*BD[i].m);
    F_apprx[i] = (BD[i].acc*BD[i].m);
  }
  
  for(i=0;i<N;i++){
    dF += sqrt((F_exact[i]-F_apprx[i])*(F_exact[i]-F_apprx[i]));//this is a dot product
    Ftot += sqrt(F_exact[i]*F_exact[i]);
  }
  return dF/Ftot;///Ftot; 
}

double calc_dt( planet BD[], double alpha ){
  
  int i;
  double ratio;
  long int exact_t;
  long int apprx_t;
 
  
  starttimer(0);
  Force(BD);
  exact_t = stoptimer(0);
  
  reset_acc(BD);

  starttimer(0);
  for( i = 0; i<N; i++ ){
    forceMagic(regions[0], BD[i], BD, alpha);
  }
  apprx_t = stoptimer(0);

  reset_acc(BD);
  ratio = apprx_t/exact_t;
  return ratio;
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

  starttimer(1);
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
  double delta_F;//stores total difference in force^2
  double ratio_times;

  double force_t;
  long int apprx_time;
  long int exact_time;
  long int pop_regions_t;
  long int resetting_t;
  long int pop_list_t;
  long int collision_t;
  long int omelyan_t;
  long int anim_t;
  long int setup;
  vector F_exact[N];

  for(i=0; i<LVL; i++)
  {
   totalIndeces += pow(8,i);
  }
  regions=(region*)malloc(sizeof(region)*totalIndeces);//FIXME Is this correct?
  //printf("!Number of planets: %d\n", N);
  // printf("!Frameskip is %d\n",fSkip);
  // printf("!Total # of indeces in region array: %d\n", totalIndeces);
    
  findKids(LVL);

  vcm.x = 0;
  vcm.y = 0;
  vcm.z = 0;
  
  double t,theta,rho=20;
  com.x = 0;//to check the COM finding of my function
  com.y = 0;
  com.z = 0;
  int q;
  int count = 0;
      
  for( i=0; i<N; i++){
    
    BD[i].num = i;
    
    BD[i].pos.x = L * 0.5 + (drand48()-0.5)*scatter; //not neg so fits planetRegion
    BD[i].pos.y = L * 0.5 + (drand48()-0.5)*scatter; //not neg so fits planetRegion
    BD[i].pos.z = L * 0.5 + (drand48()-0.5)*scatter; //not neg so fits planetRegion
    
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
    
  //printf("!Finished initializing planet locations\n");
  
  for( i=0; i<N; i++ ){
    vcm = vcm + BD[i].m*BD[i].vel;
  }
  
  com = com/MASS;
  //printf("!The center of mass is (%f,%f,%f)\n",com.x,com.y,com.z);
  
  vcm = vcm/MASS;
  
  for( i=0; i<N; i++ ){
    BD[i].vel = BD[i].vel - vcm;
  }
  
  Ei = energy(BD);
  Pi = mom(BD);
  //printf("!Pi = %f, Ei = %f\n",Pi,Ei);
  
  
  //printf("It took %f seconds to setup\n", (double) stoptimer(1)/1000000);
  //FILE * t_vs_alpha = fopen("t_vs_alpha_N5000_LVL5.txt","w");
  for( t=0; 1; t+=dt ){
    
    //starttimer(1);
    reset_region(regions[0]);
    //resetting_t = stoptimer(1);
    
    //starttimer(1);
    populateRegions(BD, totalIndeces);//find region properties
    //pop_regions_t = stoptimer(1);
    //printf("!populateRegions Finished\n");
    
    //printf("!-----Total planets cleared: %d\n", numDel);
    //starttimer(1);
    popLists(BD);//populate lists
    //pop_list_t = stoptimer(1);
    //printf("!popLists Finished\n");
    //leap(BD);
        
    /*
    if(t == 0){//only do it once per simulation
      reset_acc(BD);
      Force(BD);
      for(i=0;i<N;i++){
	F_exact[i] = BD[i].m*BD[i].acc;
      }
      for( alpha = 1; alpha<10; alpha += 0.25 ){
	//printf("%f\n",alpha);
	delta_F = calc_dF(BD, F_exact, alpha);//fuck it, make it so you pass alpha, wellSepCOM does the multiplication
	printf("%e %e\n", alpha, delta_F);
      }
      break;
    */
      /*
	for( alpha = 2; alpha<100; alpha +=2 ){
	ratio_times = calc_dt(BD, alpha);
	}
	break;
      */ 
    //}
    
    // starttimer(1);
    // Force(BD);
    //exact_time = stoptimer(1);

    /*
    for(alpha = 1; alpha < 10; alpha += 0.25){
      starttimer(1);
      for(i=0;i<N;i++){
	forceMagic(regions[0], BD[i], BD, alpha);
      }
      apprx_time = stoptimer(1);
      printf("%f %f\n",alpha, (double) apprx_time/exact_time);
      //printf("approx / exact = %ld / %ld\n",apprx_time, exact_time);
    }
    break;
    */
    //force_t = omelyan(BD, 1);
    
    //omelyan_t = stoptimer(1);
    
    //starttimer(1);
    a+=rElemNew(BD);
    //collision_t = stoptimer(1);
    
    //printf("apprx time = %ld\n",apprx_time[q]);
    //printf("tot = %ld, omel = %ld, force = %f\n",apprx_time[q], omelyan_t, force_t);
    //printf("!%ld %ld\n", apprx_time[q], pop_regions_t + resetting_t + pop_list_t  + omelyan_t + collision_t + anim_t );
    //printf("!Pop Reg = %f, Clr List = %f, Pop List = %f, Int = %f, Coll Check = %f Anim = %f\n",
    //  (double) pop_regions_t, (double) resetting_t, (double) pop_list_t,
    // (double) omelyan_t, (double) collision_t, (double) anim_t);
  

  
  
  
  
  
  
    starttimer(1);
    if( frame % fSkip == 0 ){//|| collision == 1){
    
    printf("T -0.8 0.82\nE = %f : P = %f : # Collisions = %d : t = %f\n",energy(BD),mom(BD),a,t);
    for( i=0; i<N; i++){
    printf("c3 %e %e %e %e\n",BD[i].pos.x,BD[i].pos.y,BD[i].pos.z,BD[i].r);
    }
    printf("F\n"); 
    }
    anim_t = stoptimer(1);
  
    frame++;
  
    
    
    
    //printf("!Pop Reg = %f, Clr List = %f, Pop List = %f, Int = %f, Coll Check = %f Anim = %f\n",
    //	    (double) pop_regions_t/timestep, (double) resetting_t/timestep, (double) pop_list_t/timestep,
    //	    (double) omelyan_t/timestep, (double) collision_t/timestep, (double) anim_t/timestep);
  
  
  
  
  
  
    
    //check that I am getting all the times
    //fprintf(t_vs_alpha,"%e %e\n", alpha, (double) apprx_time/exact_time);
    //printf("approx : N = %d, alpha = %f, time = %f\n",N, alpha, (double) apprx_time/1000);
    //count++;
    //printf("FORCE : Exact Time is %ld microseconds, Apprx Time is %ld microseconds\n", exact_time, apprx_time);
    //printf("t_exact/t_apprx = %f\n",exact_time/apprx_time);
    
    //break;
    
    //break;
    
  }
}
// count++;
 
  //printf("done 5000\n");
  //fprintf(dF_vs_alpha,"%e %e\n", alpha, delta_F);
  //}
  //fclose(t_vs_alpha);
  //printf("Im done\n");



