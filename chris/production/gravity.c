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
#include "omp.h"
//#include "LukeIAmYourParent.h"
#define THETA 1.35120719195966
#include "parentChild.h"
#include "timing.h"
#include "resetregion.h"

//UNITS
//Mass   : Earth Mass
//Length : AU (astronomical unit)
//Time   : YEAR
//fix planet sizes, way too big

// link for Forest-Ruth Algorithm 4th order
//http://young.physics.ucsc.edu/115/leapfrog.pdf

int LVL = 7; 
int numDel;
region *regions;
int totalIndeces;
long int exact_time, approx_time;

//GLOBAL variable declaration
int N; //number of planets    
double G  = 4*M_PI*M_PI; //newtons constant 
int a; //number of collisions
double L; //inital grid size
double dt;
vector com; //center of mass to test my functions
//double THETA = 1/(2-pow(2,1./3));
int exact,approx;
bool collide_check;
int coll1[10000],coll2[10000];
int collnum=0;



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
  int other;
  pln *otherPln;
  int regionIndex;
  region regionOfInterest;
  int neighbors[26];//there are 3^3 regions in a 3x3 cube, -1 for the center

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
  
void compactify(planet BD[])
{
  for (int i=0;i<N;i++)
  {
    if (BD[i].m==0)
    {
      N--;
      BD[i].pos=BD[N].pos;
      BD[i].vel=BD[N].vel;
      BD[i].m=BD[N].m;
      BD[i].r=BD[N].r;
    }
  }
}

int rElemNewNew( planet BD[] )
{
  int i,j;
  int ndel=0;
  double rad;
  double rho = 20;
  int other;
  pln *otherPln;
  region regionOfInterest;

  int count;

  // someone else made us a list of planets that might collide. we just need to go through it
  for (int k=0; k<collnum; k++)
  {
    i=coll1[k]; j=coll2[k];
    if(i != j){//j is not BD[i]
      rad = radius(BD[i].pos,BD[j].pos);
      if(rad < 1*(BD[i].r + BD[j].r)){
	//find new properties of collided planets
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
  if (ndel) compactify(BD); // if anyone was marked to be deleted, compactify the list
  return ndel;
}


double rElemNew( planet BD[] )
{
  int i,j;
  double rad;
  double rho = 20;
  int other;
  pln *otherPln;
  region regionOfInterest;
  double maxrad=0;

  /*for( i=0; i<N; i++ ){
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
    }*/

  int ntotest=0;
  int listtotest[27];
  
  //looping through all regions is ridiculous
  for(i=0;i<N;i++){
    if (BD[i].r > maxrad) maxrad=BD[i].r;
  }


  int count;
  for(i=0;i<N;i++){

    // Let's do a trick here: compute the region of ourself, shifted by the maximum radius of any planet plus our own radius, in each of 
    // the 27 cardinal and semicardinal directions. 
    // We then compute which lowest-level region that lies in; if it is different from our own, we worry about colliding against 
    // planets in it. Otherwise, we don't.
    count = 0;
    vector vtest; 
    int ourregion = BD[i].level[LVL-1]; // we know this one
    double rstep=BD[i].r+maxrad;

    int thisregion = getregion(BD[i].pos, LVL-1);
    ntotest=0;
    int alreadythere;


    for (vtest.x=BD[i].pos.x-rstep; vtest.x < BD[i].pos.x+rstep+1e-5; vtest.x += rstep)
      for (vtest.y=BD[i].pos.y-rstep; vtest.y < BD[i].pos.y+rstep+1e-5; vtest.y += rstep)
        for (vtest.z=BD[i].pos.z-rstep; vtest.z < BD[i].pos.z+rstep+1e-5; vtest.z += rstep)
	{
	  int thisregion = getregion(vtest, LVL-1);
	  if ((thisregion != ourregion))
	  {
            // see if it's already marked as something to test against
            alreadythere=0;
            for (int j=0; j<ntotest; j++) 
              if (listtotest[j]==thisregion) alreadythere=1;
            
            if (alreadythere == 0)
            {
              listtotest[ntotest]=thisregion;
              ntotest++;
            }
          }
        }
   if (regions[ourregion].numPln > 1)
   {
     listtotest[ntotest] = ourregion; // always have to test against our region
     ntotest++;
   }
   for (int j=0; j<ntotest; j++)
   {
   //         printf("!Testing: planet %d vs. region %d out of %d (home region %d)\n",i,listtotest[j],ntotest,ourregion);
	    regionOfInterest=regions[listtotest[j]];; 


//	        regionOfInterest = regions[BD[i].level[LVL-1]];//region our planet is in
	    if(regionOfInterest.numPln > 0){//if there is another planet to collide with
	      otherPln = regionOfInterest.planets;//point to the list of planets
	      while(otherPln != NULL){//loop through all planets in that regions list
		other = otherPln->plnNum;
		if(i != other){//other is not BD[i]
		  rad = radius(BD[i].pos,BD[other].pos);
		  if(rad < 1*(BD[i].r + BD[other].r)){
		    //find new properties of collided planets
		    BD[i].pos = (BD[i].m*BD[i].pos + BD[other].m*BD[other].pos)/(BD[i].m + BD[other].m); 
		    BD[i].vel = (BD[i].m*BD[i].vel + BD[other].m*BD[other].vel)/(BD[i].m + BD[other].m);
		    BD[i].m   = (BD[i].m + BD[other].m);
		    BD[i].r   = pow( BD[i].m/rho, 1.0/3 );

		    //shift planets
		    BD[other].pos = BD[N-1].pos;
		    BD[other].vel = BD[N-1].vel;
		    BD[other].m   = BD[N-1].m;
		    BD[other].r   = BD[N-1].r;

		    N--;
		    a++;
		  }
		}
		otherPln = otherPln -> nextPln;
	      }
       }
    }
  }
}

void get_accels(planet BD[], double a2)
{
  long int elapsed;
  static double G=4*M_PI*M_PI;
  starttimer(0);
  if (a2 == 0)
  {
    vector F;
    // compute accelerations using the naive algorithm, for performance comparison
    for(int i=0; i<N; i++ ) BD[i].acc.x=BD[i].acc.y=BD[i].acc.z=0;

    for(int i=0; i<N; i++ )
    {
      for(int j=i+1; j<N; j++ )
      {
	double rad=radius(BD[i].pos, BD[j].pos);
	F = (BD[i].pos-BD[j].pos) * G * BD[i].m * BD[j].m / (rad*rad*rad);
	BD[i].acc = BD[i].acc - F/BD[i].m;
	BD[j].acc = BD[j].acc + F/BD[j].m;
      }
    }
    elapsed=stoptimer(0);
    printf("!Accelerations computed in %ld μs (%.2f/particle) using naive method\n",elapsed,(float)elapsed/N),a2;
  }
  else
  {
//#pragma omp parallel for
  for(int i=0; i<N; i++ )
    {
      forceMagic( regions[0], BD[i], BD, a2 );
      BD[i].acc = BD[i].acc * G;
    }
    elapsed=stoptimer(0);
    printf("!Accelerations computed in %ld μs (%.2f/particle) using a2=%f\n",elapsed,(float)elapsed/N,a2);
  }
}
 
//=======================

double omelyan( planet BD[], double alpha2 ){

  int i;

  //reset our accelerations
  //foward step in position
  //force calculation
  //forward step in velocity
  //reset accelerations
  //backward step in position
  //force calculation
  //backward step in velocity
  //reset accelerations
  //backward step in position
  //force calculation
  //forward step in velocity
  //forward step in position
  reset_acc(BD);

  pos4_FWD(BD);
  get_accels(BD, alpha2);
  
  vel4_FWD(BD);

  reset_acc(BD);

  pos4_BCK(BD);
  get_accels(BD, alpha2);

  vel4_BCK(BD);

  reset_acc(BD);

  pos4_BCK(BD);

  // we are also going to check for collisions here:
  collide_check=1;
  collnum=0;
  get_accels(BD, alpha2);
  collide_check=0;

  vel4_FWD(BD);

  pos4_FWD(BD);
}


double leap( planet BD[], double alpha2){
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
    forceMagic( regions[0], BD[i],BD,alpha2 );
  }
  //FIXME: change planet struct to have .acc
  velFullStep(BD);
    
  posHalfStep(BD);
}  


//===================

double energy( planet BD[] ){
  return 0;
  starttimer(2);
  double KE = 0,PE = 0;
  static double radmin=1e10,r;
  int i,j;
  
  for( i=0; i<N; i++ ){
    KE += 0.5*BD[i].m*(BD[i].vel*BD[i].vel);
  }

  for( i=0; i<N; i++){
    for( j=i+1; j<N; j++ ){
      r=radius(BD[i].pos,BD[j].pos);
      PE += BD[i].m*BD[j].m / r;
      if (r/(BD[i].r + BD[j].r)<radmin) radmin=r/(BD[i].r + BD[j].r);
    }
  }
  printf("!Energy computed in %ld μs; minimum radius is %.2e\n",stoptimer(2),radmin);
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

  if( nParam < 9 ) {printf("\n!Usage: ./<this> <N> <dt> <frameskip> <initial velocity> <box width> <initial conditions width> <alpha2> <maxlevel> \n\n"); exit(0);}

  N          = atoi(paramList[1]);
  dt  = atof(paramList[2]);
  int fSkip  = atoi(paramList[3]);
  double v   = atof(paramList[4]);
  L          = atof(paramList[5]);
  double scatter  = atof(paramList[6]);
  double alpha2     = atof(paramList[7]);  
  LVL     = atoi(paramList[8]);  

  planet BD[N];
  vector vcm;
  int i;
  totalIndeces = 0;
  int frame = 0;
  double MASS = 0;

  for(i=0; i<LVL; i++)
  {
   totalIndeces += pow(8,i);
  }
  regions=(region*)malloc(sizeof(region)*totalIndeces);//FIXME Is this correct?
  printf("!Number of planets: %d\n", N);
  printf("!Total # of indeces in region array: %d\n", totalIndeces);
    
  findKids(LVL);

  vcm.x = 0;
  vcm.y = 0;
  vcm.z = 0;
  
  double t,theta,rho=20;
  com.x = 0;//to check the COM finding of my function
  com.y = 0;
  com.z = 0;
  
  for (i=0; i<N; i++) BD[i].num=i;
  
  for( i=0; i<N; i++){
    
    BD[i].pos.x = L * 0.5 + (drand48()-0.5)*scatter; //not neg so fits planetRegion
    BD[i].pos.y = L * 0.5 + (drand48()-0.5)*scatter; //not neg so fits planetRegion
    BD[i].pos.z = L * 0.5 + (drand48()-0.5)*scatter; //not neg so fits planetRegion

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
    starttimer(1);
    numDel=0;
    reset_region(regions[0]);
    printf("! ----- resetregions cleared %d planets\n",numDel);
    printf("! --- resetregions completed in %ld μs\n",stoptimer(1));
    starttimer(1);
    populateRegions(BD, totalIndeces);//find region properties
    printf("! --- popregions completed in %ld μs\n",stoptimer(1));
    //printf("!populateRegions Finished\n");
    numDel=0;


    for( i=0; i<totalIndeces; i++ ){//clear lists
//	printf("!Clearing list of region %d\n", i);
//      clearList(&(regions[i].planets));//if we call the pointer .list
    }    
    //printf("!-----Total planets cleared: %d\n", numDel);
    starttimer(1);
    popLists(BD);//populate lists
    printf("! --- poplists completed in %ld μs\n",stoptimer(1));
    
    //printf("!popLists Finished\n");
    //leap(BD);
    exact_time=0;
    starttimer(1);
    printf("! ------ starting integrator\n");
    omelyan(BD, alpha2);
    printf("! --- Timestep completed in %ld μs (%ld on exact): exact %d approx %d\n",stoptimer(1),exact_time,exact,approx);
    exact=approx=0;

    a+=rElemNewNew(BD);
    exit(0);
   
    if( frame % fSkip == 0 ){
    starttimer(1);
      printf("T -0.8 0.82\nE = %f : P = %f : # Collisions = %d\n",energy(BD),mom(BD),a);
      for( i=0; i<N; i++){
	printf("c3 %e %e %e %e\n",BD[i].pos.x,BD[i].pos.y,BD[i].pos.z,BD[i].r);
      }
      printf("F\n");
    printf("! --- anim output completed in %ld μs\n",stoptimer(1));
    }
    frame++;
  }
}

