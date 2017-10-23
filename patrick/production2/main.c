/*To compile, copy paste the following:

  g++ main.c forceError.c getArrayIndex.c intVec.c planetListFunctions.c populateRegions.c wellSepCoM.c vector.c parentChild.c omelyan.c resetregion.c -lm -fopenmp -o <output file name> -O4

 */
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "planet.h"
#include "pln.h"
#include "region.h"
#include "forceError.h"
#include "populateRegions.h"
#include "planetListFunctions.h"
//#include "LukeIAmYourParent.h"
#define THETA 1.35120719195966
#include "parentChild.h"
#include "omelyan.h"
//UNITS
//Mass   : Earth Mass
//Length : AU (astronomical unit)
//Time   : YEAR
//fix planet sizes, way too big

// link for Forest-Ruth Algorithm 4th order
//http://young.physics.ucsc.edu/115/leapfrog.pdf


//./FUNCTION 100 1e-3 100 2.5 20 15 1, this is the initial condition for the graphs
int LVL = 5; 
int numDel;
region *regions;
int totalIndeces;
double alpha2;
int collision = 0;
int nCollision = 0;
int parallel;


double fErr;//"Force error" - will store the total error incurred for a round of force calculations on all planets
double newtF;//"Newtonian Force"
double quickF;//"Quick Force" - our approximated force
//GLOBAL variable declaration
int N; //number of planets    
double dt; //time step
double G  = 4*M_PI*M_PI; //newtons constant in AU^3 yr^-2 M_sun^-1
double M_Sat = 1000000;
int a; //number of collisions
double L; //inital grid size
vector com; //center of mass to test my functions
//double THETA = 1/(2-pow(2,1./3));



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
/*void findNeighbors( int regionIndex, int neighbors[] ){

  int i,j,k;

//i is left right, so +-1
//j is up down, so +-n
//k is front back, so +-n^2

int n;//boxes per side on lowest level
int neighbor;

n = pow(2,LVL);

for(i=-1;i<2;i++){

for(j=-1;j<2;j++){
for(k=-1;k<2;k++){
if(i==0 && j==0 && k==0){
continue;
}

}*/

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

double rElem( planet BD[] ){//Check for collisions

	int i,j;
	double rad;
	double rho = 20;
	int other;
	pln *otherPln;
	int regionIndex;
	region regionOfInterest;
	int neighbors[26];//there are 3^3 regions in a 3x3 cube, -1 for the center
	collision = 0;
	nCollision = 0;
	for( i=0; i<N; i++ ){
		for( j=i+1; j<N; j++ ){
			rad = radius(BD[i].pos,BD[j].pos);

			if( rad < 0.8*(BD[i].r+BD[j].r) ){

				BD[i].pos = (BD[i].m*BD[i].pos + BD[j].m*BD[j].pos)/(BD[i].m + BD[j].m);//find new properties 
				BD[i].vel = (BD[i].m*BD[i].vel + BD[j].m*BD[j].vel)/(BD[i].m + BD[j].m);
				BD[i].m   = (BD[i].m + BD[j].m);
				//				BD[i].r   = BD[i].r*pow( BD[i].m/rho, 1.0/3 );

				BD[j].pos = BD[N-1].pos;//move information
				BD[j].vel = BD[N-1].vel;
				BD[j].m   = BD[N-1].m;
				BD[j].r   = BD[N-1].r;

				N--;
				a++;
				collision = 1;
				nCollision++;
			}
		}
	}
	return (nCollision);

}

double fLoop(planet *BD, int step)
{

	int i, p, thread, nThreads;
	double rad;
//	vector zero;
//	zero.x=0; zero.y=0; zero.z=0;

#pragma omp parallel for schedule(dynamic) if(parallel) num_threads(4)
	for( i=0; i<N; i++ )
	{
		/*
		   nThreads = omp_get_num_threads();
		   if(nThreads != 4)
		   printf("!nThreads = %d\n", nThreads);
		   thread = omp_get_thread_num();
		   printf("!.... Step %d .. BEGIN i=%d on thread %d\n", step, i, thread);
		 */
		//Particle-Particle
		forceMagic( regions[0], BD[i], BD, i);
		if (i==0)
                	printf("!Particle %d acc update on thread %d: acc is %e,%e,%e\n",i,omp_get_thread_num(),BD[i].acc.x, BD[i].acc.y, BD[i].acc.z);

		//		if(BD[i].acc.x != BD[i].acc.x || BD[i].acc.y != BD[i].acc.y || BD[i].acc.z != BD[i].acc.z)//If any component of .acc is nan
		//		{
		//			p = omp_in_parallel();
		//			printf("!.... ERROR: Step %d .. thread %d .. after FM .. p: %d .. BD[%d].acc = %f %f %f .. rad: %f\n", step, thread, p, i, BD[i].acc.x, BD[i].acc.y, BD[i].acc.z, rad);
		//		}


		//Saturn-Particle
		//		rad = radNew(BD[i].pos, zero);
		//		if(rad != rad)
		//			printf("!.... ERROR: Step %d .. thread %d .. rad = %f\n", step, thread, rad);
		//		BD[i].acc = BD[i].acc - G*M_Sat*BD[i].pos/(rad*rad*rad);	

		//		if(BD[i].acc.x != BD[i].acc.x || BD[i].acc.y != BD[i].acc.y || BD[i].acc.z != BD[i].acc.z)//If any component of .acc is nan
		//		{
		//			p = omp_in_parallel();
		//			printf("!.... ERROR Step %d .. thread %d .. after Sat .. p: %d .. BD[%d].acc = %f %f %f .. rad: %f\n", step, thread, p, i, BD[i].acc.x, BD[i].acc.y, BD[i].acc.z, rad);
		//		}
		//		printf("!.... Step %d .. END   i=%d on thread %d\n", step, i, thread);
	}	

}

//=======================

double omelyan( planet BD[] , double t){

	int i, p;//index variables
	double rad;
	vector zero;

	zero.x=0; zero.y=0; zero.z=0;

	r_step1(BD);

	reset_acc(BD);

	//STEP 1
	fLoop(BD, 1);
	v_step1(BD);
	r_step2(BD);
	reset_acc(BD);
	//	printf("\n");

	//STEP 2
	fLoop(BD, 2);
	v_step2(BD);
	r_step3(BD);
	reset_acc(BD);
	//	printf("\n");

	//STEP 3
	fLoop(BD, 3);
	v_step3(BD);
	r_step4(BD);
	reset_acc(BD);
	//	printf("\n");

	//STEP 4
	fLoop(BD, 4);
	v_stepf(BD);
	r_stepf(BD);
	//	printf("\n");

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
		forceMagic( regions[0], BD[i], BD, i);
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
/*
   double initStandard(planet *BD, int N)
   {
//FIXME How to pass pointer to BD, edit BD directly here?

int i;

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
return com;
}
 */

//==================

void initRing(int N, planet *BD, int rand)
{
	int i;
	double r, theta, vel;

	//Saturnogenesis
	BD[0].pos.x = 0;
	BD[0].pos.y = 0;
	BD[0].pos.z = 0;
	BD[0].vel.x = 0;
	BD[0].vel.y = 0;
	BD[0].vel.z = 0;
	BD[0].acc.x = 0;
	BD[0].acc.y = 0;
	BD[0].acc.z = 0;
	BD[0].r = 50;
	BD[0].m = M_Sat;
	BD[0].num = 0;


	if(rand==1)
	{
		for(i=1; i<N; i++)
		{
			r = 100+100*drand48();
			theta = 2*M_PI*drand48();

			BD[i].pos.x = r*cos(theta);
			BD[i].pos.y = r*sin(theta);
			BD[i].pos.z = drand48();

			vel = sqrt(G*M_Sat/r);

			BD[i].vel.x = -vel*sin(theta);
			BD[i].vel.y = vel*cos(theta);
			BD[i].vel.z = 0;
			//		printf("!i=%d  velx=%f  vely=%f\n", i, BD[i].vel.x, BD[i].vel.y);

			BD[i].acc.x = 0;
			BD[i].acc.y = 0;
			BD[i].acc.z = 0;

			BD[i].r = 1;
			BD[i].m = 1000*(.01+.001*drand48());
			//		printf("!BD[%d].m = %f\n", i, BD[i].m);
			BD[i].num = i;
		}
	}
	else
	{	
		for(i=1; i<N; i++)
		{
			r = 100+100*i/N;
			theta = 2*M_PI*i*20/N;

			BD[i].pos.x = r*cos(theta);
			BD[i].pos.y = r*sin(theta);
			BD[i].pos.z = drand48();

			vel = sqrt(G*M_Sat/r);

			BD[i].vel.x = -vel*sin(theta);
			BD[i].vel.y = vel*cos(theta);
			BD[i].vel.z = 0;
			//		printf("!i=%d  velx=%f  vely=%f\n", i, BD[i].vel.x, BD[i].vel.y);

			BD[i].acc.x = 0;
			BD[i].acc.y = 0;
			BD[i].acc.z = 0;

			BD[i].r = 1;
			BD[i].m = .001;
			//		printf("!BD[%d].m = %f\n", i, BD[i].m);
			BD[i].num = i;
		}
	}

}

//==================

int main(int nParam, char **paramList){

	omp_set_nested(0);

	//FROM CONFIG FILE: v_initial, box_width, initCond_width, alpha2
	char var[10], val[10];
	FILE *config=fopen("config.txt", "r");

	//Give initial values incase not specified in config file
	double v=1, scatter=1;
	L=1; alpha2=1;

	while( fscanf(config, "%s %s", var, val) != EOF)
	{
		if( strcmp(var, "v_initial") == 0) v=atof(val);
		if( strcmp(var, "box_width") == 0) L=atof(val);
		if( strcmp(var, "initCond_width") == 0) scatter=atof(val);
		if( strcmp(var, "alpha2") == 0) alpha2=atof(val);
	}


	//FROM CONSOLE: N, dt, tmax, frameskip, init_type, parallel, print
	if( nParam < 8 ) {printf("\n!Usage: ./<this> <N> <dt> <tmax> <frameskip> <init_type> <parallel> <print> < config.txt\n\n"); exit(0);}

	N             = atoi(paramList[1]);
	dt            = atof(paramList[2]);
	double tmax   = atof(paramList[3]);
	int fSkip     = atoi(paramList[4]);
	int init_type = atoi(paramList[5]);
	parallel      = atoi(paramList[6]);
	int print     = atoi(paramList[7]);

	printf("!IC: N=%d  dt=%f  tmax=%f  fSkip=%d  init_type=%d  parallel=%d  print=%d  v_initial=%f  scatter=%f  L=%f  alpha2=%f\n", N, dt, tmax, fSkip, init_type, parallel, print, v, scatter, L, alpha2);

	//MORE DECLARATIONS
	vector vcm;
	int i;
	totalIndeces = 0;
	int frame = 0;
	double MASS = 0;
	double Ei,Pi;
	planet BD[N];
	int nCol;
	struct timespec start, end;
	uint64_t tTot;
	double tDec;


	//INITIALIZING
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

	//Initializer used to be here. Moved to separate function.	
	initRing(N, BD, init_type);

	//	forceMagic(regions[0], BD[0], BD);

	printf("!Finished initRing\n");

	for( i=0; i<N; i++ ){
		//		printf("BD[%d].m = %f\n", i, BD[i].m);
		MASS += BD[i].m;
		//		printf("MASS = %f\n", MASS);
	}

	printf("!Finished initializing planet locations\n");

	for( i=0; i<N; i++ ){
		vcm = vcm + BD[i].m*BD[i].vel;
	}

	com = com/MASS;
	printf("!The center of mass is (%f,%f,%f)\n",com.x,com.y,com.z);

	vcm = vcm/MASS;

	//for( i=0; i<N; i++ ){
	//   BD[i].vel = BD[i].vel - vcm;
	// }

	//PREP FILES FOR WRITING DATA
	alpha2 = 900;
	Ei = energy(BD);
	Pi = mom(BD);
	printf("Pi = %f, Ei = %f\n",Pi,Ei);
	//where the magic happens
	FILE * E_vs_t = fopen("E_vs_t_alpha2_900.txt","w");
	FILE * E_vs_t_c = fopen("E_vs_t_alpha2_900_c.txt","w");//if there is a collision, print to another file
	FILE * P_vs_t = fopen("P_vs_t_alpha2_900.txt","w");
	FILE * P_vs_t_c = fopen("P_vs_t_alpha2_900_c.txt","w");
	FILE * mom_vs_t = fopen("mom_vs_t_500_1e4.txt","w");
	FILE * nColl = fopen("nCol_vs_t.txt","w");

	//NOW RUN
	clock_gettime(CLOCK_MONOTONIC_RAW, &start);//Get starting time
	for( t=0; t<tmax; t+=dt )
	{
		printf("!----------TIME: %f----------\n", t);
		populateRegions(BD, totalIndeces);//find region properties
		numDel=0;
#pragma omp parallel for schedule(dynamic)
		for( i=0; i<totalIndeces; i++ )//clear lists
		{
			//	printf("!Clearing list of region %d\n", i);
			clearList(&(regions[i].planets));//if we call the pointer .list
		}

		popLists(BD);//populate lists
		omelyan(BD, t);
		//		nCol = rElem(BD);
		printf("!Saturn.pos = %f  %f  %f\n",BD[0].pos.x, BD[0].pos.y, BD[0].pos.z);
		printf("!Saturn.vel = %f  %f  %f\n",BD[0].vel.x, BD[0].vel.y, BD[0].vel.z);
		printf("!Saturn.acc = %f  %f  %f\n",BD[0].acc.x, BD[0].acc.y, BD[0].acc.z);
		if(print==1)
		{
			if( frame % fSkip == 0 || collision == 1)
			{
				// printf("t is %f\n",t);
				//printf("!%d %f\n",collision, t);
				if(collision == 0 )
				{
					//printf("!time = %f\n",t);
					fprintf(E_vs_t,"%e %e\n",t,(Ei-energy(BD))/Ei);
					fprintf(P_vs_t,"%e %e\n",t,(Pi-mom(BD))/Pi);
				}
				else
				{
					fprintf(E_vs_t_c,"%e %e\n",t,(Ei-energy(BD))/Ei);
					fprintf(P_vs_t_c,"%e %e\n",t,(Pi-mom(BD))/Pi);
				}
				if( t>0 )
				{
					printf("T -0.8 0.82\nE = %e : P = %e : # Collisions = %d : N = %d : t = %f\n",energy(BD),mom(BD),a,N,t);
					fprintf(mom_vs_t,"%e %e\n", t, mom(BD));
					if(nCol > 0)
					{
						fprintf(nColl, "%e %d\n", t, nCol);
					}
					for( i=0; i<N; i++)
					{
						printf("c3 %e %e %e %e\n",BD[i].pos.x,BD[i].pos.y,BD[i].pos.z,BD[i].r);
					}
					//					printf("c3 0.0 0.0 0.0 50\n");
					printf("F\n\n");
				}
			}
			frame++;
		}

	}
	clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	tTot = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
	tDec = tTot;
	tDec /= 1000000;
	printf("!Total time: %fsec\n", tDec);
	fclose(E_vs_t);
	fclose(P_vs_t);
	fclose(E_vs_t_c);
	fclose(P_vs_t_c);
}
