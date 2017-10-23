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
int LVL = 5; 
int numDel;
region *regions;
int totalIndeces;
double alpha;
int collision = 0;
int collision_check;
int collision_number;
int collision_pair_1[10000];
int collision_pair_2[10000];

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

double omelyan( planet BD[] ){

	int i;
	//need to calc f_apprx and f_exact for all planets


	r_step1(BD);

	reset_acc(BD);
	for( i=0; i<N; i++ ){
		forceMagic( regions[0], BD[i], BD, alpha );
		if(BD[i].pos.x != BD[i].pos.x || BD[i].pos.y != BD[i].pos.y || BD[i].pos.z != BD[i].pos.z)
		{
			printf("--1--BD[%d]: %e %e %e\n", i, BD[i].pos.x, BD[i].pos.y, BD[i].pos.z);
		}
	}

	v_step1(BD);

	r_step2(BD);

	reset_acc(BD);
	for( i=0; i<N; i++ ){
		forceMagic( regions[0], BD[i], BD, alpha );
		if(BD[i].pos.x != BD[i].pos.x || BD[i].pos.y != BD[i].pos.y || BD[i].pos.z != BD[i].pos.z)
		{
			printf("--2--BD[%d]: %e %e %e\n", i, BD[i].pos.x, BD[i].pos.y, BD[i].pos.z);
		}
	}

	v_step2(BD);

	r_step3(BD);

	reset_acc(BD);
	for( i=0; i<N; i++ ){
		forceMagic( regions[0], BD[i], BD, alpha );
		if(BD[i].pos.x != BD[i].pos.x || BD[i].pos.y != BD[i].pos.y || BD[i].pos.z != BD[i].pos.z)
		{
			printf("--3--BD[%d]: %e %e %e\n", i, BD[i].pos.x, BD[i].pos.y, BD[i].pos.z);
		}
	}

	v_step3(BD);

	r_step4(BD);

	reset_acc(BD);
	collision_check  = 1;
	collision_number = 0;
	for( i=0; i<N; i++ ){//check for collisions during the last force step to save time
		forceMagic( regions[0], BD[i], BD, alpha );	
		if(BD[i].pos.x != BD[i].pos.x || BD[i].pos.y != BD[i].pos.y || BD[i].pos.z != BD[i].pos.z)
		{
			printf("--4--BD[%d]: %e %e %e\n", i, BD[i].pos.x, BD[i].pos.y, BD[i].pos.z);
		}
	}
	collision_check = 0;
	v_stepf(BD);

	r_stepf(BD);
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

double calc_dF( planet BD[], double alpha ){

	int i;
	double dF = 0;
	vector F_apprx[N];
	vector F_exact[N];

	reset_acc(BD);

	Force(BD);
	for(i=0; i<N; i++){
		F_exact[i] = BD[i].acc*BD[i].m;
	}

	reset_acc(BD);

	for( i = 0; i<N; i++ ){
		forceMagic(regions[0], BD[i], BD, alpha);
	}
	for(i=0;i<N;i++){
		F_apprx[i] = BD[i].acc*BD[i].m;
	}

	for(i=0;i<N;i++){
		dF += (F_exact[i]-F_apprx[i])*(F_exact[i]-F_apprx[i]);//this is a dot product
	}
	return dF;
}

/*double calc_dt( planet BD[], double alpha ){

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
}*/

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

	if( nParam != 8 ) {printf("\n!Usage: ./<this> <N> <dt> <frameskip> <initial velocity> <box width> <initial conditions width> <alpha>\n\n"); exit(0);}

	N          = atoi(paramList[1]);
	dt         = atof(paramList[2]);
	int fSkip  = atoi(paramList[3]);
	double v   = atof(paramList[4]);
	L          = atof(paramList[5]);
	double scatter  = atof(paramList[6]);
	alpha      = atof(paramList[7]);

	planet BD[N];
	vector posFirst[N];
	vector velFirst[N];
	vector posExact[N];
	vector velExact[N];
	vector posInitial[N];
	vector velInitial[N];
	double fError;
	vector dif;
	vector vcm;
	int i,j;
	totalIndeces = 0;
	int frame = 0;
	double MASS = 0;
	double Ei,Pi;
	double delta_F;//stores total difference in force^2
	double ratio_times;

	long int force_t;
	long int timestep;
	long int exact_timestep;
	long int pop_regions_t;
	long int resetting_t;
	long int pop_list_t;
	long int collision_t;
	long int omelyan_t;
	long int anim_t;

	for(i=0; i<LVL; i++)
	{
		totalIndeces += pow(8,i);
	}
	regions=(region*)malloc(sizeof(region)*totalIndeces);//FIXME Is this correct?
//	printf("!Number of planets: %d\n", N);
//	printf("!Frameskip is %d\n",fSkip);
//	printf("!Total # of indeces in region array: %d\n", totalIndeces);

	findKids(LVL);

	vcm.x = 0;
	vcm.y = 0;
	vcm.z = 0;

	double t,theta,rho=20;
	com.x = 0;//to check the COM finding of my function
	com.y = 0;
	com.z = 0;



	for( i=0; i<N; i++){

			BD[i].num = i;

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

//	printf("!Finished initializing planet locations\n");

	for( i=0; i<N; i++ ){
			vcm = vcm + BD[i].m*BD[i].vel;
	}

	com = com/MASS;
//	printf("!The center of mass is (%f,%f,%f)\n",com.x,com.y,com.z);

	vcm = vcm/MASS;

	for( i=0; i<N; i++ ){
			BD[i].vel = BD[i].vel - vcm;
	}

	Ei = energy(BD);
	Pi = mom(BD);
//	printf("!Pi = %f, Ei = %f\n",Pi,Ei);


	/*
	   FILE * E_vs_t = fopen("E_vs_t_alpha_900.txt","w");
	   FILE * E_vs_t_c = fopen("E_vs_t_alpha_900_c.txt","w");//if there is a collision, print to another file
	   FILE * P_vs_t = fopen("P_vs_t_alpha_900.txt","w");
	   FILE * P_vs_t_c = fopen("P_vs_t_alpha_900_c.txt","w");
	 */
	FILE * fE_dt = fopen("fE_dt.txt","w");
	for(i=0; i<N; i++)
	{
		posFirst[i]=BD[i].pos;
	}
	for(dt=1e-2; dt>1e-8; dt/= 2)
	{
		//printf("dt: %e\n", dt);
		fError=0;

		for(i=0; i<N; i++)//Store initial positions so we can reset after 100 small Omelyan steps
		{
			posInitial[i]=BD[i].pos;
			velInitial[i]=BD[i].vel;
			//printf("posInitial: %e %e %e\n", posInitial[i].x, posInitial[i].y, posInitial[i].z);
			//printf("BD: %e %e %e\n", BD[i].pos.x, BD[i].pos.y, BD[i].pos.z);
		}
		//Take 100 small steps
		for( t=0; t<1; t+=.01 )
		{
			reset_region(regions[0]);
			populateRegions(BD, totalIndeces);//find region properties
			for(j=0; j<totalIndeces; j++)
			{
				clearList(&(regions[j].planets));
			}
			popLists(BD);//populate lists
			omelyan(BD);
			Force(BD);
		}	
		for(i=0; i<N; i++)//Store exact positions after 1000 small Omelyan steps
		{
			posExact[i]=BD[i].pos;
			velExact[i]=BD[i].vel;
		}

		for(i=0; i<N; i++)//Reset positions to before 1000 steps
		{
			BD[i].pos=posInitial[i];
			BD[i].vel=velInitial[i];
		}

		//Then take one large Omelyan step
		dt*=100;
		reset_region(regions[0]);	
		populateRegions(BD, totalIndeces);
		for(i=0; i<totalIndeces; i++)
		{
			clearList(&(regions[i].planets));
		}
		popLists(BD);
		omelyan(BD);
		Force(BD);
		dt/=100;

		//Now add error to sum, to be averaged later
		for(i=0; i<N; i++)
		{
			//	printf("exact: %e %e %e\n", posExact[i].x, posExact[i].y, posExact[i].z);
			//	printf("approx: %e %e %e\n", BD[i].pos.x, BD[i].pos.y, BD[i].pos.z);
			dif=velExact[i]-BD[i].vel;
			//printf("sep: x=%e y=%e z=%e --- sep*sep = %e\n", sep.x, sep.y, sep.z, sep*sep);
			fError+=sqrt(dif*dif);
			//printf("sep: %e %e %e -- fError: %e\n", sep.x, sep.y, sep.z, fError);
			//	printf("fError %e\n", fError);
		}

		fError/=100;//Divide fError by the number of timesteps taken to get average error over run time
		fprintf(fE_dt,  "----%e %e----\n", dt, fError);
		printf("%e %e\n", dt, fError);
		for(i=0; i<N; i++)
		{
			BD[i].vel=velFirst[i];
			BD[i].pos=posFirst[i];
		}
	}
	//fprintf(dF_vs_alpha,"%e %e\n", alpha, delta_F);
	//}
//	fclose(dE_dt);
//	printf("Im done\n");
}


