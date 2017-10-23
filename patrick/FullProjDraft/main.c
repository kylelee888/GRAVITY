//External packages
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//Functions and classes we made:
#include "vector.h"//Keep classes separate b/c vector math we define makes file messy
#include "planet.h"
#include "region.h"
#include "getArrayIndex.h"
#include "getTotalMass.h"
#include "getRegion.h"
#include "addPlanet.h"

#define THETA 1.35120719195966
//UNITS
//Mass	 : Earth Mass
//Length : AU (astronomical unit)
//Time	 : YEAR
//fix planet sizes, way too big

// link for Forest-Ruth Algorithm 4th order
//http://young.physics.ucsc.edu/115/leapfrog.pdf


//GLOBAL variable declaration
int LVL = 6;//Lowest level of subdivision
int N; //number of planets		
double dt; //time step
double G = 4*M_PI*M_PI; //newtons constant 
int a; //number of collisions
double L; //inital grid size
vector com; //center of mass to test my functions
//double THETA = 1/(2-pow(2,1./3));
int blah, totalIndeces;
for(blah=0; blah<LVL; blah++)
{
	totalIndeces+=pow(8,blah);
} printf("Total # of indeces in region array: %d\n", totalIndeces);
region regions[totalIndeces];


//=======================

	//mass() -> getTotalMass(), included via functions.h

//=======================

	//arrayIndex() -> getArrayIndex(), included via functions.h

//=======================

	//planetRegion() -> populateRegions(), included via populateRegions.h

//=======================
double radius( vector BD1, vector BD2 ){
	return sqrt( (BD1 - BD2) * (BD1 - BD2) );
}

//=======================
vector posHalfStep( planet BD[] ){
	
	int i;
	
	for( i=0; i<N; i++ ){
		BD[i].pos = BD[i].pos + BD[i].vel*dt/2;
	}
}

//=======================
vector velFullStep( vector a[], planet BD[] ){
	
	int i;
	
	for( i=0; i<N; i++ ){
		BD[i].vel = BD[i].vel + a[i]*dt;
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
vector vel4_FWD( vector a[], planet BD[] ){

	int i;

	for( i=0; i<N; i++ ){
		BD[i].vel = BD[i].vel + THETA*a[i]*dt;
	}
}

//=======================
vector vel4_BCK( vector a[], planet BD[] ){

	int i;

	for( i=0; i<N; i++ ){
		BD[i].vel = BD[i].vel + (1-2*THETA)*a[i]*dt;
	}
}

//=======================
vector Force( planet BD[], vector a[] ){

	int i,j;
	double rad;
	vector F;
	
	for( i=0; i<N; i++ ){
		for( j=i+1; j<N; j++){//the i+1 is what makes it faster instead of j=0 and going through to N

			rad = radius(BD[i].pos,BD[j].pos);
			F	 = G * BD[j].m * BD[i].m * (BD[j].pos-BD[i].pos) / (rad*rad*rad);
			a[i] = a[i] + F/BD[i].m;
			a[j] = a[j] - F/BD[j].m;
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
	BD[i].m	 = (BD[i].m + BD[j].m);
	BD[i].r	 = pow( BD[i].m/rho, 1.0/3 );
	
	BD[j].pos = BD[N-1].pos;//move information
	BD[j].vel = BD[N-1].vel;
	BD[j].m	 = BD[N-1].m;
	BD[j].r	 = BD[N-1].r;
		
	N--;
	a++;
			}
		}
	}
}
	
//=======================

double leap( planet BD[] ){
	
	vector a[N];
	int i;

	for( i=0; i<N; i++ ){
		a[i].x = 0;
		a[i].y = 0;
		a[i].z = 0;
	}

	//	planetRegion(BD);//finds the regions a planet is in

	//regionProp(BD);//finds the number of planets, total mass, and COM of every box in every level

	//pos4_FWD(BD);
	
	// Force(BD,a);

	//vel4_FWD(a,BD);

	//pos4_BCK(BD);

	//Force(BD,a);

	//vel4_BCK(a,BD);

	//pos4_BCK(BD);

	//Force(BD,a);

	//vel4_FWD(a,BD);

	//pos4_FWD(BD);

	//rElem(BD);
	
	posHalfStep(BD);

	Force(BD,a);
		
	velFullStep(a,BD);
		
	posHalfStep(BD);
	
	rElem(BD);
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

	N = atoi(paramList[1]);
	dt = atof(paramList[2]);
	int fSkip = atoi(paramList[3]);
	double v = atof(paramList[4]);
	L = atof(paramList[5]);

	planet BD[N];//FIXME Declare globally? Can we without specifying N?
	vector vcm;
	int i;
	int frame = 0;
	
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

		BD[i].r	 = 0.1;
		BD[i].m	= rho*BD[i].r*BD[i].r*BD[i].r;

		com = com + BD[i].pos * BD[i].m;
	}

	for( i=0; i<N; i++ ){
		vcm = vcm + BD[i].m*BD[i].vel;
	}
	
	com = com/mass(BD);

	vcm = vcm/mass( BD );
	
	for( i=0; i<N; i++ ){
		BD[i].vel = BD[i].vel - vcm;
	}
			
	for( t=0; 1; t+=dt )
	{
		leap(BD);
		if( frame % fSkip == 0 ){
			printf("T -0.8 0.82\nE = %f : M = %f : P = %f : # Collisions = %d\n",energy(BD),mass(BD),mom(BD),a);
			for( i=0; i<N; i++){
	printf("c3 %e %e %e %e\n",BD[i].pos.x,BD[i].pos.y,BD[i].pos.z,BD[i].r);
			}
			printf("F\n");
		}
		frame++;
	}
}
