#include "vector.h"

/*
typedef struct//Here be planets
{
	vector pos;
	vector vel;
	double m;
	double r;
	int level[5];
} planet;

typedef struct
{
	int x;
	int y;
	int z;
} intVec;

typdef struct
{
	int level;
	int i,j,k;
	double mass;
	vector com;
} region; */

typedef struct
{
	int plnNum;
	struct el * next;
	//Whatever else we want to store
} el;


//Declare arrays of (pointers to el's) for each level
el *L0;
el *L1[8];
el *L2[64];
el *L3[512];
el *L4[4096];

el **regArr[] = {L0, L1, L2, L3, L4};//Declare globally -> won't need to be passed to functions -> no weird 2D array behavior. This array is necessary to make popLists more concise.

void add(int num, el **list)
{
	el *new;//"Make a pointer to an el, call it 'new'"
	new=(el *)malloc(sizeof(el));//"Allocate an el's worth of memory to the space pointed to by 'new'"
	new->plnNum = num;//Follow 'new' to the el it points to, set plnNum of that el to the value num
	new->next = *list;//Follow 'new' to the el it points to, set next of that el to pointer to the beginning of the old list
	*list = new;//Set the pointer to the beginning of the list to the pointer 'new'
}

void popLists(planet BD[])
{

	int i, j;
	
	for(i=0; i<N; i++)//For each planet
	{
		for(j=1; j<LVL; j++)//For each level
		{
			reg = BD[i].level[j];//Give me planet i's region in level j
			add(i, &regArr[j][reg]);//Add planet i to the [LVL][region] list
		}
	}
}
