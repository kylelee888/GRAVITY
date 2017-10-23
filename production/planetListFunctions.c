#include <math.h>
#include "pln.h"
#include "vector.h"
#include "planet.h"
#include "region.h"
#include <stdlib.h>
extern int N, LVL;
extern region *regions;//FIXME Include region.h?

void add(int num, pln **list)
{
	pln *newPln;//"Make a pointer to a pln, call it 'newPln'"
	newPln=(pln *)malloc(sizeof(pln));//"Allocate a pln's worth of memory to the space pointed to by 'newPln'"
	newPln->plnNum = num;//Follow 'newPln' to the pln it points to, set plnNum of that pln to the value num
	newPln->next = *list;//Follow 'newPln' to the pln it points to, set next of that pln to the pointer to the beginning of the old list
	*list = newPln;//Set this pointer to the pointer 'newPln'
}

void clearList(pln **list)
{
	pln *curr = *list;
	pln *next;
	while(curr!=NULL)
	{
		next = curr->next;
		free(curr);
		curr = next;
	}
}

void popLists(planet BD[])
{
	
 // FIXME: Patrick still revising this
	int i,j;
	for(i=0; i<N; i++)
	{
		for(j=1; j<LVL; j++)
		{
			add(i, &(regions[BD[i].level[j]].planets));
		}	
	}



}
