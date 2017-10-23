#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int nParam, char **paramList)
{
	if(nParam < 2) {printf("Usage: ./<this> size\n");}

	const unsigned int size = atoi(paramList[1]);
	printf("size: %d\n", size);
	double sinTable[size];

	for(int n=0; n<size; n++)
		sinTable[n] = sin(2 * M_PI * n / size);

}
