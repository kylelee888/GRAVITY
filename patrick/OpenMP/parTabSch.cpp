/*

To compile:

g++ <filename.cpp> -fopenmp -o <filename>

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int nParam, char **paramList)
{
	if(nParam < 3) {printf("Usage: ./<this> size \n");}

	const int size = atoi(paramList[1]);
	int i, y;
	#pragma omp parallel for schedule(dynamic)
	for(int n=0; n<size; n++)
	{
		y=0;
		for(i=0; i<size*1000; i++)
			y*=pow(1,1.4142);

		printf("%d\n", n);
	}

}
