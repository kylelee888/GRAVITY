/*

To compile:
g++ <filename.cpp> -fopenmp -o <filename>

To run:
./<filename> <size> <iterations> <parallel>

*/
#include <omp.h>
#include <stdint.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int nParam, char **paramList)
{
	if(nParam < 3) {printf("Usage: ./<this> <siz>e <iterations> <parallel>\n");}

	const int size = atoi(paramList[1]);
	const int its = atoi(paramList[2]);
	int parallel = atoi(paramList[3]);
	int i, j, x, p;

	struct timespec start, end;
	uint64_t dtTot, dtAvg;

	for(i=0; i<its; i++)
	{
		clock_gettime(CLOCK_MONOTONIC_RAW, &start);
		#pragma omp parallel for 
		for(int n=0; n<size; n++)
		{
			x=0;
			for(j=0; j<size*1000; j++)
				x*=pow(1,1.4142);

		//	printf("%d\n", n);
		}
		clock_gettime(CLOCK_MONOTONIC_RAW, &end);
		dtTot += (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
	}
	dtAvg=dtTot/its;
	printf("Parallel, no scheduling: %lu\n", dtAvg);

	dtTot=0;
	for(i=0; i<its; i++)
	{
		clock_gettime(CLOCK_MONOTONIC_RAW, &start);
		#pragma omp parallel for schedule(static)
		for(int n=0; n<size; n++)
		{
			x=0;
			for(j=0; j<size*1000; j++)
				x*=pow(1,1.4142);

		//	printf("%d\n", n);
		}
		clock_gettime(CLOCK_MONOTONIC_RAW, &end);
		dtTot += (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
	}
	dtAvg=dtTot/its;
	printf("Parallel, static: %lu\n", dtAvg);

	dtTot=0;
	
	for(i=0; i<its; i++)
	{
		clock_gettime(CLOCK_MONOTONIC_RAW, &start);
		#pragma omp parallel for schedule(dynamic) if(parallel)
		for(int n=0; n<size; n++)
		{
			p = omp_in_parallel();
			x=0;
			for(j=0; j<size*1000; j++)
				x*=pow(1,1.4142);

		//	printf("%d\n", n);
		}
	
		printf("parallel: %d -- p: %d\n", parallel, p);
		clock_gettime(CLOCK_MONOTONIC_RAW, &end);
		dtTot += (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
	}
	dtAvg=dtTot/its;
	printf("Parallel, dynamic: %lu\n", dtAvg);

}










