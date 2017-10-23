#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[])
{

	char var[15];
	char val[15];
	
	int N, frameskip, init_type, parallel, print;
	double dt, v_initial, box_width, initCond_width, alpha, tmax;

	FILE *in=fopen("config.txt","r");
printf("While\n");
	while (fscanf(in,"%s %s",var,val) != EOF)
	{
		printf("var: %s   val: %s\n", var, val);
		if (strcmp(var,"N") == 0) N=atoi(val);
		if (strcmp(var,"dt") == 0) dt=atof(val);
		if (strcmp(var,"frameskip") == 0) frameskip=atof(val);
		/*if (strcmp(var,"v_initial") == 0) v_initial=atoi(val);
		if (strcmp(var,"box_width") == 0) box_width=atoi(val);
		if (strcmp(var,"initCond_width") == 0) initCond_width=atoi(val);
		if (strcmp(var,"alpha") == 0) alpha=atoi(val);
		if (strcmp(var,"tmax") == 0) tmax=atoi(val);
		if (strcmp(var,"init_type") == 0) init_type=atoi(val);
		if (strcmp(var,"parallel") == 0) parallel=atoi(val);
		if (strcmp(var,"print") == 0) print=atoi(val);
		*/
	}
	
	printf("N=%d \ndt=%f \nframeskip=%d \n", N, dt, frameskip);
}
