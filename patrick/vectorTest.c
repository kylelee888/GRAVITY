#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"

int main(void)
{
	vector x, y;
	x.x=1; x.y=1; x.z=1;
	y.x=1; y.y=1; y.z=1;
	if(x==y)
	{
		printf("Yup\n");
	}
	else
	{
		printf("Nope\n");
	}

}
