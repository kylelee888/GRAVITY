#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"

int FindRegion(int i,int j,int k,int N) //i is columns, j is rows, k is depth, N is what level we are in
{
  int sum;
  int m;
  for (m=0;m<=N;m=m+1)
  {
      sum += pow(8,m);
  }
  return sum + i + j*pow(2,N) + k*pow((pow(2,N)),2);
}

//Only works if each box is divided into eighths
