#include <stdio.h>  // angle brackets: "include this thing which is part of standard C"
#include "vector.h" // quote marks: "include this thing here in this directory"
#include <math.h>

double length(vector v)
{
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

int main(void)
{
  vector A,B,C;
  printf("Enter the three components of vector A:\n");
  scanf("%lf %lf %lf",&A.x, &A.y, &A.z);
  printf("Enter the three components of vector B:\n");
  scanf("%lf %lf %lf",&B.x, &B.y, &B.z);
  
  printf("The length of vector A is %f\n",length(A));
  printf("The length of vector B is %f\n",length(B));

  C=A+B; // uses the overloaded vector addition in vector.h

  printf("A+B is equal to (%f, %f, %f)\n",C.x,C.y,C.z);
}
