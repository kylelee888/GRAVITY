#include "vector.h"
#include <math.h>

vector operator+(vector a, vector b)
{
  static vector c;
  c.x=a.x+b.x;
  c.y=a.y+b.y;
  c.z=a.z+b.z;
  return c;
}

double operator*(vector a, vector b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

vector operator^(vector a, vector b)
{
  static vector c;
  c.x=a.y*b.z - a.z*b.y;
  c.y=a.z*b.x - a.x*b.z;
  c.z=a.x*b.y - a.y*b.x;
  return c;
}

vector operator*(double a, vector b)
{
  static vector c;
  c.x = a*b.x;
  c.y = a*b.y;
  c.z = a*b.z;
  return c;
}

vector operator*(vector a, double b)
{
  static vector c;
  c.x=a.x*b;
  c.y=a.y*b;
  c.z=a.z*b;
  return c;
}

vector operator-(vector a, vector b)
{
  static vector c;
  c.x = a.x-b.x;
  c.y = a.y-b.y;
  c.z = a.z-b.z;
  return c;
}

vector operator/(vector a, double b)
{
  static vector c;
  c.x = a.x/b;
  c.y = a.y/b;
  c.z = a.z/b;
  return c;
}

vector operator/(double b,vector a)
{
    static vector c;
    c.x = b*a.x/(a.x*a.x + a.y*a.y + a.z*a.z);
    c.y = b*a.y/(a.x*a.x + a.y*a.y + a.z*a.z);
    c.z = b*a.z/(a.x*a.x + a.y*a.y + a.z*a.z);
    return c;
}

vector operator+(vector a, double b)
{
    static vector c;
    c.x = a.x + b;
    c.y = a.y + b;
    c.z = a.z + b;
    return c;
}

/*vector operator+=(vector a)
{
  this = this + a;
  return this;
}

vector operator-=(vector a)
{
  this = this - a;
  return this;
  }*/

double radius( vector BD1, vector BD2 ){
  return sqrt( (BD1 - BD2) * (BD1 - BD2) );//FIXME Does this file need to include math.h to get sqrt()?
}

  


  
