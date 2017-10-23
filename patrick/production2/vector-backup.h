typedef struct
{
  double x;
  double y;
  double z;
} vector;

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


  


  
