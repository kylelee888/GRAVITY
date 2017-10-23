typedef struct
{
  double x;
  double y;
  double z;
}vector;

vector operator+(vector a, vector b);
double operator*(vector a, vector b);
vector operator^(vector a, vector b);
vector operator*(double a, vector b);
vector operator*(vector a, double b);
vector operator-(vector a, vector b);
vector operator/(vector a, double b);
double radius( vector BD1, vector BD2 );
double radNew( vector a, vector b);
double vecLength(vector in);
