//#include "vector.h"

typedef struct
{
  vector pos;
  vector vel;
  vector acc;
  double m;
  double r;
  int level[5];//for whatever reason you can't put LVL here
} planet;

