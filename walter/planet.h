//#include "vector.h"


typedef struct
{
  vector pos;
  vector vel;
  vector acc;
  double m;
  double r;
  int level[10];//for whatever reason you can't put LVL here
               // why is this here, incidentally? -W
  int num;
} planet;

