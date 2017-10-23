//#include "vector.h"

typedef struct
{
  vector pos;
  vector vel;
  vector acc;
  vector col;
  double m;
  double r;
  int num;
  int level[6];//for whatever reason you can't put LVL here
               // why is this here, incidentally? -W
} planet;

