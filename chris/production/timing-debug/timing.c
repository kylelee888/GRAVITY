#include <time.h>

long int timer[10];

void starttimer(int t)
{
  timer[t]=clock();
}

int stoptimer(int t)
{
  return clock()-timer[t];
}

