#include <time.h>

long int timer[10];

void starttimer(int t)
{
  timer[t]=clock();
}

long int stoptimer(int t)
{
  return clock()-timer[t];
}

bool istime(int delay)
{
  static int last=0;
  if (clock() - last > delay*1000)
  {
    last=clock();
    return true;
  }
  else return false;
}
