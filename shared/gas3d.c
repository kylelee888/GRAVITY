#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
int first=1;
int N=80;
double pad=2;
int fps=60;
double dt=1e-4;
double step=1;
double vinit=0.3;
double tequilib=3;
double damp=0;
double fconst=1;
double r0=0.015;
long int timer[10];
double g=0.00;
int NG;

struct list_el {
  int val;
  struct list_el * next;
};

typedef struct list_el item;

double hypot3(double x, double y, double z)
{
  return sqrt(x*x+y*y+z*z);
}

double green(double F)
{
  static double v;
  v=(log(F)+9)*0.4; 
  if (v>0.4) return 0.4; else return v;
}

double blue(double F)
{
  static double v;
  v=(log(-F)+10)*0.4; 
  if (v>0.8) return 0.8; else return v;
}

int getzone(double x)
{
  int z;
  z=(int)((x+1)*0.5*NG);
  if (z>=NG) z=NG-1;
  if (z<0) z=0;
  return z;
}

void add(item **list, int v)
{
  //   printf("!Adding item %d\n",v);
  item *curr;
  curr=(item *)malloc(sizeof(item));
  curr->val = v;
  curr->next = *list;
  *list=curr;
}
void del (item **list, int v)
{
  item *curr;
  curr=*list;
  if (curr == NULL) {return;}
  //  printf("!Deleting item %d... list starts with %d (%d)\n",v,*list,(*list)->val);
  if (curr->val == v) {*list=curr->next; free(curr);}
  else
    while (curr -> next != NULL)
    {
      if (curr->next->val == v)
      {
	free(curr->next); curr->next=curr->next->next;
	break;
      }
      curr=curr->next;
    }

}

int check_zonelist(item *zonelist[NG][NG][NG],int xzone[N],int yzone[N],int zzone[N], double x[N],double y[N], double z[N])
{
  //    printf("!Start update\n");
  int i,j,xx,yy,zz,n=0,aaa=0;


  for (i=0;i<N;i++)
  {
    xx = getzone(x[i]);
    yy = getzone(y[i]);
    zz = getzone(z[i]);
    if (xzone[i] != xx || yzone[i] != yy || zzone[i] != zz)
    {
      //          printf("!Ball %d is in box %d,%d, should be in %d,%d\n",i,xzone[i],yzone[i],xx,yy);
      if (xzone[i] >= 0 && xzone[i] < NG && yzone[i] >= 0 && yzone[i] < NG && (zzone[i] >= 0 && zzone[i] < NG)) del(&zonelist[xzone[i]][yzone[i]][zzone[i]],i);
      add(&zonelist[xx][yy][zz],i);
    }
    xzone[i]=xx; yzone[i]=yy; zzone[i]=zz;

  }

  //  printf("!End update\n");
}

void print_zonelist(item *zonelist[NG][NG][NG],int x, int y, int z)
{
  item *curr;
  curr=zonelist[x][y][z];
  printf("!Zone list for %d,%d,%d: ",x,y,z);
  while (curr)
  {
    printf("%d ",curr->val);
    curr=curr->next;
  }
  printf("\n");
}

double V (double r)
{
  if (r > r0*pad) return V(r0*pad);
  r/=r0;
  static double r6;
  r6=r*r*r*r*r*r;
  return (4*r0*fconst*(1/(r6*r6)-1/r6));
}

void starttimer(int t)
{
  timer[t]=clock();
}

long int stoptimer(int t)
{
  return clock()-timer[t];
}

// ACTUALLY RETURNS F/r
double force (double r2)
{
  if (r2>(r0*r0*pad*pad)) return 0;
  static double f;
  r2=r2/(r0*r0);
  static double r8,r14,r6;
  r6=r2*r2*r2;
  r8=r6*r2;
  r14=r6*r8;
  f= -24 * (2/r14 - 1/r8) * fconst / r0;
  // f= -24 * (2/pow(r,13) - 1/pow(r,7)) * fconst;
  //  printf ("!force between particles at radius %e is %e\n",r,f);

  return f;
}

// f=0 when 2/r13 = 1/r7, or 2 = (r/r0)^6, r = r0 * sixth root of 2

double drnd48(void)
{
  return drand48()-0.5;
}

int main(int argc, char **argv)
{
  FILE *fp=fopen("PVNkT.ax","w");
  double thermo_interval=10,pct;
  double Ttarget=10;
  if (argc>1) N=atoi(argv[1]);
  if (argc>2) dt=atof(argv[2]);
  if (argc>3) vinit=atof(argv[3]);
  if (argc>4) pad=atof(argv[4]);
  if (argc>5) thermo_interval=atof(argv[5]);
  if (argc>6) Ttarget=atof(argv[6]); 
  printf("!Usage: <this> <N> <dt> <vinit> <pad> <thermo_interval> <T target>\n");
  printf("!Parameters read: n=%d dt=%e vinit=%e pad=%e thermo interval = %e\n",N,dt,vinit,pad,thermo_interval);
  double next_thermo=thermo_interval;
  double color;
  NG=2.0/(r0*pad);
  // NG=1;
  item *zonelist[NG][NG][NG];
  item *curr;
  int xzone[N],yzone[N],zzone[N];
  long int lastframe=0;
  int gone[N],mm;
  double x[N],y[N],z[N],vx[N],vy[N],vz[N],m[N];
  double r2;
  double xh[N],yh[N],zh[N],vxh[N],vyh[N],vzh[N],R[N],drawr[N];
  double vmult,Tframe;
  int bloc;
  int zx,zy,zz,k,l;
  int i,j,step=0;
  double F,r,t;
  double xs,ys,zs;
  double T,U;
  double interval;
  double vf=3.5,vxcom,vycom,vzcom;
  int frameskip=0;
  srand48(232444);
  int ii,jj,inter;
  int update,drawn=0;
  int time_mundane=0,time_rk2=0,time_anim=0,time_energy=0,time_zonelist=0;
  double P=0;
  double Precord=0,Trecord=0,Taccum=0;
  FILE *Tlog, *Ulog, *Elog;
  int Tsamples;
  for (i=0;i<NG;i++)   for (j=0;j<NG;j++) for (k=0;k<NG;k++)
  {
    zonelist[i][j][k]=NULL;
  }


  Tlog=fopen("Tlog","w");
  Ulog=fopen("Ulog","w");
  Elog=fopen("Elog","w");

  fprintf(Tlog,"#cm 1\n");
  fprintf(Ulog,"#cm 2\n");
  fprintf(Elog,"#cm 3\n");
  fflush(Tlog);
  fflush(Ulog);
  fflush(Elog);
  j=k=l=0;  
  interval = r0 * pow(2,1./6.);
  for (i=0;i<N;i++)
  {
    j++;
    if (j>pow(N,1./3.)) {j=0;k++;}
    if (k>pow(N,1./3.)) {k=0;l++;}

    xzone[i]=-1; yzone[i]=-1; zzone[i]=-1;
    x[i]=-0+k*interval;
    y[i]=-0+l*interval;
    z[i]=0+j*interval;;
    vx[i]=drnd48()*vinit;
    vy[i]=drnd48()*vinit;
    vz[i]=drnd48()*vinit;
    m[i]=1;
  }
  // m[0]=1000;
  check_zonelist(zonelist,xzone,yzone,zzone,x,y,z);
  //  vx[0]=vy[0]=0;
  double U0,There,Uhere;
  printf("!start main loop\n");
  int time_keep=clock();
  for (t=0;1;t+=dt)
  {
    if (clock() > time_keep + 1e6)
    {
      printf("!%d usec total: RK2 %d, animation %d, zonecheck %d, energy %d, mundane %d\n",(int)(clock()-time_keep),time_rk2,time_anim,time_zonelist,time_energy,time_mundane);
      printf("!RK2: %d interactions, %d ns per interaction\n",inter,(time_rk2*1000)/(inter));
      time_rk2=time_anim=time_zonelist=time_energy=time_mundane=inter=0;
      time_keep=clock();
   }
    frameskip++;

    if (t > next_thermo)
    {
      Precord=P/(thermo_interval*24); // P accumulates total impulse delivered by walls; walls have area 6 * 4
      P=0;
      Trecord=Taccum/Tsamples/N;
      //      printf("!At time %e: Temperature: accumulator %e, samples %d\n",t,Taccum/N,Tsamples);
      Taccum=0;
      Tsamples=0;
      next_thermo+=thermo_interval;
      fprintf(fp,"%e %e\n",t,Precord*8.0/(Trecord*N));
      fflush(fp);
    }
    update=0;
    if (clock() > lastframe + CLOCKS_PER_SEC/fps) // time to draw
    {
      update=1;
      drawn++; 
      if (drawn % 100 == 1)
      {
	starttimer(1);
	update=1;
	T=U=0;
	for (i=0;i<N;i++)
	{
	  T+=0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
	  for (j=i+1;j<N;j++)
	  {
	    r=hypot3(x[i]-x[j],y[i]-y[j],z[i]-z[j]);
	    U+=V(r);
	  }
	}
	            time_energy += stoptimer(1); 
        printf("!Energy calculated in %d usec\n",time_energy);
      }
    starttimer(3);
      printf("C 1 1 1\n");
      for (i=0;i<N;i++)
      {
	//        if (i==0) printf("C 1 1 1\n"); else printf("C 0.6 0.6 0.6\n");
	printf("c3 %e %e %e %e\n",x[i],y[i],z[i],r0/8);
      }
      printf("C 1 1 1\nT -0.9 -0.85\nt=%.3f, %d frames   %05d FPS    E=%.8e\n",t,step,frameskip*fps,T+U);
      printf("T -0.9 0.9\nkT=%.3e   P=%.3e   PV=%.3e   NkT=%.3e\n",Trecord,Precord,Precord*8,Trecord*N);
      printf("T -0.9 0.95\nKE accum=%.3e   J accum=%.3e   %.1f%% done\n",Taccum,P,100-(next_thermo-t)/(thermo_interval)*100);
      printf("T -0.9 0.85\nTarget = %e\tTemp = %e\tMult = %e\n",Ttarget,Tframe,vmult);
      double barx;
      //    for (barx=-0.9;barx<-.9+0.3*(1-(next_thermo-t)/thermo_interval); barx+=0.01) printf("l %e 1.08 %e 1.10\n",barx,barx);
      frameskip=0;
      printf("C 0.7 0.2 0.2\n"); 
  //      for (i=0;i<=NG;i+=NG)
  //      {
  //        double v;
  //        v=-1+(double)i/NG*2;
  //        printf("l3 -1 -1 %e -1  1 %e\n",v,v);
  //        printf("l3  1 -1 %e  1  1 %e\n",v,v);
  //     }

 
    printf("l3 -1 -1 -1 -1 -1  1\n");
    printf("l3 -1 -1  1 -1  1  1\n");
    printf("l3 -1  1  1 -1  1 -1\n");
    printf("l3 -1  1 -1 -1 -1 -1\n");
    printf("l3  1 -1 -1  1 -1  1\n");
    printf("l3  1 -1  1  1  1  1\n");
    printf("l3  1  1  1  1  1 -1\n");
    printf("l3  1  1 -1  1 -1 -1\n");
    printf("l3  1 -1 -1 -1 -1 -1\n");
    printf("l3  1  1 -1 -1  1 -1\n");
    printf("l3  1 -1  1 -1 -1  1\n");
    printf("l3  1  1  1 -1  1  1\n");
  
    pct=(next_thermo-t)/thermo_interval;
    printf("C 0.2 0.7 0.2\n");
    printf("l3 -%e -%e -%e -%e -%e  %e\n",1.0,1.0,pct,1.0,1.0,pct);
    printf("l3 -%e -%e  %e -%e  %e  %e\n",1.0,pct,1.0,1.0,pct,1.0);
    printf("l3 -%e  %e  %e -%e  %e -%e\n",1.0,1.0,pct,1.0,1.0,pct);
    printf("l3 -%e  %e -%e -%e -%e -%e\n",1.0,pct,1.0,1.0,pct,1.0);
    printf("l3  %e -%e -%e  %e -%e  %e\n",1.0,1.0,pct,1.0,1.0,pct);
    printf("l3  %e -%e  %e  %e  %e  %e\n",1.0,pct,1.0,1.0,pct,1.0);
    printf("l3  %e  %e  %e  %e  %e -%e\n",1.0,1.0,pct,1.0,1.0,pct);
    printf("l3  %e  %e -%e  %e -%e -%e\n",1.0,pct,1.0,1.0,pct,1.0);
    printf("l3  %e -%e -%e -%e -%e -%e\n",pct,1.0,1.0,pct,1.0,1.0);
    printf("l3  %e  %e -%e -%e  %e -%e\n",pct,1.0,1.0,pct,1.0,1.0);
    printf("l3  %e -%e  %e -%e -%e  %e\n",pct,1.0,1.0,pct,1.0,1.0);
    printf("l3  %e  %e  %e -%e  %e  %e\n",pct,1.0,1.0,pct,1.0,1.0);
    printf("F\n");
    
    lastframe = clock();
    time_anim += stoptimer(3);
    }
    starttimer(4);
    Tframe=0;
    for (i=0;i<N;i++)
    {
      Tframe+=0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]) * (2./3.); //0.25 from 2 dof
    }
    Taccum+=Tframe;

    Tframe/=N;
    vmult=exp((Ttarget/Tframe-1)*dt * .1);

    Tsamples++;

    for (i=0;i<N;i++)
    {
      vx[i]*=vmult;
      vy[i]*=vmult;
      vz[i]*=vmult;
    }
 
    

    step++;
    time_mundane+=stoptimer(4);
    starttimer(1);
    // start force loop: RK2 halfstep
    check_zonelist(zonelist,xzone,yzone,zzone,x,y,z);
    time_zonelist+=stoptimer(1);
    starttimer(2);
    for (i=0;i<N;i++)
    { 
      first=0;
      // printf("rk2 half: looking at particle %d\n",i);
      x[i] = x[i]+vx[i]*dt*0.5;
      y[i] = y[i]+vy[i]*dt*0.5;
      z[i] = z[i]+vz[i]*dt*0.5;
    }
    
for (i=0;i<N;i++)
    { 
      // loop over all adjacent zones to the zone i is in, and look for interactions
      zx=getzone(x[i]);
      zy=getzone(y[i]);
      zz=getzone(z[i]);

      for (k=zx-1;k<=zx+1;k++)       for (l=zy-1;l<=zy+1;l++)   for (mm=zz-1;mm<=zz+1;mm++)
      {
	// printf("considering interactions with stuff in zone %d,%d\n",k,l);
	if (k<0 || k>=NG || l<0 || l>=NG || mm<0 || mm>=NG) continue;
	curr=zonelist[k][l][mm];
	while (curr)
	{
	  j=curr->val; if (i<=j) {curr=curr->next; continue;}
	  inter++;
	  //     printf("!Interaction between %d and %d (halfstep)\n",i,j);
	  //          printf("C 0 1 0\nl %f %f %f %f\n",x[i],y[i],x[j],y[j]);
	  xs=x[i]-x[j];
	  ys=y[i]-y[j];
	  zs=z[i]-z[j];
	  r2=xs*xs+ys*ys+zs*zs;
	  F=-force(r2)*dt; // NOT ACTUALLY FORCE, equal to -f(r^2)/r*dt
	  vx[i] += F * (xs)/m[i];
	  vx[j] -= F * (xs)/m[j];
	  vy[i] += F * (ys)/m[i] - g*dt;
	  vy[j] -= F * (ys)/m[j];
	  vz[i] += F * (zs)/m[i];
	  vz[j] -= F * (zs)/m[j];
	  curr=curr->next;
	}

      }

    }
  for (i=0;i<N;i++)
    { 
      

      first=0;
      // printf("rk2 half: looking at particle %d\n",i);
      x[i] = x[i]+vx[i]*dt*0.5;
      y[i] = y[i]+vy[i]*dt*0.5;
      z[i] = z[i]+vz[i]*dt*0.5;
    


      // check for hitting walls, conserves energy exactly so don't need rk2 here
      if (x[i]>1 && vx[i]>0)
      {
	P+=vx[i]*m[i]*2;
	vx[i] *= -1;
      }
      if (x[i]<-1 && vx[i]<0)
      {
	P-=vx[i]*m[i]*2;
	vx[i] *= -1;
      }
      if (y[i]>1 && vy[i]>0)
      {
	P+=vy[i]*m[i]*2;
	vy[i] *= -1;
      }
      if (y[i]<-1 && vy[i]<0)
      {
	P-=vy[i]*m[i]*2;
	vy[i] *= -1;
      }
      if (z[i]>1 && vz[i]>0)
      {
	P+=vz[i]*m[i]*2;
	vz[i] *= -1;
      }
      if (z[i]<-1 && vz[i]<0)
      {
	P-=vz[i]*m[i]*2;
	vz[i] *= -1;
      }
    } // over particle i
    time_rk2+=stoptimer(2); 
  } // over t
}
