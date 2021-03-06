//rElem : a naive way to check for collisions. It would check all planets against all other planets and see if any were touching. New way is to check if any planets are close during the force calculation. 



double rElem( planet BD[] ){

  int i,j;
  double rad;
  double rho = 20;
  int other;
  pln *otherPln;
  int regionIndex;
  region regionOfInterest;
  int neighbors[26];//there are 3^3 regions in a 3x3 cube, -1 for the center
  collision = 0;
  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++ ){
      rad = radius(BD[i].pos,BD[j].pos);
      
      if( rad < 0.8*(BD[i].r+BD[j].r) ){
	
	BD[i].pos = (BD[i].m*BD[i].pos + BD[j].m*BD[j].pos)/(BD[i].m + BD[j].m);//find new properties 
	BD[i].vel = (BD[i].m*BD[i].vel + BD[j].m*BD[j].vel)/(BD[i].m + BD[j].m);
	BD[i].m   = (BD[i].m + BD[j].m);
	BD[i].r   = pow( BD[i].m/rho, 1.0/3 );

	BD[j].num = BD[N-1].num;
	BD[j].pos = BD[N-1].pos;//move information
	BD[j].vel = BD[N-1].vel;
	BD[j].m   = BD[N-1].m;
	BD[j].r   = BD[N-1].r;
	
	N--;
	a++;
	collision = 1;
      }
    }
  }

  //looping through all regions is ridiculous
  //this was an attempt to simply look at all planets in
  //the same box as the planet of interest at the lowest level
  //still not as fast as the method listed above
  /*for(i=0;i<N;i++){
    regionIndex = BD[i].level[LVL-1];
    regionOfInterest = regions[regionIndex];//region our planet is in
    //findNeighbors(regionIndex,neighbors);
    if(regionOfInterest.numPln > 1){//if there is another planet to collide with
      otherPln = regionOfInterest.planets;//point to the list of planets
      while(otherPln != NULL){//loop through all planets in that regions list
	other = otherPln->plnNum;
	if(i != other){//other is not BD[i]
	  rad = radius(BD[i].pos,BD[other].pos);
	  if(rad < 0.8*(BD[i].r + BD[other].r)){
	      //find new properties of collided planets
	      BD[i].pos = (BD[i].m*BD[i].pos + BD[other].m*BD[other].pos)/(BD[i].m + BD[other].m); 
	      BD[i].vel = (BD[i].m*BD[i].vel + BD[other].m*BD[other].vel)/(BD[i].m + BD[other].m);
	      BD[i].m   = (BD[i].m + BD[other].m);
	      BD[i].r   = pow( BD[i].m/rho, 1.0/3 );

	      //shift planets
	      BD[other].pos = BD[N-1].pos;
	      BD[other].vel = BD[N-1].vel;
	      BD[other].m   = BD[N-1].m;
	      BD[other].r   = BD[N-1].r;

	      N--;
	      a++;
	  }
	}
	otherPln = otherPln -> nextPln;
      }
    }
  }*/
