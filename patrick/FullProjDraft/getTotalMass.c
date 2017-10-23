double getTotalMass( planet BD[] ){
   
  int i;
  double mass=0.0;
  for( i=0; i<N; i++ ){
    mass += BD[i].m;
  }
  return mass;
}
