void doForce(int currPlanet, int L, planet BD[])//*****Don't need BD[] if it's global?
{
	int n = pow(2,L);//n is the number of regions in one direction. n-1 is the maximum index number these regions could have
	int a, b, c;//Arbitrary indeces to be used in for loops
	intVec curr = getIJKVec(n, BD[currPlanet].level[L]);//Find which region currPlanet is in at level L, make an intVec of the coordinates of this region
	intVec temp;//Will be used to store the coordinates of the region we want to check the separation of relative to curr
	//Should calculate xl, xr, etc. here and then pass to wellSep. More efficient but potentially small potatoes.

	if(/* Num planets in this region = 1 */ || L = 5)
	{
		//do direct calculation
	}
	else
	{
		//Loop over each region at this level
		for(a=0; a<n; a++)
		{
			for(b=0; b<n; b++)
			{
				for(c=0; c<n; c++)
				{
					temp.i=a; temp.j=b; temp.k=c;
					if(wellSep(curr, temp, n) == 1 )//Order very important here
					{
						//do CoM calcs
	
					}
					else
					{
						doForce(currPlanet, L-1, BD[]);//*****Brackets necessary?
					}	
				}
			}
		}
	}
}

void doForce2(int currPlanet, int lev, planet BD[], int regionsToCheck[])//regionsToCheck will be renamed
{//currPlanet is the planet we are doing force caluclations on, lev is the current level, regionsToCheck is an array of the region numbers for the regions we are going to do force calculations for
	int size = sizeof(regionsToCheck[])/sizeof(regionsToCheck[0]);//"Take the total number of bytes occupied by regionsToCheck, divide it by the number of bytes occupied by one element, and this is the number of elements in this array. Set size to this value." FIX ME: Don't need this because all children arrays are length 8
	int a;//For loop index
	int n = pow(2,lev);//n is the number of regions along one direction. n-1 is the maximum index number these can have
	intVec curr = getIJKVec(n, currPlanet);//Put ijk-coordinates of currPlanet's region into an intVec, call it curr
	intVec temp;//To be used later to hold the ijk-coordinates of each region in regionToCheck


	if(/* Num planets in this region = 1 */ || L = 5)
	{
		//do direct calculation
	}
	else
	{
		for(a=0; a<size; a++)//For each spot (region) in the regionsToCheck array
		{
			temp = getIJKVec(n, regionsToCheck[a]);
			if(wellSep(curr, temp, n) == 1)
			{
				//do CoM calculation
			}
			else
			{
				//get array of children, call it childArr
				doForce2(currPlanet, lev-1, BD[], childArr[]);
			}
		}

	}

}
