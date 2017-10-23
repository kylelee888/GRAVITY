typdef struct
{
    int i;
    int j;
    int k;
} intVec;

intVec getIJKVec(int n, int reg)
{
    intVec thot;
    int nSq = n*n;
    thot.k = reg/nSq;//Make nSq instead of doing reg/n/n to avoid bad rounding
    int blah = reg - nSq*thot.k;//Give me what's left of reg after subtracting k
    thot.j = blah/n;
    thot.i = reg - thot.k - thot.j;
    return(thot);
}

int wellSep(intVec curr, intVec temp, int n)
{
	/*These six variables represent the possible indeces of all regions adjacent to curr, the region containing our planet of interest.
	All 27-1=26 of the adjacent regions will have indeces that are some combination of these six numbers.*/
	int xl=curr.i-1, xr=curr.i+1;
	int yl=curr.i-1, yr=curr.j+1;
	int zl=curr.k-1, zr=curr.k+1;//l -> "left", r -> "right"
	
	if(xl<0){xl=0;}//"Is xl outside of ('left of') our space? If so, set it to zero"
	if(xr>n-1){xr=n-1;}//"Is xr outside of ('right of') our space? If so, set it to n-1, the maximum index number for the current level
	if(yl<0){yl=0;}
	if(yr>n-1){yr=n-1;}
	if(zl<0){zl=0;}
	if(zr>n-1){zr=n-1;}

	//Use the separation of curr and temp in each direction to determine if well separated
	if(abs(temp.i-xl) > 2 || /*<- This is the 'or' operator*/ abs(temp.i-xr) > 2)//If the difference is zero, curr and temp share the same i index. If the difference is 1, curr and temp may be adjacent depending on the other indeces. If the difference is 2 or more, curr and temp can not possibly be adjacent, so 2 is the minimum difference for being well separated
	{
		return 1;//Well separated
	}
	else if(abs(temp.j-yl) > 2 || abs(temp.j-yr) > 2)
	{
		return 1;
	}
	else if(abs(temp.k-zl) > 2 || abs(temp.k-zr) > 2)
	{
		return 1;
	}
	//If it makes it this far, curr and temp are adjacent and therefore not well separated, so:
	else
	{
		return 0;
	}

}
