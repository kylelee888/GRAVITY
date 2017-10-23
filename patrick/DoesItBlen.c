typedef struct
{
	int i;
	int j;
	int k;
} intVec;

int getReg(int n, intVec ijk)
{
	return(intVec.i + n*intVec.j + n*n*intVec.k);

}

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

void popAdjReg(int adjReg[], int reg, int n)
{
	int z = 0;
	int a, b, c, i, j, k;
	intVec ayy=getIJKVec(n, reg);//Put ijk coordinates of interesting box into an intVec
	intVec lmao;//For later
	for(a=-1; a<2; a++)//Start one space behind i, end one space ahead
	{
		i=ayy.i+a;
		if(i<0 || i>n-1)//Are we outside of the our space (L0)?
		{continue;}//Yes. a++ and try again
		else//No. Keep on keepin on.
		{
			for(b=-1; b<2; b++)//Start one space behind j, end one space ahead
			{
				j=ayy.j+b;
				if(j<0 || j>n-1)//Outside space?
				{continue;}//b++ and try again
				else
				{
					for(c=-1; c<2; c++)//Start one space behind k, end one space ahead
					{
						k=ayy.k+c;
						if(k<0 || k>n-1)
						{/*Insert 4 coins to */continue;}
						else
						{//--------Exclude own box?--------
							lmao.i=i;
							lmao.j=j;
							lmao.k=k;
							adjReg[z]=getReg(n, lmao);//Find which region is given by lmao's ijk coordinates, put it into the array of adjacent regions adjReg
							z++;
						}
					}
				}
			}
		}
	}
}
