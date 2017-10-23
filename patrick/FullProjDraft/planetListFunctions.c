typedef struct
{
	int plnNum;
	struct pln* next;
} pln;

pln *plnList[8+6+512+4096];//Does not include an index for the planets in level 0

void add(int num, pln **list)
{
	pln *new;//"Make a pointer to a pln, call it 'new'"
	new=(pln *)malloc(sizeof(el));//"Allocate a pln's worth of memory to the space pointed to by 'new'"
	new->plnNum = num;//Follow 'new' to the pln it points to, set plnNum of that pln to the value num
	new->next = *list;//Follow 'new' to the pln it points to, set next of that pln to pointer to the beginning of the old list
	*list = new;//Set the pointer to the beginning of the list to the pointer 'new'
}

void clearList(pln **list)
{
	pln *curr = *list;
	pln *next;
	while(curr!=NULL)
	{
		next = curr->next;
		free(curr);
		curr = next;
	}
}

void popLists(planet BD[])
{
	





}
