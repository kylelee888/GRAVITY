#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

struct list_el {
  int val;
  struct list_el * next;
};
typedef struct list_el item;


void add(item **list, int v)
{
  //   printf("!Adding item %d\n",v);
  item *curr;//"Make a pointer to an item, call it curr"
  curr=(item *)malloc(sizeof(item));//"Allocate an item's worth of memory to the place curr points to"
  curr->val = v;//"Put the value v in the space val pointed to by curr"
  curr->next = *list;//"Put the pointer to the beginning of the list in the space next pointed to by curr"
  *list=curr;//"Set the pointer to the beginning of the list to the pointer curr"
}
void del (item **list, int v)
{
  item *curr;
  curr=*list;//"Start curr at the beginning of the list"
  if (curr == NULL) {return;}
  //  printf("!Deleting item %d... list starts with %d (%d)\n",v,*list,(*list)->val);
  if (curr->val == v) {*list=curr->next; free(curr);}
  else
    while (curr -> next != NULL)//"While the thing in spot next of the item that curr points to is not NULL"
    {
      if (curr->next->val == v)
      {
	free(curr->next); curr->next=curr->next->next;
	break;
      }
      curr=curr->next;
    }
}

void clearList(item **list)
{
  item *curr;//"Make a (pointer to an item) called curr"
  curr=*list;//"Start curr at the beginning of the list"
  if (curr == NULL) {return;}//In case we accidentally call this on an empty list
  while (curr -> next != NULL)
  {
    free(curr->next);//"Free the memory pointed to by next"
    curr->next=curr->next->next;//"Jump curr past the item we just deleted"
  }
}
