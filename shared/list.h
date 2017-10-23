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

