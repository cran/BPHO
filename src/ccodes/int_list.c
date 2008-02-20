/****************************************************************************/
/* node is added to list */
void add_int_node(int_list* list, int_node* node)
{  if(list -> ll > 0){
     list->end->next = node;
     list->end = node;      
   }
   else{
       list->next = node;
       list->end = node;
   } 
   list->ll ++;
}

/****************************************************************************/
int_list* gen_int_list(int_node* node)
{  int_list* list = (int_list*) R_alloc(1,sizeof(*list)); 
   if(node == NULL)
      list->ll = 0;
   else list->ll = 1;
   list->next = node;
   list->end = node;
   return(list);
}

/****************************************************************************/
int_node* gen_int_node(int_vec* unit)
{  int_node *newnode = (int_node*) R_alloc(1,sizeof(*newnode));
   newnode -> value = unit;
   newnode -> next = NULL;
   return newnode;
}

/****************************************************************************/
int_vec* gen_int_vec(int la, int* a)
{  int_vec* unit = (int_vec*) R_alloc(1,sizeof(*unit));
   if(la == 0 & a != 0)
      printf("error in 'gen_int_vec': length of non-empty array is 0");
   else{
      if(a == 0)  a = gen_int_array(la);
      unit -> lv = la;
      unit -> v = a;
   }
   return(unit);
}
/****************************************************************************/
int* gen_int_array(int la)
{  if(la > 0) 
     return((int*) R_alloc(la,sizeof(int)));
   else return NULL;
}

/****************************************************************************/
void add_int_element(int_vec* unit, int a)
{  int* newv = (int *) R_alloc(unit->lv+1,sizeof(int));
   int i=0;
   
   for(i=0;i< unit->lv;i++)
      newv[i] = unit->v[i];
      
   newv[unit->lv] = a;
   
   unit->v = newv;
   unit->lv++;
}

/**************************************************************************/
void merge_int_vec(int_vec* unit1, int_vec* unit2)
{   int lv = unit1->lv + unit2->lv;
    int *v = (int*) R_alloc(lv, sizeof(int));
    int i;
    for(i = 0; i < unit1->lv; i++)
      v[i] = unit1->v[i];
    for(i = unit1->lv; i < lv; i++)
      v[i] = unit2->v[i - unit1->lv];
    
    unit1 -> v = v;
    unit1 -> lv = lv;    
}

/***************************************************************************/
int is_same_vec(int_vec* unit1, int_vec* unit2)
{  if(unit1 -> lv != unit2 -> lv)
      return 0;
    else
    {  int i;
       for(i=0; i < unit1->lv; i++)
         if(unit1->v[i] != unit2->v[i])
	   return 0;
       return 1;
    }    
}

/***************************************************************************/
void print_int_vec(int_vec* unit,int items_line)
{  int i;
   for(i = 0; i < unit -> lv; i++){
      printf("%d ",unit -> v[i]);
      if(i%items_line == items_line - 1) printf("\n");
   }
   printf("\n");
}

/**************************************************************************/
void print_int_list(int_list * list, int items_line)
{ /*int i,k;*/
  int_node* p = list->next;
  
  printf("Totally %d vectors in list\n", list->ll);
  int i;
  for(i = 0; i < list->ll; i++){
     printf("[[%d]]:\n", i);
     int k;
     print_int_vec(p->value,items_line);          
     p = p->next;     
  }  
}

/**************************************************************************/
void fwrite_int_list(int_list *list, FILE* fp)
{   int_node* node;
    int i;       
    node = list -> next;
    
    if(fwrite(&(list->ll), sizeof(int), 1,fp) !=1)
       printf("error in writing the length of list\n");
    
    while(node != NULL){
       if(fwrite(&(node->value->lv), sizeof(int), 1,fp) !=1)
          printf("error in writing the length of a vector\n");
       
       if(fwrite(node->value->v,sizeof(int),node->value->lv,fp)
          != node->value->lv)
	  printf("error in writting values of integer list\n"); 	         
       node = node->next;
    }

}

/***************************************************************************/
void fread_int_list(int_list *list, FILE* fp)
{   int_node* node;    
    int* v;
    int lv;
    int i;
    int ll;
    if(fread(&ll,sizeof(int),1,fp)!=1)
       printf("error in reading the length of a list\n");
    for(i=0; i < ll;i++){
       if(fread(&lv,sizeof(int),1,fp) != 1)
         printf("error in reading the length of a vector\n"); 
       v = gen_int_array(lv);
       if(fread(v,sizeof(int),lv,fp) != lv)
         printf("error in reading value of a vector\n");        
       node = gen_int_node(gen_int_vec(lv,v));       
       add_int_node(list,node);         
    }    
}
/***************************************************************************/
void copy_array(int la, int *a, int *b)
{  int i;
   for(i = 0; i < la; i++){
      a[i] = b[i];
   }
   
}

/***************************************************************************/
int is_all_same_array(int la, int* a)
{  int i;
   for(i=1; i < la; i++){
     if(a[i] != a[0]) 
        return(0);
   }
   return(1);
}

/***************************************************************************/
int_vec* seq(int start,int end)
{  int_vec* q = gen_int_vec(end-start+1,NULL);
   int i;
   for(i = 0; i < q -> lv; i++)
   {  q -> v[i] = i + start;
   }
   return q;
}
/***************************************************************************/
void free_int_list(int_list* list){
   int_node *node, *node_next;
   node = list-> next;
   while(node!=0){
       node_next = node -> next;
       free(node->value->v);
       free(node->value);
       free(node);
       node = node_next;
   }

}
/***************************************************************************/
void set_value_int_vec(int_vec* unit, int value){
    int i;
    for(i = 0; i < unit->lv; i++){
       unit->v[i] = value;
    }    
}
/***************************************************************************/
int sum_int_vec(int_vec* unit){
   int i, s;
   for(i = 0; i < unit->lv; i++){
      s += unit->v[i];
   }
   return(s);
}
/**************************************************************************/
int_vec* sub_int_vec(int_vec* unit, int_vec* indice_sub){
    int_vec* sub_unit = gen_int_vec(unit->lv - indice_sub->lv,NULL);
    int out[unit->lv];
    int i;
    for(i = 0; i< unit->lv; i++) out[i] = 0;
    for(i = 0; i < indice_sub->lv; i++) out[indice_sub->v[i]] = 1;
    int j = 0;
    
    for(i = 0; i < unit->lv; i++) 
       if(!out[i]){
         sub_unit->v[j] = unit->v[i]; 
	 j++;
       }
    return(sub_unit);
}

/************************************************************************/
void free_int_vec(int_vec* unit){
    free(unit->v);
    free(unit);
}

void set_value_array(int la, int a[], int value){
    int i;
    for(i = 0; i < la; i++){
       a[i] = value;
    }
}




