/*******************************************************************************************/
int_list_list* gen_int_list_list(int_list_node* node)
{  int_list_list* list = (int_list_list*) R_alloc(1,sizeof(*list)); 
   if(node == NULL)
      list->ll = 0;
   else list->ll = 1;
   list->next = node;
   list->end = node;
   return(list);
}

/*******************************************************************************************/
/* node is added to list */
void add_int_list_node(int_list_list* list, int_list_node* node)
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

/*******************************************************************************************/
int_list_node* gen_int_list_node(int_list* unit)
{  int_list_node *newnode = (int_list_node*) R_alloc(1,sizeof(*newnode));
   newnode -> value = unit;
   newnode -> next = NULL;
   return newnode;
}
/*******************************************************************************************/
void merge_int_list(int_list* list1, int_list* list2){
   if(list1->ll > 0){
     list1->end->next = list2->next;
     list1->end = list2->end;
   }
   else{
      list1->next = list2->next;
      list1->end = list2->end;
   }
   list1->ll += list2->ll; 
}
/*******************************************************************************************/
void fwrite_int_list_list_into_int_list(int_list_list* list,int lv, 
                                       FILE* fp){
   int_list_node* list_node;
   int_node* node;
   int i;       
   int ll; 
   list_node = list -> next;
    
   if(fwrite(&(list->ll), sizeof(int),1,fp) !=1)
       printf("error in writing the length of list\n");
    
   while(list_node != NULL){
      ll = (list_node -> value-> ll)*lv;
      if(fwrite(&(ll), sizeof(int), 1,fp) != 1)
       printf("error in writing the length of list_node\n");
      node = list_node->value->next;
      while(node != NULL){         
         if(fwrite(node->value->v,sizeof(int),node->value->lv,fp)
            != node->value->lv)
	    printf("error in writting values of integer list\n"); 	         
         node = node->next;
      }  
      list_node = list_node->next; 
   }
}
