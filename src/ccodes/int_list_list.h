typedef struct int_vec_list_node{int_list* value;
                              struct int_vec_list_node *next;
			      } int_list_node;
			      			 
typedef struct int_vec_list_list{int ll;
                              int_list_node* next;
			      int_list_node* end;
			      } int_list_list;
			 
			
int_list_list* gen_int_list_list(int_list_node* node);

void add_int_list_node(int_list_list* list, int_list_node* node);

int_list_node* gen_int_list_node(int_list* unit);
void fwrite_int_list_list_into_int_list(int_list_list* list,int lv, FILE* fp);
