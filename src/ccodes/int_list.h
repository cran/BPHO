/*This file contains the definition of int vector and int list as well as some
basice operations on them*/

typedef struct int_vector { int lv;
                            int *v;
	                  } int_vec;


typedef struct int_vec_node { int_vec* value;
	                      struct int_vec_node *next; 
	                    } int_node;
		  
typedef struct int_vec_list { int ll; 
                           int_node *next;
		           int_node *end;
	                 } int_list;
			 
void add_int_node(int_list* list, int_node* node);

int_list* gen_int_list(int_node* node);

int_node* gen_int_node(int_vec* unit);

int_vec* gen_int_vec(int la, int* a);

int* gen_int_array(int la);

void merge_int_lists(int_list* list1, int_list* list2);

void add_int_element(int_vec* unit, int a);

void merge_int_vec(int_vec* unit1, int_vec* unit2);

int is_same_vec(int_vec* unit1, int_vec* unit2);

void print_int_vec(int_vec* unit,int items_line);

void print_int_list(int_list * list, int items_line);

void fwrite_int_list(int_list *list, FILE* fp);

void fread_int_list(int_list *list, FILE* fp);

void copy_array(int la, int *a, int *b);

int is_all_same_array(int la, int* a);

int_vec* seq(int start,int end);

void free_int_list(int_list* list);

void set_value_int_vec(int_vec* unit, int value);

int sum_int_vec(int_vec* unit);

int_vec* sub_int_vec(int_vec* unit, int_vec* sub);

void free_int_vec(int_vec* unit);
void print_expression(int_vec* unit,int items_line);
