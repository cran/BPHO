typedef struct double_vector { int lv;
                            double *v;
	                  } double_vec;


typedef struct double_vec_node { double_vec* value;
	                      struct double_vec_node *next; 
	                    } double_node;
		  
typedef struct double_vec_list { int ll; 
                              double_node *next;
		              double_node *end;
	                    } double_list;
			 
double* gen_double_array(int la);

double_vec* gen_double_vec(int la, double* a);

void max_double_vec(double_vec* unit, double* max, int *ix_m);

void free_double_vec(double_vec* unit);

void set_value_double_vec(double_vec* unit, int value);
