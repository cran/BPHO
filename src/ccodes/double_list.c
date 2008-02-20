
/****************************************************************************/
double* gen_double_array(int la)
{  if(la > 0) 
     return((double*) R_alloc(la,sizeof(double)));
   else return NULL;
}

/****************************************************************************/
double_vec* gen_double_vec(int la, double* a)
{  double_vec* unit = (double_vec*) R_alloc(1,sizeof(*unit));
   if(la == 0 & a != 0)
      printf("error in 'gen_double_vec': length of non-empty array is 0");
   else{
      if(a == NULL) a = gen_double_array(la);
      unit -> lv = la;
      unit -> v = a;
   }
   return(unit);
}
			 
/****************************************************************************/
void max_double_vec(double_vec* unit, double* max, int *ix_m)
{   int i;
    *max = unit->v[0];
    *ix_m = 0;
    for(i = 1 ; i < unit->lv; i++){     
       if(unit->v[i] > *max){
          *max = unit->v[i];
	  *ix_m = i;
       }
    }
}

/***************************************************************************/
void free_double_vec(double_vec* unit){
    free(unit->v);
    free(unit);
}

/****************************************************************************/
void set_value_double_vec(double_vec* unit, int value){
    int i;
    for(i = 0; i < unit->lv; i++){
       unit->v[i] = value;
    }    
}

/****************************************************************************/
 
void print_double_array(int la, double a[],int item_line){
    int i;
    for(i = 0; i < la; i++){
       printf("%f ",a[i]);
       if( i%item_line == item_line - 1 ) printf("\n");
    }
    
}
