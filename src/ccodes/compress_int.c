
/****************************************************************************/
void compress_cls(int order,int n, int p, int features[n][p],
                  int* nos_fth,int no_cases_ign, char* file, int quiet,
              int do_comp){
    FILE *fp;
    int_list* cases_g = gen_int_list(NULL);
    int_list_list* ptns_g = gen_int_list_list(NULL);
   
    int_vec* m_cases = seq(0,n-1);
    int_vec* m_fths = seq(0,p-1);
    int i;
    int max_nos_fth = 0;
    
    for(i = 1; i < p; i++){
       if(nos_fth[i] > max_nos_fth) max_nos_fth = nos_fth[i];
    } 
    int_vec* NOs_class = gen_int_vec(max_nos_fth,NULL);
    
    int* first_ptn = gen_int_array(2*p+2);
    for(i=0; i < p+1; i++) first_ptn[i] = 0;
    for(i=p+1; i < 2*(p+1) - 1; i++) first_ptn[i] = 1;
    first_ptn[i] = p;
    if(!quiet) printf("Diverging ...\n");

    diverge_cls(n,p,features,m_cases,m_fths,nos_fth,cases_g,ptns_g,
            order,first_ptn,NOs_class,no_cases_ign,do_comp);
            
    if(!quiet){
              printf("Finished diverging,found %d patterns\n",cases_g->ll); 
              printf("Merging patterns with the same expression ... \n"); 
    }       
    /*removing the patterns with same expression*/
    if(do_comp == 1) {
       trim_ptns(cases_g,ptns_g);
       if(!quiet) 
          printf("Finished merging, found %d patterns\n",cases_g->ll);
    } 
    
    int sigmas_g[ptns_g->ll][order+1];
    set_value_array((ptns_g->ll)*(order+1), &sigmas_g[0][0], 0);     
    find_sigmas_interact(ptns_g,order,p,sigmas_g);
    if(!quiet) 
       printf("Number of expressed patterns: %d \n", 
         sum_int((ptns_g->ll)*(order+1), &sigmas_g[0][0]));
       
    /*write the data relavant to the patterns to file_out[0]*/
    if(!quiet) printf("Writing the cases_ptns to file %s ... \n",file);
           
    fp = fopen_chk(file,"wb"); 
    int sequence=0;
    fwrite_chk(&sequence,sizeof(int),1,fp,file);    
    fwrite_chk(&order,sizeof(int),1,fp,file);     
    fwrite_chk(&n,sizeof(int),1,fp,file);        
    fwrite_chk(&p,sizeof(int),1,fp,file);
    fwrite_int_list(cases_g,fp);
    fwrite_int_list_list_into_int_list(ptns_g,2*p+2,fp);
    fwrite_chk(sigmas_g,sizeof(int),(ptns_g->ll)*(order+1), fp,file);      
    
    fclose(fp);        
}

void diverge_cls(int n, int p,int features[n][p], 
            int_vec* m_cases,int_vec* m_fths,int* nos_fth, 
            int_list* cases_g, int_list_list* ptns_g, int order, 
            int* m_ptn, int_vec* NOs_class, int no_cases_ign,
            int do_comp){
   
  int i,j,id_fth,f;
  int* new_ptn;
  int_vec* sub_cases;
  int_vec* sub_fths;  
  double_vec* entropy;double max_entropy;
  int_vec *ix_m = gen_int_vec(1,NULL);
  
  /*calculating entropies for all features in m_fths */ 
  entropy = entr_sub_sub(n,p,features,m_cases,m_fths,NOs_class);    
  /*find the biggest entropy and the associated indice*/    
  max_double_vec(entropy,&max_entropy,ix_m->v);
  /*grouping patterns*/
  if(max_entropy < 1e-10 & do_comp == 1)
  {      
      new_ptn = gen_int_array(2*p+2);copy_array(2*p+2,new_ptn,m_ptn);
      for(f = 0; f < m_fths->lv; f++){
         id_fth = m_fths->v[f];
         new_ptn[id_fth] = features[m_cases->v[0]]                                   [id_fth];
         new_ptn[p+1+id_fth] = 0;                      
      }
      new_ptn[p+1+p] -= m_fths->lv; 
      add_int_list_node(ptns_g,gen_int_list_node(gen_int_list(
                        gen_int_node(gen_int_vec(2*p+2,new_ptn)))));
      add_int_node(cases_g,gen_int_node(m_cases));
  }
  /*if no features with entropy 0*/                   
  else 
  {                      
     id_fth = m_fths->v[ix_m->v[0]];      
     sub_fths =  sub_int_vec(m_fths,ix_m);  
     
     //pass on the mother pattern with a id_fth removed from unconsidered
     if(sub_fths->lv == 0)
     {    
       add_int_list_node(ptns_g,gen_int_list_node(gen_int_list(
                         gen_int_node(gen_int_vec(2*p+2,m_ptn)))));
       add_int_node(cases_g,gen_int_node(m_cases));
     }
     else 
       diverge_cls(n,p,features,m_cases,sub_fths,nos_fth,
                   cases_g,ptns_g,order,m_ptn,NOs_class,
                   no_cases_ign,do_comp);
     
     //fixing id_fth in pattern to a particular value, and diverge
     if(m_ptn[p] < order)
     {
       for(i=1; i <= nos_fth[id_fth]; i++)
       {
         sub_cases = take_cases(m_cases,n,p,features,id_fth,i);
         if(sub_cases->lv > no_cases_ign)
         {                           
           new_ptn = gen_int_array(2*p+2);copy_array(2*p+2,new_ptn,m_ptn); 
           new_ptn[id_fth] = i; new_ptn[p]++;

           if(sub_fths->lv == 0 || new_ptn[p] == order)
           {
             add_int_node(cases_g,gen_int_node(sub_cases));              
             add_int_list_node(ptns_g,gen_int_list_node(
                               gen_int_list(gen_int_node(
                               gen_int_vec(2*p+2,new_ptn)))));
           }
           else
           {
            diverge_cls(n,p,features,sub_cases,
                          sub_fths,nos_fth,cases_g,ptns_g,order,
                          new_ptn,NOs_class,no_cases_ign,do_comp);
           }                 
         }
       }
     }
  }
}


/******************************************************************************/
double calc_entropy(int_vec* NOs_class, int N){
   int i;
   double raw_ent = 0;   
   for(i = 0; i < NOs_class->lv; i++){
      if(NOs_class->v[i])
         raw_ent += (NOs_class->v[i]) * log(NOs_class->v[i]);      
   }
   if( N == 0 ) N = sum_int_vec(NOs_class);
   
   return( (- raw_ent/N + log(N))); }

/******************************************************************************/
/*this function calculates the NEGATIVE entropy of features `sub_fths'
restricted on `sub_cases', return the entropy values in form of ``double_vec*''
*/

double_vec* entr_sub_sub(int n, int p, int features[n][p], 
                      int_vec* sub_cases,int_vec* sub_fths, 
                      int_vec* NOs_class){
   int f, k;
   
   double_vec* entropy = gen_double_vec(sub_fths->lv,NULL);
   if(sub_cases->lv > 1){
     for(f = 0; f < sub_fths->lv; f++){
        set_value_int_vec(NOs_class,0);
        for(k = 0; k < sub_cases->lv; k++){
            NOs_class->v
            [ features[ sub_cases->v[k] ][ sub_fths->v[f]] - 1 ]++;
        }
        entropy->v[f] = calc_entropy(NOs_class, sub_cases->lv); 
     }
   }
   else set_value_double_vec(entropy,0);
      
   return(entropy);   
}

/******************************************************************************/
void trim_ptns(int_list *cases_g, int_list_list *ptns_g)
{  int i;
   int_node *ps_cases,*pm_cases,*pm_cases_next;
   int_list_node *ps_ptns,*pm_ptns,*pm_ptns_next;
   ps_cases = cases_g -> next;
   ps_ptns = ptns_g -> next;
      
   while(ps_cases ->next != NULL) { 
       pm_cases = ps_cases;
       pm_ptns = ps_ptns;
       while(pm_cases -> next != NULL){
          pm_cases_next = pm_cases -> next;
          pm_ptns_next = pm_ptns -> next;
          if(is_same_vec(pm_cases_next->value, ps_cases->value)){
             pm_cases -> next = pm_cases_next -> next;
             pm_ptns -> next = pm_ptns_next -> next;
             cases_g -> ll --;
             ptns_g -> ll --;
             merge_int_list(ps_ptns->value,pm_ptns_next->value);
                   
          }
          else{
              pm_cases = pm_cases -> next; 
              pm_ptns = pm_ptns -> next;
          }
       }
       cases_g -> end = pm_cases;
       ptns_g -> end = pm_ptns;
       if(ps_cases -> next != NULL){ 
          ps_cases = ps_cases -> next;
          ps_ptns = ps_ptns -> next;
       }                       
   }    
}
            
/****************************************************************************/
void add_sigmas_patterns(int id_g,int no_g, int order, int p, 
                        int sigmas_g[][order+1], 
                        int pattern[][p+1]) {
    int no_nofixed = p - pattern[1][p];
    int no_nozero = pattern[0][p];    
    int no_nozero_lack = imin2(order - no_nozero, no_nofixed);
    int i;

    for(i = 0; i <= no_nozero_lack; i++){
      sigmas_g[id_g][i+no_nozero] += choose(no_nofixed, i);
    }

}

/****************************************************************************/
void find_sigmas_interact(int_list_list* ptns_g,int order,int p,
                  int sigmas_g[][order+1]){
    int id_g = 0;
    int_list_node* list_node = ptns_g->next;
    int_node* node;
    while(list_node != NULL){
        node = list_node -> value -> next;
        while(node != NULL){

             add_sigmas_patterns(id_g,ptns_g->ll,order,p,
                                 sigmas_g,node->value->v);
             node = node -> next;
        }
        list_node = list_node -> next;
        id_g ++;
    }   
}


