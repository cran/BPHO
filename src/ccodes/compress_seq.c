
/*************************************************************************/
void compress_seq(int n, int p, int features[n][p],int* nos_fth,
                  int no_cases_ign, char* file, int quiet, int do_comp)
{
    int_list* cases_g = gen_int_list(NULL);
    int_list* ptns_g = gen_int_list(NULL);   
    int_vec* m_cases = seq(0,n-1);
    int* first_ptn = gen_int_array(p+1);
    FILE *fp;

    set_value_array(p+1,first_ptn,0);    
    if(!quiet) printf("Diverging ...\n");
    diverge_seq(n,p,features,nos_fth,m_cases,first_ptn,no_cases_ign,
	    cases_g,ptns_g,do_comp);       
    if(!quiet)
       printf("Finished diverging,found %d patterns\n",cases_g->ll); 
              
    int sigmas_g[ptns_g->ll][p+1] ;
         
    find_sigmas_seq(ptns_g,p,sigmas_g);    
    if(!quiet){ 
       printf("Number of expressed patterns: %d\n", 
         sum_int((ptns_g->ll)*(p+1), &sigmas_g[0][0]));
       printf("Writing the compression information to file %s ... \n",
         file);
    }
    fp = fopen_chk(file,"wb"); 
    int sequence=1;
    fwrite_chk(&sequence,sizeof(int),1,fp,file);
    /*write p as depth*/
    fwrite_chk(&p,sizeof(int),1,fp,file);     
    fwrite_chk(&n,sizeof(int),1,fp,file);        
    fwrite_chk(&p,sizeof(int),1,fp,file);           
    fwrite_int_list(cases_g,fp);
    fwrite_int_list(ptns_g,fp);
    fwrite_chk(sigmas_g,sizeof(int),(ptns_g->ll)*(p+1), fp,file); 
    fclose(fp);        
}

/**************************************************************************/
void diverge_seq(int n,int p,int features[n][p], int nos_fth[p], 
            int_vec* m_cases,int m_ptn[p+1],int no_cases_ign,
	    int_list* cases_g, int_list* ptns_g, int do_comp)
{   int i_fth,class,*new_ptn;
    int_vec* sub_cases;
    i_fth = m_ptn[p];
    if(do_comp == 1){
       while(is_unique(n,p,features,m_cases,i_fth) & i_fth < p)
       {  m_ptn[i_fth] = features[m_cases->v[0]][i_fth];
          i_fth ++;
       }
    }
    /*else i_fth++;*/
    
    add_int_node(cases_g, gen_int_node(m_cases)); 
    add_int_node(ptns_g, gen_int_node(gen_int_vec(p+1,m_ptn)));
    if(i_fth < p){
      for(class = 1; class < nos_fth[i_fth]+1; class++){
         sub_cases = take_cases(m_cases,n,p,features,i_fth,class);
	 if(sub_cases->lv > no_cases_ign){
	    new_ptn = gen_int_array(p+1);
	    copy_array(p+1,new_ptn,m_ptn);
	    new_ptn[i_fth] = features[sub_cases->v[0]][i_fth];
	    new_ptn[p] = i_fth + 1;
	    diverge_seq(n,p,features,nos_fth,sub_cases,new_ptn,no_cases_ign,
	            cases_g,ptns_g,do_comp);
	 }
      } 
    }  
}

/**************************************************************************/
void find_sigmas_seq(int_list* ptns_g,int p,int sigmas_g[][p+1])
{  
    int_node* node_ptn;
    int* ptn;
    node_ptn = ptns_g -> next;
    int i_g, order;
    set_value_array((ptns_g->ll)*(p+1), &sigmas_g[0][0], 0);
    i_g = 0;
    while(node_ptn != 0){
       ptn = node_ptn -> value -> v;
       
       for(order = ptn[p]; order < p+1; order++){
           sigmas_g[i_g][order]++;
           if(ptn[order]==0) break;
       }
       i_g++;
       node_ptn = node_ptn -> next;         
    }    
}
  




