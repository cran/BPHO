#include "compress_seq.c"
#include "compress_int.c"

/**************************************************************************/
/* this function compresses the patterns in ``features''. */
void R_compress(int *R_sequence,int *R_order,
                    int *R_n,int *R_p, int features[*R_n][*R_p],
                    int *R_nos_fth,int *R_no_cases_ign,
                    char** files, int* R_quiet, int* R_do_comp)
{  if(R_sequence[0] ==1)
      compress_seq(*R_n,*R_p,features,R_nos_fth,
                   *R_no_cases_ign,files[0],*R_quiet,*R_do_comp);
   else
      compress_cls(*R_order,*R_n,*R_p,features,R_nos_fth,
                   *R_no_cases_ign,files[0],*R_quiet,*R_do_comp);

}

/*************************************************************************/
int_vec* take_cases(int_vec* m_cases,int n,int p,int features[n][p],
                    int id_fth,int value)
{  int i,nos_sub=0;
   int* ind_pickup = (int*) R_alloc(m_cases->lv,sizeof(int));

   for(i=0; i< m_cases->lv; i++)
   { if(features[m_cases->v[i]][id_fth] == value)
     {  ind_pickup[i] =1;
        nos_sub ++;
     }
     else
        ind_pickup[i] =0;
   }
   int_vec* sub_cases = gen_int_vec(nos_sub,NULL);
   if(nos_sub > 0){
      int id_sub=0;
      for(i=0; i < m_cases -> lv; i++)
      {  if(ind_pickup[i]==1)
           sub_cases->v[id_sub++] = m_cases->v[i];
      }
   }

   return sub_cases;
}


/***************************************************************************/
int is_unique(int n,int p,int features[n][p],
             int_vec* sub_cases,int id_fth)
{
   int i;
   int a = features[sub_cases->v[0]][id_fth];
   for(i=1; i < sub_cases->lv; i++){
      if(features[sub_cases->v[i]][id_fth] != a)
         return(0);
   }
   return(1);

}

/****************************************************************************/
/*this function loads the information about compressed parameters from 'file'
to objects defined outside */
void load_compress
      (int_list* cases_g, int_list* ptns_g,
       int* P_sequence,int* P_order,int* P_n, int* P_p,
       int** P_sigmas_g,char* file, int quiet)
{
    FILE *fp;
    fp = fopen_chk(file,"rb");
    fread_chk(P_sequence,sizeof(int),1,fp,file);
    fread_chk(P_order,sizeof(int),1,fp,file);
    fread_chk(P_n,sizeof(int),1,fp,file);
    fread_chk(P_p,sizeof(int),1,fp,file);
    fread_int_list(cases_g,fp);
    fread_int_list(ptns_g,fp);
    *P_sigmas_g = gen_int_array(ptns_g->ll*(*P_order+1));
    fread_chk(*P_sigmas_g,sizeof(int),ptns_g->ll*(*P_order+1),fp,file);

    fclose(fp);

}

/****************************************************************************/
/*this function is to be called by R for displaying the
  pattern information for specified indice of groups*/
void R_display_compress
     (char** files, int* no_show, int* gids, int info[6])
{
    /*the blocks of loading compression information from `file'*/
    int n,p,order,sequence,no_g,*sigmas_g,no_ptns;
    int_list *cases_g=gen_int_list(NULL);
    int_list *ptns_g=gen_int_list(NULL);

    load_compress(cases_g,ptns_g,&sequence,&order,&n,&p,
                  &sigmas_g,files[0],0);

    no_g = cases_g ->ll;
    no_ptns = sum_int(no_g*(order+1),sigmas_g);

    int id_g = 0;
    int id_list = 0;
    int_node* ptn_node = ptns_g -> next;
    int_node* cse_node = cases_g -> next;
    while(ptn_node!=NULL & id_g < *no_show){
        if(id_list==gids[id_g]){

        printf(
"\n*********************** Information of Group %d ****************************\n",
        gids[id_g] + 1);//adding 1 for R user interface 

       if(sequence==1){
        printf("Superpatterns:\n");
        printf("(Last number indicates the smallest order in this group)\n\n");
       }
       else{
        printf("Superpatterns:(Every two lines is a superpattern)\n");
        printf("Last number in 1st line indicates the smallest order in this group.\n");
        printf("Last number in 2nd line indicates number of fixed positions.\n\n");
       }
        print_int_vec(ptn_node->value,p+1);

        printf("Expression:(Indice of cases expressing this superpattern)\n\n");
        print_expression(cse_node->value, 15); printf("\n");

        printf("number of Sigmas for each order from 0 to %d:\n\n",order);
        print_int_array(order+1,sigmas_g+gids[id_g]*(order+1),order+1);
            id_g ++;
        }
        id_list ++;
        ptn_node = ptn_node -> next;
        cse_node = cse_node -> next;
   }

   info[0] = sequence;info[1] = order;
   info[2] = no_g; info[3] = no_ptns;
   info[4] = n; info[5] = p;


}
