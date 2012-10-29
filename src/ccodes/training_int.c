
double dstable(double x, double mu, double width, int alpha)
{   if(alpha == 1)
        return(dcauchy(x,mu,width,1));
    if(alpha == 2)
        return(dnorm(x,mu,sqrt(width),1));
}

/**************************************************************************/
/*y_tr is assumed to coded by 0,1,...,no_cls-1 */
void R_training
     (char** mc_file,char** cases_ptns_file, int y_tr[],
      int no_cls[], int iters_mc[],int iters_ss[], int iters_hp[],
      double w_ss[], double w_hp[], int m_ss[], int m_hp[],
      int alpha[],double a_sigmas[], double b_sigmas[],
      double log_sigmas[]){
  training(
       mc_file[0],cases_ptns_file[0],y_tr,
       no_cls[0],iters_mc[0],iters_ss[0],iters_hp[0],
       w_ss[0],w_hp[0],m_ss[0],m_hp[0],alpha[0],
       a_sigmas,b_sigmas,log_sigmas);

}
/***************************************************************************/
/*y_tr is assumed to coded by 0,1,...,no_cls-1 */
void training
     (char mc_file[],char cases_ptns_file[], int y_tr[],
      int no_cls, int iters_mc,int iters_ss, int iters_hp,
      double w_ss, double w_hp, int m_ss, int m_hp,
      int alpha, double a_sigmas[], double b_sigmas[],
      double log_sigmas[]) { 
    int n,p,no_g,sequence,order,no_hp,i_mc, i_ss, i_hp, i,i_case,
        i_g, i_cls, i_sgm,J,K,next_iter_mc=0, cls_begin=0;
    int_list* cases_g = gen_int_list(NULL);
    int_list* ptns_g = gen_int_list(NULL);
    int* sigmas_g;
    int_vec* cases;
    FILE *fp_mc;

    double beta_o,log_sigma_o,U,z,beta_n,log_sigma_n,beta_L,beta_R,
           sigma_L,sigma_R,eval_ss,eval_hp;
    /*loading compression information*/
    load_compress(cases_g,ptns_g,&sequence,&order,&n,&p,&sigmas_g,
                  cases_ptns_file,1);
    no_g = cases_g -> ll;
    /*requesting space for markov chaing sampling*/
    double
         betas[no_g][no_cls],widths_g[no_g],widths_g_p[no_g],lv[n][no_cls],
         log_sum_exp_lv[n],sum_log_likev[1],sum_log_prior_betas[1],
         log_prior_sigmas[order+1],sum_log_prior_sigmas[1],
	 log_postv[1],sum_log_prior_sigmas_p,sum_log_prior_betas_p,
         sum_log_likev_p, lv_p[n][no_cls];    
    long size_each_iter=sizeof(int)+sizeof(double)*(no_g*(no_cls+1)+order+7);
    
    if(no_cls == 2) cls_begin = 1;
    no_hp = 0;
    for(i_sgm = 0; i_sgm < order + 1; i_sgm ++) {
        if(a_sigmas[i_sgm] > 1e-10) no_hp ++;
    }
    
    //creating new log file or read from the existing log file
    fp_mc = fopen(mc_file,"rb");
    if(fp_mc == NULL){
        fp_mc = fopen_chk(mc_file,"wb");
        /*creating log file */
        fwrite_chk(&no_cls,sizeof(int),1,fp_mc,mc_file);
        fwrite_chk(&no_g,sizeof(int),1,fp_mc,mc_file);
        fwrite_chk(&order,sizeof(int),1,fp_mc,mc_file);
	fwrite_chk(&alpha,sizeof(int),1,fp_mc,mc_file);
	fwrite_chk(&(next_iter_mc),sizeof(int),1,fp_mc,mc_file);
        set_double_array(no_g*no_cls,betas,0);
        fclose(fp_mc);
    }
    else{
        //if log file exists, read the log_sigmas and betas in
        fseek(fp_mc,-size_each_iter,SEEK_END);
        fread_chk(betas,sizeof(double),no_g*no_cls,fp_mc,mc_file);
        fread_chk(log_sigmas,sizeof(double),order+1,fp_mc,mc_file);
        fseek(fp_mc,(6+no_g)*sizeof(double),SEEK_CUR);
        fread_chk(&next_iter_mc,sizeof(int),1,fp_mc,mc_file);
        fclose(fp_mc);
    }

    /*initialize the markov chain state*/
    initialize_mc_state(alpha,cls_begin,no_cls,order,n,no_g,
       widths_g,log_sigmas,betas,y_tr,cases_g,a_sigmas,b_sigmas,sigmas_g,
       sum_log_prior_betas,log_prior_sigmas,sum_log_prior_sigmas,
       log_postv,lv,sum_log_likev,log_sum_exp_lv);
       
    /*starting mc updating with slice sampling*/
    GetRNGstate();
    U = runif(0,1);
    PutRNGstate();
    
    iters_mc += next_iter_mc;
    
    while(next_iter_mc < iters_mc) {
       next_iter_mc++;
       eval_ss = 0;
       eval_hp = 0;
       //updating betas
       i_ss = 0;
       for(i_ss = 0; i_ss < iters_ss; i_ss++){
            //Rprintf("I am starting updating betas\n");
           for(i_g = 0; i_g < no_g; i_g++) {
                cases = find_cases_gid(cases_g, i_g);
                for(i_cls = cls_begin; i_cls < no_cls; i_cls++){
                   GetRNGstate();
                   U = runif(0,1);
                   PutRNGstate();
                   if(U == 1) break;
                   z =  log_postv[0] + log(U);

                   /*  placing the first interval */
                   GetRNGstate();
                   U = runif(0,1);
                   PutRNGstate();
                   beta_L = betas[i_g][i_cls] - (w_ss) * U;
                   beta_R = beta_L + (w_ss);
                   
                   GetRNGstate();
                   U = runif(0,1);
                   PutRNGstate();
                   J = floor((m_ss) * U);
                   K = (m_ss-1) - J;

                   beta_o = betas[i_g][i_cls];
                   //remove terms in sum_log_likev involving this beta
                   sum_log_likev_p = sum_log_likev[0];
                   for(i = 0; i < cases->lv; i++)
                   {   i_case = cases->v[i];
                       if(i_cls == y_tr[i_case])
                          sum_log_likev_p -= lv[i_case][i_cls];
                       sum_log_likev_p += log_sum_exp_lv[i_case];
                       lv_p[i_case][i_cls] = lv[i_case][i_cls] - 
                                             betas[i_g][i_cls];
                   }
                   //remove terms in sum_log_prior_betas involving this beta
                   sum_log_prior_betas_p = sum_log_prior_betas[0] -
                       dstable(betas[i_g][i_cls],0,widths_g[i_g],alpha);

                   
                   /*stepping left*/
                   while((J--) > 0){
                      betas[i_g][i_cls] = beta_L;
                      update_by_beta(i_g,i_cls,n,no_cls,alpha,
                         log_postv,sum_log_likev,sum_log_prior_betas,
                         sum_log_prior_sigmas,lv,log_sum_exp_lv,
                         sum_log_prior_betas_p,lv_p,sum_log_likev_p,
                         betas,widths_g,y_tr,cases);
                      eval_ss++;
                      if(log_postv[0] <= z) break;
                      beta_L -= w_ss;
                   }
                   /*stepping right*/
                   /*printf("stepping right\n");*/
                   while((K--) > 0){
                      betas[i_g][i_cls] = beta_R;
                      update_by_beta(i_g,i_cls,n,no_cls,alpha,
                         log_postv,sum_log_likev,sum_log_prior_betas,
                         sum_log_prior_sigmas,lv,log_sum_exp_lv,
                         sum_log_prior_betas_p,lv_p,sum_log_likev_p,
                         betas,widths_g,y_tr,cases);
                      eval_ss++;
                      if(log_postv[0] <= z ) break;
                      beta_R += w_ss;
                   }
                   //drawing new points
                   while(1){
                      //Rprintf("I am sampling for beta\n");
                      GetRNGstate();
                      betas[i_g][i_cls] = runif(beta_L,beta_R);
                      PutRNGstate();
                      update_by_beta(i_g,i_cls,n,no_cls,alpha,
                         log_postv,sum_log_likev,sum_log_prior_betas,
                         sum_log_prior_sigmas,lv,log_sum_exp_lv,
                         sum_log_prior_betas_p,lv_p,sum_log_likev_p,
                         betas,widths_g,y_tr,cases);
                      eval_ss++;
                      if(log_postv[0] >= z) break;
                      //Rprintf("draw beta ");
                      //Rprintf("z is %f, logpost is %f\n",z,log_postv[0]);
                      if(betas[i_g][i_cls] > beta_o) beta_R = betas[i_g][i_cls];
                      else beta_L = betas[i_g][i_cls]; 
                      //avoid infinit loop
                      if(beta_R-beta_L <= 1e-10) break;
                   }
               }
           }
            
	  for(i_hp = 0;i_hp < iters_hp; i_hp++) {
             //Rprintf("I am starting updating sigma\n");
             for(i_sgm = 0; i_sgm < order+1; i_sgm++){
 		 if(a_sigmas[i_sgm] > 1e-10){
                      GetRNGstate();
                      U = runif(0,1);
                      PutRNGstate();
                      if(U == 1) break;
                      z = log_postv[0] + log(U);
                      /*  placing the first interval */
                      GetRNGstate();
                      U = runif(0,1);
                      PutRNGstate();
                      sigma_L = log_sigmas[i_sgm] - (w_hp) * U;
                      sigma_R = sigma_L + (w_hp);
                        
                      GetRNGstate();
                      U = runif(0,1);
                      PutRNGstate();
                      J = floor((m_hp) * U);
                      K = (m_hp-1) - J;

                      log_sigma_o = log_sigmas[i_sgm]; 
                      //remove terms in sum_log_likev involving this sigma
                      sum_log_prior_sigmas_p = sum_log_prior_sigmas[0]-
                                               log_prior_sigmas[i_sgm];
                      sum_log_prior_betas_p = sum_log_prior_betas[0];
                      for(i_g = 0; i_g < no_g; i_g++){
                          if(*(sigmas_g + i_g * (order+1) + i_sgm) > 0){
                             widths_g_p[i_g] = widths_g[i_g] - 
                                  exp(alpha * log_sigmas[i_sgm])*
                                  (*(sigmas_g + i_g * (order+1) + i_sgm)); 
                             if(widths_g_p[i_g] < 0) widths_g_p[i_g] = 0;
                             for(i_cls = cls_begin; i_cls < no_cls; i_cls++){
                               sum_log_prior_betas_p -= 
                                  dstable(betas[i_g][i_cls],0,
                                          widths_g[i_g],alpha); 
                             }
                          }
                      }
                                              
                      /*stepping left*/
                      while((J--) > 0){
                        //if(sigma_L < -10) {
                        //   sigma_L = -10;
                        //   break;
                        //}
                        log_sigmas[i_sgm] = sigma_L;
                        update_by_sigma(i_sgm, log_sigmas,widths_g,
                          log_prior_sigmas, sum_log_prior_sigmas,
                          sum_log_prior_betas,sum_log_likev,log_postv,
                          sum_log_prior_sigmas_p,widths_g_p,
                          sum_log_prior_betas_p,
                          no_cls,no_g,order,sigmas_g, 
                          a_sigmas,b_sigmas,alpha,cls_begin,betas);
                        eval_hp++;
                        if(log_postv[0] <= z) break;
                        sigma_L -= w_hp;
                      }
                      /*stepping right*/
                      /*printf("stepping right\n");*/
                      while(K-- > 0){
                        //if(sigma_R > 10) {
                        //   sigma_R = 10;
                        //   break;
                        //}
                        log_sigmas[i_sgm] = sigma_R;
                        update_by_sigma(i_sgm,log_sigmas,widths_g,
                          log_prior_sigmas, sum_log_prior_sigmas,
                          sum_log_prior_betas,sum_log_likev,log_postv,
                          sum_log_prior_sigmas_p,widths_g_p,
                          sum_log_prior_betas_p,
                          no_cls,no_g,order,sigmas_g, 
                          a_sigmas,b_sigmas,alpha,cls_begin,betas);
                        eval_hp++;
                        if(log_postv[0] <= z ) break;
                        sigma_R += w_hp;
                       }
                       
                       while(1){
                         //Rprintf("I am sampling for beta\n");
                         GetRNGstate();
                         log_sigmas[i_sgm] = runif(sigma_L,sigma_R);
                         PutRNGstate();
                         update_by_sigma(i_sgm, log_sigmas,widths_g,
                          log_prior_sigmas, sum_log_prior_sigmas,
                          sum_log_prior_betas,sum_log_likev,log_postv,
                          sum_log_prior_sigmas_p,widths_g_p,
                          sum_log_prior_betas_p,
                          no_cls,no_g,order,sigmas_g, 
                          a_sigmas,b_sigmas,alpha,cls_begin,betas);
                         eval_hp++;
                         if(log_postv[0] >= z) break;
                         if(log_sigmas[i_sgm] > log_sigma_o) 
                             sigma_R = log_sigmas[i_sgm];
                         else sigma_L = log_sigmas[i_sgm];
                         //avoid infinit loop
                         if(sigma_R-sigma_L <= 1e-10) break;
                      }
                 }
             }
          }
       }
      
      /*Rprintf("write markov chain state to a file\n");
      /*write markov chain state to a file*/
       eval_ss /= (no_cls - cls_begin)*no_g*iters_ss;
       eval_hp /= no_hp*iters_hp*iters_ss;
       //writing markov chain iterations to file 
       fp_mc = fopen_chk(mc_file,"ab");
       fwrite_chk(betas,sizeof(double),no_g * no_cls,fp_mc,mc_file);
       //printf("\n\niter %d, beta3103 = %f, beta3538 = %f\n\n",
       //       next_iter_mc,*(&betas[0][0]+3103),*(&betas[0][0]+3538));
       fwrite_chk(log_sigmas,sizeof(double), order + 1,fp_mc,mc_file);
       //print_double_array(order+1,log_sigmas,order+1);
       fwrite_chk(sum_log_likev,sizeof(double),1,fp_mc,mc_file);
       fwrite_chk(sum_log_prior_betas,sizeof(double),1,fp_mc,mc_file);
       fwrite_chk(sum_log_prior_sigmas,sizeof(double),1,fp_mc,mc_file);
       fwrite_chk(log_postv,sizeof(double),1,fp_mc,mc_file);
       fwrite_chk(&eval_ss, sizeof(double),1,fp_mc,mc_file);
       fwrite_chk(&eval_hp, sizeof(double),1,fp_mc,mc_file);
       fwrite_chk(widths_g, sizeof(double),no_g,fp_mc,mc_file);
       fwrite_chk(&next_iter_mc,sizeof(int),1,fp_mc,mc_file);
       fclose(fp_mc);
    }   
    Rprintf("There are %d of samples in file '%s'.\n", next_iter_mc,mc_file);
}
/*****************************************************************************/
int_vec* find_cases_gid(int_list *cases_g, int gid){

    int_node *node;
    node = cases_g->next;
    if( gid >= cases_g -> ll){
       Rprintf("You are searching a node not existing");
       exit(1);
    }
    int id_node = 0;
    while(id_node++ < gid & node != 0){
        node = node -> next;
    }
    return(node->value);
}

/*****************************************************************************/
/* NOTE: this function assumes that the response, y_tr, is
  coded by 0,...,no_cls - 1  */
void update_by_beta
   (int i_g,int i_cls,int n,int no_cls,int alpha,
    double log_postv[1],double sum_log_likev[1],double sum_log_prior_betas[1],
    double sum_log_prior_sigmas[1], 
    double lv[][no_cls],double log_sum_exp_lv[],
    double sum_log_prior_betas_p,double lv_p[][no_cls],double sum_log_likev_p,
    double betas[][no_cls],double widths_g[],int y_tr[],int_vec* cases)
{
    int i,i_case;
    // updating sum_log_likev 
    sum_log_likev[0] = sum_log_likev_p;
    for(i = 0; i < cases->lv; i++){
       i_case = cases->v[i];
       //updating lv
       lv[i_case][i_cls] = lv_p[i_case][i_cls] + betas[i_g][i_cls];
       log_sum_exp_lv[i_case] = log_sum_exp(no_cls,&lv[i_case][0]);
       //adding new value
       if(i_cls == y_tr[i_case]) 
          sum_log_likev[0] += lv[i_case][i_cls]; 
       sum_log_likev[0] -= log_sum_exp_lv[i_case]; 
    } 
    //updating prior of betas   
    sum_log_prior_betas[0] = sum_log_prior_betas_p +
        dstable(betas[i_g][i_cls],0,widths_g[i_g],alpha);

    log_postv[0] = sum_log_likev[0] + sum_log_prior_betas[0] + 
                   sum_log_prior_sigmas[0]; 
}

/*****************************************************************************/
void update_by_sigma
    (int i_sgm, double log_sigmas[], double widths_g[],
     double log_prior_sigmas[], double sum_log_prior_sigmas[1],
     double sum_log_prior_betas[1],double sum_log_likev[1],double log_postv[1],
     double sum_log_prior_sigmas_p,double widths_g_p[],
     double sum_log_prior_betas_p,
     int no_cls,int no_g, int order, int sigmas_g[][order+1], 
     double a_sigmas[],double b_sigmas[],int alpha,int cls_begin,
     double betas[no_g][no_cls]
    )
// widths_g_p is the width for each group without adding sigma[i_sgm]  
//sum_log_prior_betas_p is without the log prior of those betas involving 
//sigma[i_sgm]  
//sum_log_prior_sigmas_p is without the log prior of sigma[i_sgm]
{
    int i_g, i_cls;
    sum_log_prior_betas[0] = sum_log_prior_betas_p;
    for(i_g = 0; i_g < no_g; i_g++){
        if(sigmas_g[i_g][i_sgm] > 0){
            widths_g[i_g] = widths_g_p[i_g] + sigmas_g[i_g][i_sgm] * 
                            exp(alpha * log_sigmas[i_sgm]); 
            if(widths_g[i_g] <= 0 ){
                log_postv[0] = log(0);
                return;
            }
            
            for(i_cls = cls_begin; i_cls < no_cls; i_cls++){
                sum_log_prior_betas[0] += 
                   dstable(betas[i_g][i_cls],0,widths_g[i_g],alpha); 
            }
        }
    }

    log_prior_sigmas[i_sgm] =
        dnorm(log_sigmas[i_sgm],b_sigmas[i_sgm],a_sigmas[i_sgm],1);
    sum_log_prior_sigmas[0] = sum_log_prior_sigmas_p + log_prior_sigmas[i_sgm];

    log_postv[0] = sum_log_likev[0] + sum_log_prior_betas[0] +
                   sum_log_prior_sigmas[0];
}


/*****************************************************************************/

void initialize_mc_state
     (int alpha, int cls_begin,int no_cls,int order,int n,
      int no_g, double widths_g[], double log_sigmas[],
      double betas[no_g][no_cls], int y_tr[],
      int_list* cases_g,
      double a_sigmas[],double b_sigmas[],
      int sigmas_g[][order+1],
      double sum_log_prior_betas[],
      double log_prior_sigmas[],
      double sum_log_prior_sigmas[],
      double log_postv[],
      double lv[][no_cls],
      double sum_log_likev[],
      double log_sum_exp_lv[]
    )
{
    int i_sgm,i_cls,i_g,i_case,i;
    int_vec* cases;
    
    set_double_array(order+1,log_prior_sigmas,0);
    for(i_sgm = 0; i_sgm < order + 1; i_sgm++){
        if(a_sigmas[i_sgm] > 1e-10)
           log_prior_sigmas[i_sgm] =
              dnorm(log_sigmas[i_sgm],b_sigmas[i_sgm],a_sigmas[i_sgm],1);
    }
    sum_log_prior_sigmas[0] = sum_double(order+1,log_prior_sigmas);

    set_double_array(no_g,widths_g,0);
    for(i_g = 0; i_g < no_g; i_g++){
       for(i_sgm = 0; i_sgm < order + 1; i_sgm ++){
         if(sigmas_g[i_g][i_sgm] > 0)
           widths_g[i_g] += sigmas_g[i_g][i_sgm] *
                            exp(alpha*log_sigmas[i_sgm]);
       }
    }
    set_double_array(n*no_cls,lv,0.0);
    sum_log_prior_betas[0] = 0;
    for(i_g = 0; i_g < no_g; i_g++){
      cases = find_cases_gid(cases_g, i_g);
      for(i_cls = cls_begin; i_cls < no_cls; i_cls ++){
         sum_log_prior_betas[0] += 
              dstable(betas[i_g][i_cls],0,widths_g[i_g],alpha);
         for(i = 0; i < cases->lv; i++)
         {  i_case = cases->v[i];
            lv[i_case][i_cls] += betas[i_g][i_cls];
         }
      }
    }
    sum_log_likev[0] = 0;
    for(i_case = 0 ; i_case < n; i_case++){
       log_sum_exp_lv[i_case] = log_sum_exp(no_cls, &lv[i_case][0]);
       sum_log_likev[0] += lv[i_case][y_tr[i_case]] - log_sum_exp_lv[i_case];
    }
    
    log_postv[0] = sum_log_likev[0] + sum_log_prior_betas[0] +
                   sum_log_prior_sigmas[0];
}

