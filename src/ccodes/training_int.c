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
      double log_sigmas[])
{   /*srand(time(0)); */
    int n,p,no_g,sequence,order,no_hp;
    int_list* cases_g = gen_int_list(NULL);
    int_list* ptns_g = gen_int_list(NULL);
    int* sigmas_g;
    
    /*loading compression information*/
    load_compress(cases_g,ptns_g,&sequence,&order,&n,&p,&sigmas_g,
                  cases_ptns_file,1);
    no_g = cases_g -> ll;
    /*requesting space for markov chaing sampling*/
    double
         betas[no_g][no_cls],widths_g[no_g],widths_g_o[no_g],lv[n][no_cls],
         log_sum_exp_lv[n],sum_log_likev[1],sum_log_prior_betas[1],
         log_prior_sigmas[order+1],sum_log_prior_sigmas[1],
	 log_postv[1],log_postv_o;

    int i_mc, i_ss, i_hp, i_g, i_cls, i_sgm;
    double beta_o,log_sigma_o,U,z,beta_n,log_sigma_n,beta_L,beta_R,
           log_sigma_L,log_sigma_R,eval_ss,rej_hp;
    long size_each_iter=sizeof(int)+sizeof(double)*(no_g*(no_cls+1)+order+7);
    int J,K;
    FILE *fp_mc;
    int next_iter_mc = 0;
    int cls_begin = 0;
    if(no_cls == 2) cls_begin = 1;
    
    //count how many hyperparameters to be updated with MCMC
    no_hp = 0;
    for(i_sgm = 0; i_sgm < order + 1;i_sgm++) {
        if(a_sigmas[i_sgm] < 1e10) no_hp++;
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
    i_mc = 0;
    while(i_mc ++ < iters_mc) {
       next_iter_mc++;
       eval_ss = 0;
       rej_hp = 0;
       i_ss = 0;
       //updating betas
       while(i_ss++ < iters_ss) {
            //Rprintf("I am starting updating betas\n");
            for(i_g = 0; i_g < no_g; i_g++) {
                for(i_cls = cls_begin; i_cls < no_cls; i_cls++){
                   beta_o = betas[i_g][i_cls];
                   GetRNGstate();
                   z =  log_postv[0] + log(runif(0,1));
                   PutRNGstate();

                   /* stepping out */
                   GetRNGstate();
                   U = runif(0,1);
                   PutRNGstate();
                   beta_R = beta_o - (w_ss) * U;
                   beta_L = beta_R + (w_ss);
                   GetRNGstate();
                   U = runif(0,1);
                   PutRNGstate();
                   J = floor((m_ss) * U);
                   K = (m_ss) - J;
                   /*stepping left*/
                   /*printf("stepping left\n");*/
                   do{
                      beta_L -= w_ss;
                      update_by_beta(beta_L,beta_o,i_g,i_cls,
                        n,no_cls,lv,log_sum_exp_lv,log_postv,
                        sum_log_prior_betas,sum_log_likev,
                        widths_g,y_tr,cases_g,alpha);
                      beta_o = beta_L;
                      J --;
                      eval_ss++;
                   }  while(log_postv[0] > z & J > 0);
                   /*stepping right*/
                   /*printf("stepping right\n");*/

                   do{
                      beta_R += w_ss;
                      update_by_beta(beta_R,beta_o,i_g,i_cls,
                          n,no_cls,lv,log_sum_exp_lv,log_postv,
                          sum_log_prior_betas,sum_log_likev,
                          widths_g,y_tr,cases_g,alpha);
                      beta_o = beta_R;
                      K --;
                      eval_ss++;
                   }  while(log_postv[0] > z & K > 0);
                   /*drawing new points*//*
                   Rprintf("drawing new points\n");
		   Rprintf("betaL %f betaR %f\n",beta_L,beta_R);*/
                   do{
                      //Rprintf("I am sampling for beta\n");
                      GetRNGstate();
                      beta_n = runif(beta_L,beta_R);
                      PutRNGstate();
                      update_by_beta(beta_n,beta_o,i_g,i_cls,
                          n,no_cls,lv,log_sum_exp_lv,log_postv,
                          sum_log_prior_betas,sum_log_likev,
                          widths_g,y_tr,cases_g,alpha);
                      beta_o = beta_n;

		      eval_ss++;
                      if(beta_n > betas[i_g][i_cls]) beta_R = beta_n;
                      else beta_L = beta_n;
                    } while(log_postv[0] < z & fabs(beta_L-beta_R) > 1e-10);
		   /*Rprintf("betaL %f betaR %f betaN %f\n",beta_L,
		   beta_R,beta_n);
                   */
		   betas[i_g][i_cls] = beta_n;
               }
            }

            /*updating hyperparameters*/
            /*Rprintf("updating hyperparameters\n");
            */
            //updating hyperparmeters
	    i_hp = 0;
            while(i_hp++ < iters_hp) {
                //Rprintf("I am starting updating sigma\n");
                for(i_sgm = 0; i_sgm < order+1; i_sgm++){
		   if(a_sigmas[i_sgm] < 1e10){
		       update_by_sigma(i_sgm,log_sigmas,
                             widths_g,widths_g_o,w_hp,no_cls,no_g,order,
			     sigmas_g,sum_log_likev,sum_log_prior_betas,
                             log_prior_sigmas,sum_log_prior_sigmas,
                             log_postv,a_sigmas,b_sigmas,betas,alpha,
                             cls_begin, &rej_hp);
                   }
                }
            }

      }
      /*Rprintf("write markov chain state to a file\n");
      /*write markov chain state to a file*/
      eval_ss /= (no_cls - cls_begin)*no_g*iters_ss;
      rej_hp /= no_hp*iters_hp*iters_ss;
      //writing markov chain iterations to file 
      fp_mc = fopen_chk(mc_file,"ab");
      fwrite_chk(betas,sizeof(double),no_g * no_cls,fp_mc,mc_file);
      fwrite_chk(log_sigmas,sizeof(double), order + 1,fp_mc,mc_file);
      fwrite_chk(sum_log_likev,sizeof(double),1,fp_mc,mc_file);
      fwrite_chk(sum_log_prior_betas,sizeof(double),1,fp_mc,mc_file);
      fwrite_chk(sum_log_prior_sigmas,sizeof(double),1,fp_mc,mc_file);
      fwrite_chk(log_postv,sizeof(double),1,fp_mc,mc_file);
      fwrite_chk(&eval_ss, sizeof(double),1,fp_mc,mc_file);
      fwrite_chk(&rej_hp, sizeof(double),1,fp_mc,mc_file);
      fwrite_chk(widths_g, sizeof(double),no_g,fp_mc,mc_file);
      fwrite_chk(&next_iter_mc,sizeof(int),1,fp_mc,mc_file);
      fclose(fp_mc);
    }

    Rprintf("Markov chain sampling is done; %d iterations exist in file %s\n",
             next_iter_mc,mc_file);
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
void update_by_beta(double beta_n,double beta_o,int i_g,int i_cls,
                    int n,int no_cls,double lv[][no_cls],
                    double log_sum_exp_lv[],double log_postv[],
                    double sum_log_prior_betas[],
                    double sum_log_likev[],
                    double widths_g[],int y_tr[],
                    int_list* cases_g,int alpha)
{
    int i,i_case;
    double lv_new,log_likev_new,log_sum_exp_lv_new,diff;
    int_vec *cases;
    /*printf("old beta %f\n",beta_o);
    printf("new beta %f\n",beta_n);*/

    cases = find_cases_gid(cases_g, i_g);
    log_postv[0] -= sum_log_likev[0];

    for(i = 0; i < cases->lv; i++){
       i_case = cases->v[i];
       lv_new = lv[i_case][i_cls] - beta_o + beta_n;
       diff = lv_new - lv[i_case][i_cls];

       lv[i_case][i_cls] = lv_new;
       if(i_cls == y_tr[i_case]){
          sum_log_likev[0] += diff;
       }
       /*printf("sumloglikev %f\n",sum_log_likev[0]);
         */
       log_sum_exp_lv_new = log_sum_exp(no_cls,&lv[i_case][0]);
       diff = log_sum_exp_lv_new - log_sum_exp_lv[i_case];
       log_sum_exp_lv[i_case] = log_sum_exp_lv_new;
       sum_log_likev[0] -= diff;
       /*printf("sumloglikev %f\n",sum_log_likev[0]);
         */
    }
    log_postv[0] += sum_log_likev[0];

    if(alpha == 1)
        diff = dcauchy(beta_n,0,widths_g[i_g],1)
              - dcauchy(beta_o,0,widths_g[i_g],1);
    if(alpha == 2)
        diff = dnorm(beta_n,0,sqrt(widths_g[i_g]),1)
              - dnorm(beta_o,0,sqrt(widths_g[i_g]),1);
    sum_log_prior_betas[0] += diff;

    log_postv[0] += diff;
    /*some debuging statement*/
    /*printf("lv new %f\n",lv_new);
         */
    /*log_sum_exp_lv_new = log_minus_add_exp(log_sum_exp_lv[i_case],
             lv[i_case][i_cls],lv_new);*/
    /*printf("logsumexplvnew %f\n",log_sum_exp_lv_new);*/

}

/*****************************************************************************/
void update_by_sigma
     (int i_sgm,double log_sigmas[],double widths_g[],double widths_g_o[],
      double w_hp,
      int no_cls,int no_g,int order,int sigmas_g[][order+1],
      double *sum_log_likev, double *sum_log_prior_betas,
      double log_prior_sigmas[],double sum_log_prior_sigmas[],
      double *log_postv,double a_sigmas[],double b_sigmas[],
      double betas[][no_cls],int alpha,int cls_begin, int rej_hp[]
     )
{
    int i_g, i_cls;
    double log_sigma_o,log_prior_sigma_o,sum_log_prior_betas_o,log_postv_o,
           log_sigma_n,sum_log_prior_sigmas_o;
/*    Rprintf("sigma_n %f, b_sigma %f a_sigma %f\n",
                 sigma_n,b_sigmas[i_sgm],a_sigmas[i_sgm]);
*/
    log_sigma_n = rnorm(log_sigmas[i_sgm],w_hp);

    if(log_sigma_n < -10 || log_sigma_n > 10){ 
    	rej_hp[0]++;
    	return;
    }
    else
    {
	//save the old values
	for(i_g = 0; i_g < no_g; i_g++){
		if(sigmas_g[i_g][i_sgm] > 0){
		   widths_g_o[i_g] = widths_g[i_g];
		}
	}
	log_prior_sigma_o = log_prior_sigmas[i_sgm];
	sum_log_prior_sigmas_o = sum_log_prior_sigmas[0];
	sum_log_prior_betas_o = sum_log_prior_betas[0];
	log_postv_o = log_postv[0];
	log_sigma_o = log_sigmas[i_sgm];

	//start updating with new sigma
	log_sigmas[i_sgm] = log_sigma_n;
	log_prior_sigmas[i_sgm] =
	     dnorm(log_sigmas[i_sgm],b_sigmas[i_sgm],1.0/a_sigmas[i_sgm],1);
	sum_log_prior_sigmas[0] += log_prior_sigmas[i_sgm] - log_prior_sigma_o;

	for(i_g = 0; i_g < no_g; i_g++){
	    if(sigmas_g[i_g][i_sgm] > 0){
		widths_g[i_g] += sigmas_g[i_g][i_sgm] *
				(exp(alpha*log_sigmas[i_sgm]) -
				 exp(alpha*log_sigma_o));
		if(widths_g[i_g] < 0) { 
		    widths_g[i_g] = 0;
		    Rprintf("Note: width became 0\n");
		}

		if(alpha == 1){
		   for(i_cls = cls_begin; i_cls < no_cls; i_cls++)
		      sum_log_prior_betas[0] +=
			dcauchy(betas[i_g][i_cls],0,widths_g[i_g],1)
			- dcauchy(betas[i_g][i_cls],0,widths_g_o[i_g],1);
		}
		if(alpha == 2){
		   for(i_cls = cls_begin; i_cls < no_cls; i_cls++)
		      sum_log_prior_betas[0] +=
			dnorm(betas[i_g][i_cls],0,sqrt(widths_g[i_g]),1)
			- dnorm(betas[i_g][i_cls],0,sqrt(widths_g_o[i_g]),1);
		}
	   }
	}
	log_postv[0] = sum_log_likev[0] + sum_log_prior_betas[0] +
		       sum_log_prior_sigmas[0];
	if(log(runif(0,1)) >= log_postv[0] - log_postv_o){
	       //rejecting the proposed value and retrieve everything
	       //retrieve the old values
	       rej_hp[0]++;
	       for(i_g = 0; i_g < no_g; i_g++){
		  if(sigmas_g[i_g][i_sgm] > 0){
		      widths_g[i_g] = widths_g_o[i_g];
		  }
	       }
	       log_prior_sigmas[i_sgm] = log_prior_sigma_o;
	       sum_log_prior_sigmas[0] = sum_log_prior_sigmas_o;
	       sum_log_prior_betas[0] = sum_log_prior_betas_o;
	       log_postv[0] = log_postv_o;
	       log_sigmas[i_sgm] = log_sigma_o;
	}
   }
}
/*****************************************************************************/

void initialize_mc_state
     (int alpha, int cls_begin,int no_cls,int order,int n,
      int no_g, double widths_g[], double log_sigmas[],
      double betas[no_g][no_cls], double y_tr[],
      int_list* cases_g,
      double a_sigmas[],double b_sigmas[],
      int sigmas_g[][order+1],
      double sum_log_prior_betas[],
      double log_prior_sigmas[],
      double sum_log_prior_sigmas[],
      double log_postv[],
      double lv[][no_cls+1],
      double sum_log_likev[],
      double log_sum_exp_lv[]
    )
{
    int i_sgm,i_cls,i_g;
    sum_log_prior_betas[0] = 0;
    set_double_array(order+1,log_prior_sigmas,0);
    for(i_sgm = 0; i_sgm < order + 1; i_sgm++){
        if(a_sigmas[i_sgm] < 1e10)
           log_prior_sigmas[i_sgm] =
              dnorm(log_sigmas[i_sgm],b_sigmas[i_sgm],1.0/a_sigmas[i_sgm],1);
    }
    sum_log_prior_sigmas[0] = sum_double(order+1,log_prior_sigmas);
    /*print_double_array(p+1,log_prior_sigmas,p+1);
    */

    set_double_array(no_g,widths_g,0);
    for(i_g = 0; i_g < no_g; i_g++){
       for(i_sgm = 0; i_sgm < order + 1; i_sgm ++){
         if(sigmas_g[i_g][i_sgm] > 0)
           widths_g[i_g] += sigmas_g[i_g][i_sgm] *
                            exp(alpha*log_sigmas[i_sgm]);
       }
    }

    for(i_g = 0; i_g < no_g; i_g++){
      if(alpha == 1)
         sum_log_prior_betas[0] +=
            (no_cls-cls_begin)*dcauchy(0,0,widths_g[i_g],1);
      if(alpha == 2)
         sum_log_prior_betas[0] +=
            (no_cls-cls_begin)*dnorm(0,0,R_pow(widths_g[i_g],0.5),1);
    }
    /*print_double_array(p+1,sigmas,p+1);
    */
    set_double_array(n*no_cls,lv,0.0);
    set_double_array(n,log_sum_exp_lv, log(no_cls));
    sum_log_likev[0] = - n * log_sum_exp_lv[0];
    log_postv[0] = sum_log_likev[0] + sum_log_prior_sigmas[0] +
                   sum_log_prior_betas[0];
    /*Rprintf("logpriorsigma %f, loglike %f\n",
    sum_log_prior_sigmas[0],sum_log_likev[0]);
    */
    for(i_g = 0; i_g < no_g; i_g ++)
      for(i_cls = cls_begin; i_cls < no_cls; i_cls ++)
         update_by_beta(betas[i_g][i_cls],0.0,i_g,i_cls,
            n,no_cls,lv,log_sum_exp_lv,log_postv,
            sum_log_prior_betas,sum_log_likev,widths_g,y_tr,
            cases_g,alpha);
}


/*

log_sigma_o = log_sigmas[i_sgm];
                     GetRNGstate();
                     z = log_postv[0] + log(runif(0,1));
                     PutRNGstate();
                     // stepping out
                     GetRNGstate();
                     U = runif(0,1);
                     PutRNGstate();
                     log_sigma_R = log_sigma_o - (w_hp) * U;
                     log_sigma_L = log_sigma_R + (w_hp);
                     GetRNGstate();
                     U = runif(0,1);
                     PutRNGstate();
                     J = floor((m_hp) * U);
                     K = (m_hp)-J;
                     do{
                        log_sigma_L -= w_hp;
                        update_by_sigma(log_sigma_L,log_sigma_o,i_sgm,
                             widths_g,no_cls,no_g,order,sigmas_g,
                             sum_log_likev,sum_log_prior_betas,
                             log_prior_sigmas,sum_log_prior_sigmas,
                             log_postv,a_sigmas,b_sigmas,betas,alpha,
                             cls_begin);
                             log_sigma_o = log_sigma_L;
                             J --;
                             eval_hp ++;
                     } while (log_postv[0] > z & J > 0);
                     do{

                        log_sigma_R += w_hp;
                        update_by_sigma(log_sigma_R,log_sigma_o,i_sgm,
                             widths_g,no_cls,no_g,order,sigmas_g,
                             sum_log_likev,sum_log_prior_betas,
                             log_prior_sigmas,sum_log_prior_sigmas,
                             log_postv,a_sigmas,b_sigmas,betas,alpha,
                             cls_begin);
                        log_sigma_o = log_sigma_R;
                        K --;
                        eval_hp ++;
                     } while (log_postv[0] > z & K > 0);
		     do{
		         Rprintf("I am sampling sigmas\n");
                         GetRNGstate();
                         log_sigma_n = runif(log_sigma_L,log_sigma_R);
                         PutRNGstate();
                         update_by_sigma(log_sigma_n,log_sigma_o,i_sgm,
                             widths_g,no_cls,no_g,order,sigmas_g,
                             sum_log_likev,sum_log_prior_betas,
                             log_prior_sigmas,sum_log_prior_sigmas,
                             log_postv,a_sigmas,b_sigmas,betas,alpha,
                             cls_begin);
                             log_sigma_o = log_sigma_n;
                             eval_hp ++;

                        if(log_sigma_n > log_sigmas[i_sgm])
                            log_sigma_R = log_sigma_n;
                        else log_sigma_L = log_sigma_n;
                     } while (log_postv[0] < z  &
                                 fabs( log_sigma_L - log_sigma_R) > 1e-10);
		     log_sigmas[i_sgm] = log_sigma_n;


 */