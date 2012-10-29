/****************************************************************************/
void find_sigmas_test_seq
        (int_list* ptns_g,int test[],int order,int sigmas_test[][order+1],
         int no_splits[1])
{
    int_node *node;
    int *ptn,expressed,o,t;
    node = ptns_g -> next;
    int i_sgm,i_g,no_g = ptns_g -> ll;    
    set_value_array((no_g)*(order+1),&sigmas_test[0][0],0);
                                          
    set_value_array(order+1,&sigmas_test[no_g][0],1); 
    i_g = 0;
    while(node != NULL){
       ptn = node -> value -> v;
       expressed = 1;
       for(o = 0; o < ptn[order]; o ++){
          if(ptn[o] != test[o]){
              expressed = 0;
              break;
          }
       } 
       if(expressed==1){
         sigmas_test[i_g][ptn[order]]++;
         sigmas_test[no_g][ptn[order]]--;
         for(t = ptn[order]; t < order; t++){
           if(ptn[t] == 0) break;
           if(ptn[t] != test[t]) {
              no_splits[0]++;
              break;
           }
           sigmas_test[i_g][t+1]++;
           sigmas_test[no_g][t+1]--;
         }
       }
       node = node->next;
       i_g ++;
    }
}
/*****************************************************************************/
void add_sigmas_test_cls
       (int no_g,int id_g,int no_ptn,int p,int ptn[][p+1],
        int test[],int order,int sigmas_test[][order+1]){
   int no_nofixed,no_nozero,no;
   int i_ptn,i_fth,no_nozero_lack,i,exist;
   for(i_ptn = 0; i_ptn < no_ptn; i_ptn++){
      no_nozero = ptn[2*i_ptn][p];
      no_nofixed = p - ptn[2*i_ptn+1][p];
      exist = 1;
      for(i_fth = 0; i_fth < p; i_fth++){
         if(ptn[2*i_ptn+1][i_fth] == 1){
            if(ptn[2*i_ptn][i_fth] != 0 & 
                (ptn[2*i_ptn][i_fth] != test[i_fth])) {
               exist = 0; 
               break; 
            }
         }
         else{
             if(ptn[0][i_fth] != test[i_fth]){
                no_nofixed --; 
             }
         }
     }
     if(exist){
       no_nozero_lack = imin2(order - no_nozero, no_nofixed);
       for(i = 0; i <= no_nozero_lack; i++){
           no = choose(no_nofixed, i);
           sigmas_test[id_g][i+no_nozero] += no;
           sigmas_test[no_g][i+no_nozero] -= no;
       }
     }
  }
}

/***************************************************************************/
void find_sigmas_test_cls
        (int_list* ptns_g,int test[],int order,int p,
	int sigmas_test[][order+1])
{
    int_node *node;
    node = ptns_g -> next;
    int i_sgm,i_g, no_g = ptns_g -> ll;    
    set_value_array((no_g+1)*(order+1),&sigmas_test[0][0],0);      
    for(i_sgm = 0; i_sgm < order+1; i_sgm++)  
        sigmas_test[no_g][i_sgm] = choose(p,i_sgm);
    i_g = 0;
    while(node != NULL){
        add_sigmas_test_cls
	  (no_g,i_g,(node->value->lv)/(2*(p+1)),
	   p,node->value->v,test,order,sigmas_test);
        node = node->next;
        i_g ++;
    }
}
/***************************************************************************/

void find_widths_g
       (int no_g, int order, int sigmas_g[no_g][order+1],
        double log_sigmas[order+1],double widths_g[no_g], int alpha)
{   set_double_array(no_g, widths_g, 0.0);
    int i_g, i_sgm;
    for(i_g = 0; i_g < no_g; i_g ++)
        for(i_sgm = 0; i_sgm < order + 1; i_sgm ++){
            if(sigmas_g[i_g][i_sgm] > 0) 
                widths_g[i_g] += sigmas_g[i_g][i_sgm] *
		                 exp(alpha*log_sigmas[i_sgm]);    
        }
}            

/****************************************************************************/
void R_pred
        (int *R_no_test,int *R_p,int tests[][R_p[0]],
         int *R_no_cls, double prediction[][*R_no_cls],
         char** mc_files,char** cases_ptns_files, 
         int *R_iter_b,int *R_forward,int *R_iter_n)	 
{  
   int *sigmas_g,n,p,no_g,order,no_cls,alpha,no_iters,sequence;
   int i_pred,i_g,i_ts,i_cls, no_splits=0;

   int_list *cases_g = gen_int_list(NULL);
   int_list *ptns_g = gen_int_list(NULL);
   FILE *fp_mc, *fp_ptns;
   
   fp_mc = fopen_chk(mc_files[0],"rb");
   fread_chk(&no_cls,sizeof(int),1,fp_mc,mc_files[0]);
   fread_chk(&no_g,sizeof(int),1,fp_mc,mc_files[0]);
   fread_chk(&order,sizeof(int),1,fp_mc,mc_files[0]);
   fread_chk(&alpha,sizeof(int),1,fp_mc,mc_files[0]); 
   fseek(fp_mc,-sizeof(int),SEEK_END);
   fread_chk(&no_iters,sizeof(int),1,fp_mc,mc_files[0]);
   fclose(fp_mc);

   /*loading compression information*/
   load_compress(cases_g,ptns_g,&sequence,&order,&n,&p,&sigmas_g,
                 cases_ptns_files[0],1);

   /*request space for each iteration of mc*/
   int sigmas_test[no_g+1][order+1];   
   double (*betas)[no_g][no_cls],log_sigmas[R_iter_n[0]][order+1],
          widths_g[R_iter_n[0]][no_g],log_probs_pred[R_iter_n[0]][no_cls],
          betas_test[no_g][no_cls],widths_test[no_g+1];
   betas =  R_alloc (R_iter_n[0], sizeof (*betas) );	  

   int cls_begin = (no_cls > 2)?0:1;
   /*calcuate size of bytes for each iteration of mc and size in mc file head*/
   long size_head_mc = sizeof(int) * 4;
   long size_each_iter = sizeof(int)+sizeof(double)*(no_g*(no_cls+1)+order+7);
   long skip_each_pred = (R_forward[0]-1) * size_each_iter;

   // read Markov chain samples 
   fp_mc = fopen_chk(mc_files[0],"rb");
   fseek(fp_mc,size_head_mc + R_iter_b[0] * size_each_iter,SEEK_SET); 
   for(i_pred = 0; i_pred < R_iter_n[0]; i_pred ++)
   {
      fseek(fp_mc,sizeof(int),SEEK_CUR);
      fread_chk(&betas[i_pred][0][0],sizeof(double),no_g*no_cls,fp_mc,
                mc_files[0]); 
      fread_chk(&log_sigmas[i_pred][0],sizeof(double),order+1,fp_mc,mc_files[0]);
      fseek(fp_mc,6*sizeof(double),SEEK_CUR);
      fread_chk(&widths_g[i_pred][0],sizeof(double),no_g,fp_mc,mc_files[0]);
      fseek(fp_mc, skip_each_pred, SEEK_CUR);
   }
   fclose(fp_mc);

   for(i_ts = 0; i_ts < *R_no_test; i_ts++){
      set_double_array(R_iter_n[0]*no_cls,&log_probs_pred[0][0],0); 
      if(sequence == 1)
         find_sigmas_test_seq(ptns_g,&tests[i_ts][0],order,sigmas_test, 
                              &no_splits);
      else
         find_sigmas_test_cls(ptns_g,&tests[i_ts][0],order,p,sigmas_test);
      for(i_pred = 0; i_pred < R_iter_n[0]; i_pred ++)
      {
	 find_widths_g(no_g+1,order,sigmas_test,&log_sigmas[i_pred][0],
                       widths_test,alpha);
         for(i_cls = cls_begin; i_cls < no_cls; i_cls ++) {
             for(i_g = 0; i_g < no_g; i_g ++) 
	     {
                 split(&betas_test[i_g][i_cls], &betas[i_pred][i_g][i_cls],
                       &widths_test[i_g], &widths_g[i_pred][i_g], &alpha); 
                 log_probs_pred[i_pred][i_cls] += betas_test[i_g][i_cls];
             }
             GetRNGstate();
             if(alpha == 1)
                log_probs_pred[i_pred][i_cls]+=rcauchy(0,widths_test[no_g]);
             else
                log_probs_pred[i_pred][i_cls]+=rnorm(0,sqrt(widths_test[no_g]));
             PutRNGstate();
         }
         /*normalizing latent value to make probabilities*/
         norm_log_array(no_cls,&log_probs_pred[i_pred][0]);
     }
     /*average over iterations*/
     for(i_cls = 0; i_cls < no_cls; i_cls ++) {
         for(i_pred = 0; i_pred < R_iter_n[0]; i_pred ++){
	    prediction[i_ts][i_cls] += exp(log_probs_pred[i_pred][i_cls]);
         }
	 prediction[i_ts][i_cls] /= R_iter_n[0]; 
     }
   }
   if( sequence == 1)
   Rprintf("Fraction of test cases for which splitting occures is %0.4f\n",
           (no_splits+0.0)/R_no_test[0] );
}

/***************************************************************************/

void split
      (double *beta_p, double *beta_s, double *width_p, double *width_s,
       int *alpha)
{   
    double widthratio;
    int zero=0;
    if(*width_p < 1e-10) {
       *beta_p = 0;       
       return;
    }
    if(*width_s - *width_p < 1e-10 ){
       
       *beta_p = *beta_s;
       
       return;
    }
    
    if(*alpha == 2)
    {  GetRNGstate();
       widthratio = width_p[0] / width_s[0];
       beta_p[0] = rnorm( beta_s[0] * widthratio, 
                        sqrt(width_p[0]*(1-widthratio)));
       PutRNGstate();       
       return;
    }
    
    if(*alpha == 1)
    {  split_cauchy(beta_p,beta_s,width_p,width_s,&zero);       
       return;
    }
}    
    
                                     
/*************************************************************************/

void split_cauchy
     (double *x, double *s, double *sigma1, double *sigmasum, int *debug)
{
    /*define maximum search steps and the tolerable error*/
    int MAX_ITER = 50;
    int iter = 0;
    double TOL = 1e-5;
    double* sigma2 = gen_double_array(1);
    double C,r,sigma_sq_sum,sigma_sq_diff,p0,ps;
    double L,U,yL,yU,p,error,new_root,new_y, denom;
    int CUR_sign;
        
      
     sigma2[0] = sigmasum[0] - sigma1[0];
     sigma_sq_sum = sq(sigma1[0]) + sq(sigma2[0]);
     sigma_sq_diff = sq(sigma1[0]) - sq(sigma2[0]);
     denom = quan(s[0]) + 2 * sigma_sq_sum * sq(s[0]) + sq(sigma_sq_diff);
    
     if(denom > 1e-10 ){
     
    	p0 = (sq(s[0]) - sigma_sq_diff) / denom / sigma1[0];
    	ps = (sq(s[0]) + sigma_sq_diff) / denom / sigma2[0] ;       
    	C = M_PI * sigmasum[0] / sigma1[0] / sigma2[0] /
	    ( sq(s[0]) + sq(sigmasum[0]) );
    	r = s[0]/ denom;
       
    	/*draw from uniform(0,1)*/
    	GetRNGstate();
    	p = runif(0,1);
    	PutRNGstate();
    
    	/*using illinois methods to find the x[0] with F(x[0]) = U*/
    	L = -M_PI/2; 
	yL = 0.0 - p;
    	U = M_PI/2;
	yU = 1.0 - p;
	CUR_sign = 1;
    	do{
       	   new_root = (L*yU - U*yL)/(yU-yL);
       	   x[0] = tan(new_root);      
       	   new_y =  /*compute the CDF at x[0]*/
             1.0/C * (r*log( sq(x[0]) + sq(sigma1[0]) ) - 
	              r*log( sq(x[0]-s[0])+ sq(sigma2[0]) ) + 
		      p0*(atan(x[0]/sigma1[0]) + M_PI/2) + 
		      ps*(atan((x[0]-s[0])/sigma2[0]) + M_PI/2)
		      ) 
		  - p;
	   if(*debug) Rprintf("iter=%d,new_y=%0.7f\n",iter,new_y);
	    
	   if(new_y > 0){
	      U = new_root;
	      yU = new_y;
	      if(CUR_sign == 1) yL /= 2.0;
	      CUR_sign = 1;
	   }
	   else {
	      L = new_root;
	      yL = new_y;
	      if(CUR_sign == -1) yU /= 2.0;
	      CUR_sign = -1;
	   }
           iter++;
    	} 
    	while(fabs(new_y) > TOL & iter < MAX_ITER);
     }	
     else{
        GetRNGstate();
        x[0] = rt(3) * sigma1[0]/sqrt(3);
	PutRNGstate();
     }
}


