
void find_widths_g
       (int no_g, int p, int sigmas_g[no_g][p+1],
        double sigmas[p+1],double widths_g[no_g], int alpha);

void split_cauchy
     (double *x, double *s, double *sigma1, double *sigmasum, int *debug);        

void R_pred
        (int *R_no_test, int *R_p,int tests[][R_p[0]],
         int *R_no_cls, double prediction[][*R_no_cls],
         char** mc_files,char** cases_ptns_files, 
         int *R_iter_b,int *R_forward,int *R_iter_n);                       
void split
      (double *beta_p, double *beta_s, double *width_p, double *width_s,
       int *alpha);        
void find_sigmas_test_cls
        (int_list* ptns_g,int test[],int order,int p,
	int sigmas_test[][order+1]);

void add_sigmas_test_cls
       (int no_g,int id_g,int no_ptn,int p,int ptn[][p+1],
        int test[],int order,int sigmas_test[][order+1]);	
void find_sigmas_test_seq
        (int_list* ptns_g,int test[],int order,int sigmas_test[][order+1],
         int no_splits[1]);
