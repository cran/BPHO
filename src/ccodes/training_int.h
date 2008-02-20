void R_training
     (char** mc_file,char** cases_ptns_file, int y_tr[],
      int no_cls[], int iters_mc[],int iters_ss[], int iters_hp[],
      double w_ss[], double w_hp[], int m_ss[], int m_hp[],
      int alpha[],double a_sigmas[], double b_sigmas[],
      double log_sigmas[]);

void training
     (char mc_file[],char cases_ptns_file[], int y_tr[],
      int no_cls, int iters_mc,int iters_ss, int iters_hp,
      double w_ss, double w_hp, int m_ss, int m_hp,
      int alpha, double a_sigmas[], double b_sigmas[],
      double log_sigmas[]);

int_vec* find_cases_gid(int_list *cases_g, int gid);

void update_by_beta(double beta_n,double beta_o,int i_g,int i_cls,
                    int n,int no_cls,
                    double lv[][no_cls],double log_sum_exp_lv[],
		    double log_postv[], double sum_log_prior_betas[],
		    double sum_log_likev[],
		    double widths_g[],int y_tr[],
		    int_list* cases_g,int alpha);

void initialize_mc_state
     (int alpha, int cls_begin,int no_cls,int depth,int n,
      int no_g, double widths_g[], double log_sigmas[],
      double betas[no_g][no_cls], double y_tr[],
      int_list* cases_g,
      double a_sigmas[],double b_sigmas[],
      int sigmas_g[][depth+1],
      double sum_log_prior_betas[],
      double log_prior_sigmas[],
      double sum_log_prior_sigmas[],
      double log_postv[],
      double lv[][no_cls+1],
      double sum_log_likev[],
      double log_sum_exp_lv[]
    );

void update_by_sigma
     (int i_sgm,double log_sigmas[],double widths_g[],double widths_g_o[],
      double w_hp,
      int no_cls,int no_g,int order,int sigmas_g[][order+1],
      double *sum_log_likev, double *sum_log_prior_betas,
      double log_prior_sigmas[],double sum_log_prior_sigmas[],
      double *log_postv,double a_sigmas[],double b_sigmas[],
      double betas[][no_cls],int alpha,int cls_begin, int rej_hp[]
      );
