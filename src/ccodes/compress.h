void compress_cls(int order,int n, int p, int features[n][p],
                  int* nos_fth,int no_cases_ign, char* file, int quiet,
		  int do_comp);


void diverge_cls(int n, int p,int features[n][p], 
            int_vec* m_cases,int_vec* m_fths,int* nos_fth, 
	    int_list* cases_g, int_list_list* ptns_g, int order, 
	    int* m_ptn, int_vec* NOs_class, int no_cases_ign,
	    int do_comp);
	    
double calc_entropy(int_vec* NOs_class, int N);

	    
double_vec* entr_sub_sub(int n, int p, int features[n][p], 
                      int_vec* sub_cases,int_vec* sub_fths, 
                      int_vec* NOs_class);
		      
void trim_ptns(int_list *cases_g, int_list_list *ptns_g);

void add_sigmas_patterns(int id_g,int no_g, int order, int p, 
                        int sigmas_g[][order+1], 
			int pattern[][p+1]);		      

void find_sigmas_interact(int_list_list* ptns_g,int order,int p,
                  int sigmas_g[][order+1]);
		  
/****************************************************************************/
		  
void compress_seq(int n, int p, int features[n][p],int* nos_fth,
                  int no_cases_ign, char* file, int quiet, int do_comp);
		  
void diverge_seq(int n,int p,int features[n][p], int nos_fth[p], 
            int_vec* m_cases,int m_ptn[p+1],int no_cases_ign,
	    int_list* cases_g, int_list* ptns_g, int do_comp);
void find_sigmas_seq(int_list* ptns_g,int p,int sigmas_g[][p+1]);
/****************************************************************************/
void R_compress(int *R_sequence,int *R_order,
                    int *R_n,int *R_p, int features[*R_n][*R_p], 
                    int *R_nos_fth,int *R_no_cases_ign,
                    char** files, int* R_quiet, int* R_do_comp);

int_vec* take_cases(int_vec* m_cases,int n,int p,int features[n][p],
                    int id_fth,int value);

int is_unique(int n,int p,int features[n][p],
             int_vec* sub_cases,int id_fth);
	     
void load_compress
      (int_list* cases_g, int_list* ptns_g,
       int* P_sequence,int* P_order,int* P_n, int* P_p, 
       int** P_sigmas_g,char* file, int quiet);
       
void R_display_compress
     (char** files, int* no_show, int* gids, int info[6]);
