
double get_create(int i, double betas[], double sigma, int alpha, 
                  int* no_betas,char* betas_file)
{
   FILE *fp;
   if(betas[i] == 0 & i > 0)
   {
      GetRNGstate();	 
      if(alpha == 1)
         betas[i] = rcauchy(0,sigma);
      if(alpha == 2)
         betas[i] = rnorm(0,sigma);
      PutRNGstate();
      no_betas[0] ++;
      fp = fopen(betas_file,"a");
      fprintf(fp,"%d,%f\n",i,betas[i]);
      fclose(fp);       
   }
   return betas[i];
      
}

/**************************************************************************/

void add_lv(int i_beta, int i_fth, int o, int p, int power[],
            double lv[1], double betas[], int features[], 
	    double sigmas[], int alpha, int order, int* no_betas,
	    char* betas_file)
{
    if(i_fth == p - 1){
       lv[0] += get_create(i_beta,betas,sigmas[o],alpha,no_betas,
                           betas_file);
       
       if( o < order )
       lv[0] += get_create(i_beta+power[i_fth]*features[i_fth],
                           betas,sigmas[o+1],alpha,no_betas,
			   betas_file);
    }
    else{
       add_lv(i_beta,i_fth+1,o,p,power,lv,
              betas,features,sigmas,alpha,
              order,no_betas,betas_file);
       
       if( o < order)
       add_lv(i_beta+power[i_fth]*features[i_fth],
              i_fth + 1, o + 1, p, power,lv,betas,
	      features,sigmas,alpha,order,no_betas,
	      betas_file);    
    }
    
} 

/**************************************************************************/

void gen_bin_ho(int *n, int *p, int *order, int *alpha,
                int features[*n][*p], int nos_fth[*p],
                double sigmas[*order+1], int *no_betas, 
		int y[*n], double beta0[1],char** betas_file)
{
     int i,power[*p];
     double lv[1],prob,u;
     FILE* fp;
     
     power[0] = 1;
     for(i = 1; i < *p; i++)
        power[i] = power[i-1] * (nos_fth[i-1]+1);
     
     double betas[power[*p-1]*(nos_fth[*p-1]+1)];
     
     set_double_array(power[*p-1]*(nos_fth[*p-1]+1), betas, 0);
     
     betas[0] = *beta0;
     no_betas[0] = 1;
     fp = fopen(betas_file[0],"a");
     fprintf(fp,"%d,%f\n",0,betas[0]);
     fclose(fp);     
     
     for(i = 0; i < *n; i++){        
	lv[0] = 0;	
        add_lv(0,0,0,*p,power,lv,betas,&features[i][0],sigmas,*alpha,
               *order,no_betas,betas_file[0]);	
	prob = 1.0/(1.0 + exp(-lv[0]));	
	GetRNGstate();	
	y[i] = ((runif(0,1) < prob)?1:0) + 1;
	PutRNGstate();
     }
     
}

/**************************************************************************/

void gen_X( int* n, int* p, int* K,int X[*n][*p])
{
     int i,j;
     for(i = 0; i < *n; i++)
        for(j = 0; j < *p; j++)
	{
	   GetRNGstate();
	   X[i][j] = floor(rn_unif(1,*K+1));
	   PutRNGstate();
	}
     
}
