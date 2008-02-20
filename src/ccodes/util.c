double sq(double x)
{  return R_pow_di(x,2);
}
/***********************************************************************/
double cb(double x)
{  return R_pow_di(x,3);
}
/***********************************************************************/
double quan(double x)
{  return R_pow_di(x,4);
}
/***********************************************************************/
double find_double_max(int la,double a[]){
    int i;
    int m = a[0];
    for(i = 1; i < la; i++){
       if(a[i] > m) m = a[i];
    }
    return(m);
}
/***********************************************************************/
double sum_double(int len,double A[len])
{  int i; double out=0;
   for(i = 0; i < len; i++) out += A[i];
   return(out);  
} 
/***********************************************************************/

int sum_int(int len,int A[len])
{  int i, out=0;
   for(i = 0; i < len; i++) out += A[i];
   return(out);  
} 
/***********************************************************************/
double avg_double(int len,double A[len])
{  return(sum_double(len,A)/len);  
} 

/***********************************************************************/
double log_sum_exp(int la, double a[])
{  int i;
   double m,out=0;
   m = find_double_max(la,a);
   for(i = 0; i < la; i++){
      out += exp(a[i] - m);
   }
   out = log(out) + m;
   return(out);
}

/***********************************************************************/
void subtract_double_array(int la, double a[], double value)
{  
    int i;
    for(i = 0; i < la; i++)
       a[i] -= value;
}
/***********************************************************************/

void norm_log_array(int la, double a[])
{   
   subtract_double_array(la, a, log_sum_exp(la,a));
}
/***********************************************************************/

/*loga must greater than logb*/
double log_diff_exp(double loga, double logb){
   double m;
   m = loga;
   return(log(exp(loga-m) - exp(logb-m)) + m);
}
/***********************************************************************/

double log_add_exp(double loga, double logb){
   double m;
   if(loga > logb) m = loga;
   else m = logb;
   return(log(exp(loga-m) + exp(logb-m)) + m);
}
/***********************************************************************/
double log_minus_add_exp(double logs,double loga, double logb){
   Rprintf("logs %f loga %f logb %f\n",logs,loga,logb);
   if(loga > logs) {
     exit(1);
   }
   return(logs+log(1.0-exp(loga-logs)+exp(logb-logs)));
}

/***********************************************************************/

void set_double_array(int la, double a[], double value){
   int i;
   for(i = 0; i < la; i++)
      a[i] = value;
    
}
/***********************************************************************/

/* a wrapper function that calls fwrite and does checking */
void fwrite_chk
(const void *buffer, size_t num_bytes, size_t count, 
 FILE *fp, char file[]){
   int MAX_TRY = 50;
   
   size_t no = 0; 
   int time = 0;
   
   do{
      if(no > 0) fseek(fp,-no*num_bytes,SEEK_CUR);
      no = fwrite(buffer,num_bytes,count,fp);
      time++;      
   }
   while( no!= count & time < MAX_TRY);
    
   if(time == MAX_TRY){
      Rprintf("I give up writing into file %s\n", file);
      exit(1);
   }   
}
/***********************************************************************/
void fread_chk
(void *buffer, size_t num_bytes, size_t count, 
 FILE *fp, char file[]){
   int MAX_TRY = 50;
   
   size_t no = 0; 
   int time = 0;
   
   do{
      if(no > 0) fseek(fp,-no*num_bytes,SEEK_CUR);
      no = fread(buffer,num_bytes,count,fp);
      time++;      
   }
   while( no!= count & time < MAX_TRY);
    
   if(time == MAX_TRY){
      Rprintf("I give up reading from file %s\n", file);
      exit(1);
   }   
}
/***********************************************************************/

FILE *fopen_chk(char file[],const char *mode){
    FILE *fp;
    int time=0,MAX_TRY = 50;
     do{
       fp = fopen(file,mode);
       time++;
     } while(fp == NULL & time < MAX_TRY);
    if(time == MAX_TRY){
       Rprintf("I give up opening file %s\n", file);
    }
    return(fp);
}
/***********************************************************************/

double rn_unif(double a, double b){
   return(a + (rand()+1.0)/(RAND_MAX + 2.0) * (b-a) );
}
/***********************************************************************/

void copy_double_array(int la,double *a1, double *a2){
    int i;
    for(i = 0; i < la; i++)
       a1[i] = a2[i];
}

void R_ldgamma(double q[],double alpha[], double scale[], double p[])
{
    p[0] = dgamma(q[0],alpha[0],scale[0],1);
}

/***********************************************************************/
void print_int_array(int la, int a[],int items_line)
{  int i;
   for(i = 0; i < la; i++){
      printf("%d ",a[i]);
      if(i%items_line == items_line - 1) printf("\n");
   }
   printf("\n");
   
}
