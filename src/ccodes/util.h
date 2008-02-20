double find_double_max(int la,double a[]);

double sums(int len,double A[len]);

double sum_double(int len,double A[len]);

double log_sum_exp(int la, double a[]);

double log_add_exp(double loga, double logb);

double log_diff_exp(double loga, double logb);

void set_double_array(int la, double a[], double value);

int sum_int(int len,int A[len]);

void fwrite_chk
(const void *buffer, size_t num_bytes, size_t count, 
 FILE *fp, char file[]);
 
void fread_chk
(void *buffer, size_t num_bytes, size_t count, 
 FILE *fp, char file[]);
 
double rn_unif(double a, double b);

void copy_double_array(int la,double *a1, double *a2);

FILE *fopen_chk(char file[],const char *mode);

void subtract_double_array(int la, double a[], double value);
void norm_log_array(int la, double a[]);

double log_minus_add_exp(double logs,double loga, double logb);
void print_int_array(int la, int a[],int items_line);
