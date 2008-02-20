void R_read_mc(char **group, int *ix, char **mc_file, int *iter_b,
                int *forward, int *n, double out[]){
    read_mc(group[0],ix[0],mc_file[0],iter_b[0],
              forward[0],n[0],out);

}
/**************************************************************************/

void read_mc
   (char* group, int ix, char *mc_file, int iter_b, int forward, int n,
    double out[]){
    FILE *fp_mc;
    int no_iters, no_cls, order, no_g, alpha, known=0;
    long head_in_iter,size_each_iter;
    fp_mc = fopen_chk(mc_file,"rb");
    fseek(fp_mc, - sizeof(int), SEEK_END);
    fread_chk(&no_iters,sizeof(int),1,fp_mc,mc_file);

    while(iter_b > no_iters - 1 ){
        printf("Your first iteration index is out of range\n");
	printf("Please input a new one:\n");
	scanf("%d",&iter_b);
    }

    rewind(fp_mc);
    fread_chk(&no_cls,sizeof(int),1,fp_mc,mc_file);
    fread_chk(&no_g,sizeof(int),1,fp_mc,mc_file);
    fread_chk(&order,sizeof(int),1,fp_mc,mc_file);
    fread_chk(&alpha,sizeof(int),1,fp_mc,mc_file);

    size_each_iter = sizeof(int)+sizeof(double)*(no_g*(no_cls+1)+order+7);

    while((known++) == 0){
       head_in_iter = sizeof(int);
       if(strcmp(group,"betas") == 0){
           while(!(ix < no_g*no_cls)){
	       printf("Index should be from %d to %d\n",
	              0,no_g*no_cls-1);
	       printf("Please input new index:\n");
	       scanf("%d",&ix);
	   }
	   head_in_iter += ix * sizeof(double);
       }
       else{
           head_in_iter += no_g * no_cls * sizeof(double);
	   if(strcmp(group,"lsigmas") == 0){
	      while(!(ix < order+1)){
	         printf("Index should be from %d to %d\n",
	              0,order);
	         printf("Please input new index:\n");
	         scanf("%d", &ix);
	       }
	      head_in_iter +=  ix * sizeof(double);
	   }
	   else{
	       head_in_iter += sizeof(double) * (order+1);
	       if(strcmp(group,"lprobs")==0){
	          while(!(ix < 4)){
	             printf("Index should be from 0 to 3\n");
	             printf("Please input new index:\n");
	             scanf("%d",&ix);
	          }
		  head_in_iter += ix * sizeof(double);
	       }
	       else{
	          head_in_iter += sizeof(double) * 4;
		  if(strcmp(group,"evals")==0){
		     while(!(ix < 2)){
	             printf("Index should be from 0 to 1\n");
	             printf("Please input new index:\n");
	             scanf("%d",&ix);
	             }
		     head_in_iter += ix * sizeof(double);
		  }
		     else{
		         printf("Do not know what you want to read\n");
			 printf("Group name should be one from below:\n");
			 printf("`betas',`lsigmas',`lprobs',`evals'\n");
			 printf("Please select one from above:\n");
			 scanf("%s",group);
		         known = 0;
		     }
	       }

	  }
      }
   }
   fseek(fp_mc, iter_b * size_each_iter + head_in_iter, SEEK_CUR);
   int i = 0;
   while(iter_b+i*forward < no_iters & i < n){
       fread_chk(&out[i],sizeof(double),1,fp_mc,mc_file);
       fseek(fp_mc, forward*size_each_iter - sizeof(double),SEEK_CUR);
       i++;
   }
   fclose(fp_mc);
}
/**************************************************************************/

void display_mc(char** mc_files,int out[5])
{
   FILE* fp_mc;
   fp_mc = fopen_chk(mc_files[0],"rb");
   fseek(fp_mc, - sizeof(int), SEEK_END);
   /*reading no_iters*/
   fread_chk(&out[0],sizeof(int),1,fp_mc,mc_files[0]);
   rewind(fp_mc);
   /*reading no_cls*/
   fread_chk(&out[1],sizeof(int),1,fp_mc,mc_files[0]);
   /*reading no_g*/
   fread_chk(&out[2],sizeof(int),1,fp_mc,mc_files[0]);
   /*reading order*/
   fread_chk(&out[3],sizeof(int),1,fp_mc,mc_files[0]);
   /*reading alpha*/
   fread_chk(&out[4],sizeof(int),1,fp_mc,mc_files[0]);
   fclose(fp_mc);
}



