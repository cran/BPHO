
############################################################################

gen_hmm <- function(n,p,no_h,no_o, prob_h_stay, prob_o_stay){
   
   mod <- function(x,div) { x - floor(x/div)*div }
   probs_h <- matrix(runif(no_h^2),no_h,no_h)   
   diag(probs_h) <- prob_h_stay
   for(i in 1:no_h){
      probs_h[i,-i] <- probs_h[i,-i]/sum(probs_h[i,-i])*(1- prob_h_stay)
      
   }
   
   probs_h <- probs_h[sample(1:no_h),]
   
   d <- matrix(0,n,p+1)
   d[,1] <- sample(1:no_h,n,replace=TRUE)
   
   for(i in 1:n) {
      for(j in 2: (p+1)){
         d[i,j] <- sample(1:no_h,1,TRUE,probs_h[d[i,j-1],])
      }
   }
   
   probs_o <- matrix(runif(no_o^2),no_o,no_o)
   diag(probs_o) <- prob_o_stay
   
   for(i in 1:no_o)
   {
      probs_o[i,-i] <- (probs_o[i,-i]/sum(probs_o[i,-i])) * (1- prob_o_stay)
   }
   
   for(i in 1:n) {
       for(j in 1:(p+1)){
           d[i,j] <- sample(1:no_o, 1, TRUE, probs_o[mod(d[i,j],no_o)+1,] )
       }
   }

   list(X=d[,p:1],y=d[,p+1],probs_h = probs_h,probs_o=probs_o)
}




############################################################################


gen_bin_ho <- function(n,p,order,alpha,sigmas,nos_features,beta0)
{
   X <- matrix(0, n, p)
   
   betas_file="dfjsdlfjlsdfdkjfkldwo3xlaKAALJLSDSLJLSlwiwe1xfddfdldkfldfdfdfd"
     
   for(i in 1:p)
       X[,i] <- sample(1:nos_features[i],nrow(X),replace = TRUE)
   if(file.exists(betas_file)) file.remove(betas_file)
   out = .C("gen_bin_ho",
      as.integer(n),as.integer(p),as.integer(order),
      as.integer(alpha),as.integer(t(X)),as.integer(nos_features),
      c(0,sigmas), no_betas= as.integer(0),y=integer(n),
      beta0,betas_file,PACKAGE="BPHO")[8:9]
   betas <- read.table(betas_file,header=FALSE,sep=",")
   file.remove(betas_file)
   c(list(X=X,betas=betas),out)
}

############################################################################
text_to_3number <- function(p,file){
   vtext <- scan(file,what="character",quiet=TRUE)
   ctext <- c()
   for(i in 1:length(vtext)){
     ctext <- c(ctext,strsplit(vtext[i],split="")[[1]]," ")
   }

   
   vowl = c(strsplit("aeiou",split="")[[1]],
            strsplit("AEIOU",split="")[[1]])
   cons = c(strsplit("bcdfghjklmnpqrstvwxyz",split="")[[1]],
            strsplit("BCDFGHJKLMNPQRSTVWXYZ",split="")[[1]])
   ntext <- rep(1,length(ctext))	    
   for(i in 1:length(ctext))
   {  if(sum(vowl==ctext[i]))
         ntext[i] <- 2
      else 
         if(sum(cons==ctext[i]))
	    ntext[i] <- 3	 
   }
   
   ntext2 <- numeric(0)
   i <- 1
   while(i < length(ntext) + 1)
   {  
       ntext2 <- c(ntext2,ntext[i])
       i <- i + 1
       if(ntext[i-1] == 1){
         while(i < length(ntext)+1  )
	   if( ntext[i] == 1 )
             i <- i + 1
	   else break
       } 
   }
   
   ntext3 <- matrix(0,length(ntext2)-p,p+1)
   for(i in 1:nrow(ntext3))
      ntext3[i,] <- ntext2[i:(i+p)]
  
   list(X=ntext3[,p:1],y=ntext3[,p+1])
}

text_to_number <- function(p,file){
   wtext <- scan(file,what="character",quiet=TRUE,blank.lines.skip=TRUE)
   ctext <- c()
   for(i in 1:length(wtext)){
     ctext <- c(ctext,strsplit(wtext[i],split="")[[1]]," ")
   }
   lowercase <- strsplit("abcdefghijklmnopqrstuvwxyz",split="")[[1]]
   uppercase <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ",split="")[[1]]
   convert_one_char <- function(char)
   {  number <- c(27,(1:26)[lowercase==char],(1:26)[uppercase==char])
      number[length(number)]
   }
   ntext <- sapply(ctext,convert_one_char)
   X <- matrix(0,length(ntext)-p,p)
   y <- rep(0,length(ntext) - p)
   for(i in seq(1,length(y)) ) {
      y[i] <- ntext[i + p]
      X[i,] <- ntext[seq(i,i+p-1)]
   }
   list(X=X,y=y)
}

############################################################################

gen_X <- function(n,p,K)
{
   X <- integer(n*p)
   matrix(.C("gen_X",as.integer(n),
             as.integer(p),as.integer(K),X=t(X),PACKAGE="BPHO")$X,
             n,p,byrow=TRUE)
}

