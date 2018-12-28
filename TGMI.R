
################### DO NOT MODIFY THE CODE IN THIS FILE, PLEASE CHANGE PARAMETERS IN "run_TGMI.R" script. ##############

########################################################################################################################

registerDoMC(cores=ncore)
################### DO NOT MODIFY THE CODE IN THIS FILE, PLEASE CHANGE PARAMETERS IN "run_TGMI.R" script. ##############

########################################################################################################################

create_empty_table <- function(num_rows, num_cols) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  return(frame)
}

expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}


mi3<-function(y1,y2,x){
  
  y1<-discretize(y1)
  y2<-discretize(y2)
  x<-discretize(x)
  
  enty1 <- entropy(y1)
  
  enty2 <- entropy(y2)
  
  entx <- entropy(x)
  
  
  enty1x <- condentropy(y1,x, method="emp")
  
  enty2x <- condentropy(y2,x, method="emp")
  
  enty1_y2x <- condentropy(y1,data.frame(x,y2), method="emp")
  
  enty2_y1x <- condentropy(y2,data.frame(x,y1), method="emp")
  
  entx_y1y2 <- condentropy(x,data.frame(y2,y1), method="emp")
  
  #calculates MI y1;y2
  Iy1y2 <- condinformation(y1, y2, method="emp")
  
  #calculates MI y1;y2|x
  Iy1y2_x <- condinformation(y1, y2,x, method="emp")
  
  Iy1x_y2<- condinformation(y1,x,y2, method="emp")
  
  Iy2x_y1 <- condinformation(y2,x,y1, method="emp")
  
  Iy1y2x <- Iy1y2 - Iy1y2_x
  I3<-abs(Iy1y2x/(enty1_y2x+enty2_y1x+entx_y1y2))
  if(I3>0){
    return(I3)
  }else{
    return(0)
  }
}

generate_permute<-function(y1,y2,x){
  vec<-c()
  for(i in 1:1000){
    vec<-c(vec,mi3(y1,y2,sample(x)))
  }
  return(vec)
}


lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

selectTFs <- function(DF,alpha) {
  selectTF<-c()
  #correcting pvalues
  #https://www.statisticssolutions.com/bonferroni-correction/
  altered_p<-alpha/dim(DF)[1]
  alcritical = 1 - (1 - altered_p)^dim(DF)[1]
  
  for(ind in seq(1,dim(DF)[2],2)){
    tempDF<-DF[,c(ind,ind+1)]
    tempDF<-tempDF[tempDF[,2] < alcritical,]
    selectTF<-c(selectTF,as.character(tempDF[1:50,1]))
  }
  selectTF<-unique(selectTF)
  return(selectTF)
}




func<-function(pw1,pw2,tf){
  y1<-pw1
  y2<-pw2
  x <-tf
  y1<-discretize(y1)
  y2<-discretize(y2)
  x<-discretize(x)

  enty1 <- entropy(y1)
  #if(enty1==0) return(NULL)
  enty2 <- entropy(y2)
  #if(enty2==0) return(NULL)
  entx <- entropy(x)
  #if(entx==0) return(NULL)

  enty1x <- condentropy(y1,x, method="emp")
  #if(enty1x==0) return(NULL)

  enty2x <- condentropy(y2,x, method="emp")
  #if(enty2x==0) return(NULL)

  enty1_y2x <- condentropy(y1,data.frame(x,y2), method="emp")
  #if(enty1_y2x==0) return(NULL)

  enty2_y1x <- condentropy(y2,data.frame(x,y1), method="emp")
  #if(enty2_y1x==0) return(NULL)

  entx_y1y2 <- condentropy(x,data.frame(y2,y1), method="emp")
  #if(entx_y1y2==0) return(NULL)

  #calculates MI y1;y2
  Iy1y2 <- condinformation(y1, y2, method="emp")
  #if(Iy1y2==0) return(NULL)

  #calculates pvalue for MI y1;y2

  #calculates MI y1;y2|x
  Iy1y2_x <- condinformation(y1, y2,x, method="emp")
  #if(Iy1y2_x==0) return(NULL)

  Iy1x_y2<- condinformation(y1,x,y2, method="emp")
  #if(Iy1x_y2==0) return(NULL)

  Iy2x_y1 <- condinformation(y2,x,y1, method="emp")
  #if(Iy2x_y1==0) return(NULL)

  Iy1y2x <- Iy1y2 - Iy1y2_x
  #if(Iy1y2x <= 0) return(NULL)
  #I3 <-0.105713486623665
  I3<-Iy1y2x/(enty1_y2x+enty2_y1x+entx_y1y2)
  z_score<-(I3-mean(randomized_3GI_2PW_1TF))/sd(randomized_3GI_2PW_1TF)
  I3_pval <- pnorm(z_score,lower.tail=FALSE)
  #I3_pval <- generate_permut_pval(y1,y2,x,I3)

  if(enty1==0||enty2==0||entx==0||enty1x==0||enty2x==0||enty1_y2x==0||
     enty2_y1x==0||entx_y1y2==0||Iy1y2==0||Iy1y2_x==0||Iy1x_y2==0||Iy2x_y1==0||
     Iy1y2x <= 0){
    return(NULL)
  }else{
  return(c(enty1,
           enty2,
           entx,
           enty1_y2x,
           enty2_y1x,
           entx_y1y2,
           Iy1x_y2,
           Iy2x_y1,
           Iy1y2_x,
           Iy1y2x,
           I3,
           I3_pval
         ))
  }
}

parallel_approach1<-function(Y,X,gn,tfn)
{

  n<-dim(Y)[2]
  m<-dim(X)[2]
  nsample<-dim(X)[1]
    result=
    foreach(i=1:(n-1),.combine="rbind") %dopar%
    {
      foreach(j =(i+1):n,.combine="rbind") %dopar%
      {
        #match("COMT1",colnames(pw.exp))
        #k=10
        #i=25
        #j=27
        tempout=data.frame(Y1=character(0),Y2=character(0),X=character(0),
                           H_Y1=numeric(0),H_Y2=numeric(0),H_X=numeric(0),
                           S1=numeric(0),S2=numeric(0),S3=numeric(0),
                           S4=numeric(0),S5=numeric(0),S6=numeric(0),S7=numeric(0),
                           S7_div_123=character(0),S7_div_123_pval=character(0))

        pw1<-as.numeric(Y[,i])
        pw2<-as.numeric(Y[,j])
        for(k in 1:m)
        {
          tf<-as.numeric(X[,k])
          infotemp<-func(pw1,pw2,tf)
          if(length(infotemp)==0){
            next
          }
          out<-cbind(Y1=gn[i],Y2=gn[j],X=tfn[k],H_Y1=infotemp[1],H_Y2=infotemp[2],H_X=infotemp[3],S1=infotemp[4],S2=infotemp[5],
                     S3=infotemp[6],S4=infotemp[7],S5=infotemp[8],S6=infotemp[9],S7= infotemp[10],S7_div_123=infotemp[11],
                     S7_div_123_pval=infotemp[12])
          tempout=rbind(tempout,out)
        }
        tempout
      }
    }

return(result)
}

generate_permut_pval<-function(x,y,z,I3){
  vec<-c()
  for(i in 1:1000){
    vec<-c(vec,mi3(sample(x),sample(y),sample(z)))
  }

  p.value<-sum(vec>I3)/length(vec)

  z_score<-(I3-mean(randomized_3GI_2PW_1TF))/sd(randomized_3GI_2PW_1TF)
  #pnorm(z_score,lower.tail=FALSE)
  I3_pval <- pnorm(z_score,lower.tail=FALSE)

  return(I3_pval)
}

randomized_3GI_2TF_1PW<-c()
i2 <- floor(runif(1,1,dim(tf.data)[2]))
j2 <- floor(runif(1,1,dim(tf.data)[2]))
k2 <- floor(runif(1,1,dim(pathway.data)[2]))
for(i in 1:1000){
  randomized_3GI_2TF_1PW<-c(randomized_3GI_2TF_1PW,mi3(sample(tf.data[,i2]),sample(tf.data[,j2]),sample(pathway.data[,k2])))
}

randomized_3GI_2PW_1TF<-c()
i1 <- floor(runif(1,1,dim(pathway.data)[2]))
j1 <- floor(runif(1,1,dim(pathway.data)[2]))
k1 <- floor(runif(1,1,dim(tf.data)[2]))
for(i in 1:1000){
  randomized_3GI_2PW_1TF<-c(randomized_3GI_2PW_1TF,mi3(sample(pathway.data[,i1]),sample(pathway.data[,j1]),sample(tf.data[,k1])))
}



