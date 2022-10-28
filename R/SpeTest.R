

########################################################################
###################### Specification test ##############################



SpeTest<-function(eq,type="icm",rejection="bootstrap",norma="no",
                  boot="wild",nboot=50,para=FALSE,ker="normal",knorm="sd",
                  cch="default",hv="default",nbeta="default",
                  direct="default",alphan="default"){
  
  if (inherits(eq,"lm") || inherits(eq,"nls")){
    
    
    ############################################################################
    ############# Central weighting Matrix W / "kernels" matrix ################
    
    ##### Standardized kernel smoothing function used to compute the central matrices
    
    k_sd <- function(x,ker='normal',knorm="sd")
    {
      
      # Densities such as the integral of the square of the density is one
      # if knorm is sq or such that the sd is one if knorm is sd.
      
      if (ker=='normal') # Normal pdf
      {
        s <- switch(knorm,
                    sd = 1,
                    sq = 1/(2*sqrt(pi)))
        out <- exp(- x^2 /(2*s^2)) / (sqrt(2*pi)*s)
      }
      
      if (ker=='triangle') # Triangle pdf
      {
        s <- switch(knorm,
                    sd = sqrt(6),
                    sq = (2/3))
        out <- (s-abs(x))
        out <- out * (out >0)/s^2
      }
      
      if (ker=='logistic') # Logistic pdf
      {
        s <- switch(knorm,
                    sd = 1/sqrt(2),
                    sq = 1/4)
        out <- exp(-abs(x)/s)/(2*s)
      }
      
      if (ker=="sinc") # Sine-cardinal function
      {
        out <- sqrt(3/2)*(sin(x*pi)/(x*pi))^2
        out <- apply(out, 1, function(x) replace(x,list=which(x=="NaN"),values=1))
      }
      
      return(out)
      
    }
    
    ### Compute column i of the matrix for a single dimension of x
    
    WB_sd <- function(x, i ,ker='normal',knorm="sd",remov=F, h){
      
      #   ICM and smoothing matrix
      #   The principle diagonal is replaced with zeros if remov = TRUE.
      
      n <- dim(x)[1]
      k <- dim(x)[2]
      
      W_pre <- x-x[i]
      
      # kernel smoothing
      
      W <-  k_sd(x=W_pre / h,ker=ker,knorm=knorm)/h
      
      # principle diagonal of the matrix replaced with zeros
      # "leave-one-out"
      
      if (remov==TRUE) 
        W[i] <-  0
      return(W)
      
    }
    
    
    ##### Bierens (1982) ICM and Zheng (1996) tests central matrices
    
    
    ### Compute the column i of the central matrix
    
    WB <- function(x,i,ker="normal",knorm="sd",remov=F, h){
      
      n <- dim(x)[1]
      k <- dim(x)[2]
      
      if (k > 1){
        W<-WB_sd(x=x[,1],i=i,ker=ker,knorm=knorm,remov=remov,h=h)
        for (j in 2:k){
          W <- W * WB_sd(x=x[,j],i=i,ker=ker,knorm=knorm,remov=remov,h=h)
        }
      } else {
        W <- WB_sd(x=x,i=i,ker=ker,knorm=knorm,remov=remov,h=h)
      }
      return(W)
    }
    
    
    
    
    
    ##### Escanciano (2006) test central matrix
    
    # Similar to the ICM weighting matrix except that
    # kernel functions are not used but different functions
    # which simplify integration over the unit hypersphere
    # See Escanciano (2006) for more details
    
    ### Compute element Aijr
    
    WE_el<-function(xi,xj,xr){
      
      k<-length(xr)
      
      if (k>1){
        if (prod(xi==xr) || prod(xj==xr)) {
          Aijr<-pi
        } else {
          Aij0<-abs(pi-acos(t(xi-xr)%*%(xj-xr)/(norm(xi-xr,type="O")*norm(xj-xr,type="O"))))
          Aijr<-Aij0*pi^(k/2-1)/gamma(k/2+1)
        }
      } else if (k==1){
        if (xi==xr || xj==xr) {
          Aijr<-pi
        } else {
          Aij0<-abs(pi-acos((xi-xr)*(xj-xr)/(norm(xi-xr,type="O")*norm(xj-xr,type="O"))))
          Aijr<-Aij0*pi^(k/2-1)/gamma(k/2+1)
        }
      }
      
      return(Aijr)
      
    }
    
    ### Compute element Aij
    
    WE_sel<-function(i,j,x){
      
      k<-dim(x)[2]
      
      if(k>1){
        Aij<-mean(apply(x,1, function(e) WE_el(xi=x[i,],xj=x[j,],xr=as.matrix(e))))
      } else if (k==1){
        Aij<-mean(sapply(x, function(e) WE_el(xi=x[i],xj=x[j],xr=as.matrix(e))))
      }
      
      return(Aij)
    }
    
    ### Compute column l of the central matrix
    
    WE_col<-function(x,l){
      n<-dim(x)[1]
      W<-as.matrix(sapply(seq(l,n), function(e) WE_sel(i=l,j=e,x=x)))
      return(W)
    }
    
    WE<-function(x){
      
      n<-dim(x)[1]
      
      W_pre<-sapply(seq(1,n), function(i) WE_col(x=x,l=i))
      W<-matrix(0,n,n)
      
      for (i in 1:(n)){
        W[i,i:n]<-W_pre[[i]]
      }
      
      W<-W+t(W)
      diag(W)<-diag(W)/2
      
      return(W)
    }
    
    
    
    ##### Lavergne and Patilea (2008) test weighting matrix
    
    # Similar to SICM except that the statistic
    # uses the random beta which maximizes it
    
    ### Draw a random beta in the hypersphere of norm 1
    
    r_beta<-function(k){
      
      b1<-as.matrix(rnorm(k))
      b2<-norm(b1,type="F")
      b3<-b1/b2
      return(b3)
      
    }
    
    ### Compute column i of the central matrix given a beta
    
    WP<-function(x,i,beta,ker='normal',knorm="sd",remov=T,h){
      
      k<-dim(x)[2]
      n<-dim(x)[1]
      
      if (k==1){
        x<-x/sd(x)
      } else {
        x<-apply(x[,1:k],2,function(e) e/sd(e))
      }
      
      if (k>1){
        x_beta<-x%*%as.matrix(beta)
      } else if (k==1){
        x_beta<-x*beta
      }
      
      W<-WB(x=x_beta,i=i,ker=ker,knorm=knorm,remov=remov,h=h)
      return(W)
      
    }
    
    ### Compute the criterion for maximization given a beta
    
    CP<-function(uhat,x,beta,beta0,alphan,ker='normal',knorm="sd",remov=T,h){
      
      n <- dim(x)[1]
      
      WU<-sapply(1:n, function(e) t(uhat)%*%WP(x=x,i=e,beta=beta,ker=ker,knorm=knorm,remov=remov,h=h))
      C<- abs(t(uhat)%*%WU/((n-1)*sqrt(h))-alphan*prod(beta!=beta0))
      
      return(C)
      
    }
    
    ##### Lavergne an Patilea (2012) SICM test central matrix
    
    # Here instead of integrating over the whole hypersphere
    # sums are used and a direction can be given
    
    ### Finds a random beta in the hypersphere norm 1
    ### which may follow the direction given initially
    
    r_beta_direct<-function(direct){
      
      Ok<-F
      
      while (Ok==F){
        
        b1<-as.matrix(rnorm(length(direct)))
        b2<-norm(b1,type="F")
        b3<-b1/b2
        
        s<-sum(sign(b3))
        sd<-sum(direct)
        sd0<-sum(direct==0)
        
        if ( s <= sd+sd0 & s >=sd-sd0){
          
          Ok=T
          bf<-sortr(b3,direct)
          
          
        } else if ( -s<= sd+sd0 & -s>=sd-sd0){
          
          Ok=T
          bf<-sortr(-b3,direct)
          
        }
      }
      
      return(as.numeric(bf))
      
    }
    
    ### Reorders the random beta so that it follows
    ### the initial direction
    
    sortr<-function(x,direct){
      
      d1<-which(direct==1)
      ld1<-length(d1)
      dm1<-which(direct==-1)
      ldm1<-length(dm1)
      
      x1<-which(x>0)
      lx1<-length(x1)
      xm1<-which(x<0)
      lxm1<-length(xm1)
      
      b0<-rep(0,length(direct))
      
      if (length(dm1)==0 & length(d1)==0){
        b0<-x
      } else if (length(dm1)==0){
        b0[d1]<-x[x1[1:ld1]]
        b0[b0==0]<-x[-x1[1:ld1]]
      } else if (length(d1)==0){
        b0[dm1]<-x[xm1[1:ldm1]]
        b0[b0==0]<-x[-xm1[1:ldm1]]
      } else {
        b0[dm1]<-x[xm1[1:ldm1]]
        b0[d1]<-x[x1[1:ld1]]
        b0[b0==0]<-x[-c(xm1[1:ldm1],x1[1:ld1])]
      }
      
      return(b0)
      
    }
    
    ### Compute column l of the central matrix
    
    WS<-function(x,hyper,l,ker="normal",knorm="sd",remov=F,h){
      
      n<-dim(x)[1]
      k<-dim(x)[2]
      
      if (k>1){
        W<-rowMeans(apply(hyper,2,function(e) WP(x,i=l,beta=e,ker=ker,knorm=knorm,remov=remov,h=h)))
      } else if (k==1){
        W<-WP(x,i=l,beta=hyper,ker=ker,knorm=knorm,remov=remov,h=h)
      }
      
      return(W)
      
    }
    
    
    #########################################################################
    ########### Non-parametric estim of variance of errors ################
    
    # Yin, Geng, Li and Wang (2010) estimator
    
    f_ratio<-function(a,b){
      S<-as.numeric(t(a)%*%b/sum(b))
      return(S)
    }
    
    f_u_bar<-function(uhat,x,i,ker="normal",knorm="sd",remov=T,h){
      u_bar<-f_ratio(uhat,WB(x,i,ker=ker,knorm=knorm,remov=remov,h=h))
      return(u_bar)
    }
    
    est_cova<-function(uhat,x,i,ker="normal",knorm="sd",remov=T,h){
      
      vec_cov<-sapply(i, function(e) f_ratio((uhat-f_u_bar(uhat=uhat,x=x,i=e,ker=ker,knorm=knorm,remov=remov,h=h))^2,
                                             WB(x=x,i=e,ker=ker,knorm=knorm,remov=remov,h=h)))
      return(vec_cov)
    }
    
    
    
    ################################################################
    #################### Individual Statistics #####################
    
    ### Bierens (1982)
    
    icm<-function(uhat,x,norma=F,cova="naive",ker,knorm,hv){
      
      n<-dim(x)[1]
      
      if (norma==F){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WB(x,e,ker=ker,knorm=knorm,remov=F,h=1))
        S<-as.numeric(t(uhat)%*%WU)/n
        
      } else  if (norma==T & cova=="naive"){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WB(x,e,ker=ker,knorm=knorm,remov=F,h=1))
        S_num<-as.numeric(t(uhat)%*%WU)
        WU2<-sapply(1:n, function(e) t(uhat^2)%*%(WB(x,e,ker=ker,knorm=knorm,remov=F,h=1)^2))
        S_den<-sqrt(as.numeric(t(uhat^2)%*%WU2)*2)
        S<-S_num/S_den
        
      } else if (norma==T & cova=="np"){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WB(x,e,ker=ker,knorm=knorm,remov=F,h=1))
        S_num<-as.numeric(t(uhat)%*%WU)
        uhatb<-est_cova(uhat,x,1:n,ker=ker,knorm=knorm,remov=F,h=hv)
        WU2<-sapply(1:n, function(e) t(uhatb)%*%(WB(x,e,ker=ker,knorm=knorm,remov=F,h=1)^2))
        S_den<-sqrt(as.numeric(t(uhatb)%*%WU2)*2)
        S<-S_num/S_den
        
      }
      
      return(S)
      
    }
    
    ### Zheng (1996)
    
    zheng<-function(uhat,x,norma=F,cova="naive",ker,knorm,cch,hv){
      
      n<-dim(x)[1]
      
      if (norma==F){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WB(x,e,ker=ker,remov=T,knorm=knorm,h=cch))
        S<-as.numeric(t(uhat)%*%WU)/n
        
      } else  if (norma==T & cova=="naive"){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WB(x,e,ker=ker,remov=T,knorm=knorm,h=cch))
        S_num<-as.numeric(t(uhat)%*%WU)
        WU2<-sapply(1:n, function(e) t(uhat^2)%*%(WB(x,e,ker=ker,remov=T,knorm=knorm,h=cch)^2))
        S_den<-sqrt(as.numeric(t(uhat^2)%*%WU2)*2)
        S<-S_num/S_den
        
      } else if (norma==T & cova=="np"){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WB(x,e,ker=ker,remov=T,knorm=knorm,h=cch))
        S_num<-as.numeric(t(uhat)%*%WU)
        uhatb<-est_cova(uhat,x,1:n,ker=ker,knorm=knorm,remov=T,h=hv)
        WU2<-sapply(1:n, function(e) t(uhatb)%*%(WB(x,e,ker=ker,remov=T,knorm=knorm,h=cch)^2))
        S_den<-sqrt(as.numeric(t(uhatb)%*%WU2)*2)
        S<-S_num/S_den
        
      }
      
      
      return(S)
      
    }
    
    ### Escanciano (2006)
    
    esca<-function(uhat,W_mat,norma=F,cova="naive",ker,knorm,hv){
      
      n<-dim(W_mat)[1]
      
      if (norma==F){
        
        WU<-W_mat%*%uhat
        S<-as.numeric(t(uhat)%*%WU/n)
        
      } else if (norma==T & cova=="naive"){
        
        WU<-W_mat%*%uhat
        S_num<-as.numeric(t(uhat)%*%WU)
        WU2<-(WU)^2
        S_den<-sqrt(as.numeric(t(uhat^2)%*%WU2)*2)
        S<-S_num/S_den
        
      } else if (norma==T & cova=="np"){
        
        WU<-W_mat%*%uhat
        S_num<-as.numeric(t(uhat)%*%WU)
        uhatb<-est_cova(uhat,x,1:n,ker=ker,knorm=knorm,remov=F,h=hv)
        WU2<-(W_mat%*%sqrt(uhatb))^2
        S_den<-sqrt(as.numeric(t(uhatb)%*%WU2)*2)
        S<-S_num/S_den
        
      }
      
      return(S)
      
    }
    
    ### Lavergne and Patilea (2008)
    
    pala<-function(uhat,x,hyper,norma=F,cova="naive",ker,knorm,cch,hv,
                   direct,alphan){
      
      n<-dim(x)[1]
      k<-dim(x)[2]
      
      if (k==1){
        betam<-hyper[which.max(sapply(hyper, function(e) CP(uhat=uhat,x=x,beta=e,beta0=direct,alphan=alphan,ker=ker,knorm=knorm,remov=T,h=cch)))]
      } else {
        betam<-hyper[,which.max(apply(hyper,2, function(e) CP(uhat=uhat,x=x,beta=e,beta0=direct,alphan=alphan,ker=ker,knorm=knorm,remov=T,h=cch)))]
      }
      
      if (norma==F){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WP(x=x,i=e,beta=betam,ker=ker,knorm=knorm,remov=T,h=cch))
        S<-as.numeric(t(uhat)%*%WU)/n
        
      } else if (norma==T & cova=="naive"){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WP(x=x,i=e,beta=betam,ker=ker,knorm=knorm,remov=T,h=cch))
        S_num<-as.numeric(t(uhat)%*%WU)
        WU2<-sapply(1:n, function(e) t(uhat^2)%*%(WP(x,i=e,beta=betam,ker=ker,knorm=knorm,remov=T,h=cch)^2))
        S_den<-sqrt(as.numeric(t(uhat^2)%*%WU2)*2)
        S<-S_num/S_den
        
      } else if (norma==T & cova=="np"){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WP(x=x,i=e,beta=betam,ker=ker,knorm=knorm,remov=T,h=cch))
        S_num<-as.numeric(t(uhat)%*%WU)
        uhatb<-est_cova(uhat,x,1:n,ker=ker,knorm=knorm,remov=T,h=hv)
        WU2<-sapply(1:n, function(e) t(uhatb)%*%(WP(x,i=e,beta=betam,ker=ker,knorm=knorm,remov=T,h=cch)^2))
        S_den<-sqrt(as.numeric(t(uhatb)%*%WU2)*2)
        S<-S_num/S_den
        
      }
      
      return(S)
      
    }
    
    ### Lavergne and Patilea (2012)
    
    sicm<-function(uhat,x,hyper,norma=F,cova="naive",ker,knorm,cch,hv){
      
      n<-dim(x)[1]
      
      if (norma==F){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WS(x=x,hyper=hyper,l=e,ker=ker,knorm=knorm,remov=T,h=cch))
        S<-as.numeric(t(uhat)%*%WU)/n
        
      } else if (norma==T & cova=="naive"){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WS(x=x,hyper=hyper,l=e,ker=ker,knorm=knorm,remov=T,h=cch))
        S_num<-as.numeric(t(uhat)%*%WU)
        WU2<-sapply(1:n, function(e) t(uhat^2)%*%(WS(x=x,hyper=hyper,l=e,ker=ker,knorm=knorm,remov=T,h=cch)^2))
        S_den<-sqrt(as.numeric(t(uhat^2)%*%WU2)*2)
        S<-S_num/S_den
        
      } else if (norma==T & cova=="np"){
        
        WU<-sapply(1:n, function(e) t(uhat)%*%WS(x=x,hyper=hyper,l=e,ker=ker,knorm=knorm,remov=T,h=cch))
        S_num<-as.numeric(t(uhat)%*%WU)
        uhatb<-est_cova(uhat,x,1:n,ker=ker,knorm=knorm,remov=T,h=hv)
        WU2<-sapply(1:n, function(e) t(uhatb)%*%(WS(x=x,hyper=hyper,l=e,ker=ker,knorm=knorm,remov=T,h=cch)^2))
        S_den<-sqrt(as.numeric(t(uhatb)%*%WU2)*2)
        S<-S_num/S_den
        
      }
      
      return(S)
      
    }
    
    
    
    #############################################################################
    ############################# Bootstrap methods #############################
    
    # Fixed design wild bootstrap weights
    
    w_boot<-function(n){
      
      pr <- 0.72360679774998   # (5+sqrt(5))/10
      c1 <-  - 0.6180339887499   # (1-sqrt(5))/2
      c2 <-  1.6180339887499    # (1+sqrt(5))/2
      vb <- runif(n)
      d1 <- (vb < pr)
      wei <- (c1*d1)+(c2*(1-d1))
      
      return(wei)
    }
    
    # Compute one bootstrap statistics
    
    s_boot <- function(uhat,w,x,x_sd,fit,eq,W_mat,hyper,type,norma,cova,ker,knorm,cch,hv,direct,alphan){
      
      # Create a vector of simulated outcomes
      
      ynew <- fit + w*uhat
      
      # Compute residuals for the new simulated outcomes
      
      if (inherits(eq,"lm")){
        
        uhatb <- lm(ynew~x)$residuals
        
      } else if (inherits(eq,"nls")){
        
        df<-x
        df[[as.character(eq$m$formula())[2]]]<-ynew
        uhatb <- as.numeric(resid(nls(formula = eq$m$formula(),data=df,
                                      start = eval(str2expression(as.character(eq$call)[4])),
                                      control = eq$control,
                                      trace = eval(str2expression(as.character(eq$call)[7])))))
      }
      
      return(switch(type,"icm"=icm(uhat=uhatb,x=x_sd,norma=norma,cova=cova,ker=ker,knorm=knorm,hv=hv),
                    "zheng"=zheng(uhat=uhatb,x=x_sd,norma=norma,cova=cova,ker=ker,knorm=knorm,cch=cch,hv=hv),
                    "esca"=esca(uhat=uhatb,W_mat=W_mat,norma=norma,cova=cova,ker=ker,knorm=knorm,hv=hv),
                    "pala"=pala(uhat=uhatb,x=x_sd,hyper=hyper,norma=norma,cova=cova,ker=ker,
                                knorm=knorm,cch=cch,hv=hv,direct=direct,alphan=alphan),
                    "sicm"=sicm(uhat=uhatb,x=x_sd,hyper=hyper,norma=norma,cova=cova,ker=ker,
                                knorm=knorm,cch=cch,hv=hv)))
      
    }
    
    # Obtain a vector of bootstrapped statistics
    
    v_boot <- function(uhat,x,x_sd,eq,W_mat,hyper,type,norma,cova,boot,nboot,ker,knorm,cch,hv,direct,alphan)
    {
      
      n<-dim(x)[1]
      k<-dim(x)[2]
      
      if (inherits(eq,"lm")){
        fit<-eq$fitted.values
      } else if (inherits(eq,"nls")){
        fit<-as.numeric(fitted(eq))
      }
      
      # Wild bootstrap or smooth conditional moment bootstrap
      # What is redrawn is either the residual itself or the estimated
      # standard deviation
      
      if (boot=="smooth"){
        remov<-switch(type,"icm"=F,"zheng"=T,"esca"=F,"pala"=T,"sicm"=T)
      }
      uhatb<-switch(boot,"wild"=uhat,"smooth"=sqrt(est_cova(uhat=uhat,x=x_sd,i=1:n,ker=ker,
                                                            knorm=knorm,
                                                            remov=remov,h=hv)))
      
      # Using a matrix to hold all the bootstrap weights
      # is too computationally intensive
      
      b <- sapply(1:nboot, function(e) s_boot(uhat=uhatb,w=w_boot(n),x=x,x_sd=x_sd,fit=fit,
                                              eq=eq,W_mat=W_mat,hyper=hyper,type=type,norma=norma,
                                              cova=cova,ker=ker,knorm=knorm,
                                              cch=cch,hv=hv,direct=direct,alphan=alphan))
      return(b)
      
    }
    
    
    
    #####################################################################
    ######################## Preliminary preparation ####################
    
    if (inherits(eq,"lm")){
      
      x<-as.matrix(eq$model[,-1])
      uhat<-as.numeric(eq$residuals)
      
    } else if (inherits(eq,"nls")){
      
      all_names<-names(eq$m$getEnv())
      para_names<-names(eq$m$getAllPars())
      x_names<-all_names[all_names!=as.character(eq$m$formula())[2] & all_names!=".swts" & apply(sapply(para_names, function(e) all_names!=e),1,prod)]
      x<-sapply(1:length(x_names), function(e) eq$m$getEnv()[[x_names[e]]])
      uhat<-as.numeric(resid(eq))
      
    }
    
    n<-dim(x)[1]
    k<-dim(x)[2]
    
    if (k==1){
      x_sd<-as.matrix(x/sd(x))
    } else if (k>1){
      x_sd<-as.matrix(apply(x,2,function(e) e/sd(e)))
    }
    
    ##### Directly compute weighting matrix if type = "esca"
    
    if (type=="esca"){
      W_mat<-WE(x=x_sd)
    }
    
    ##### Normalization
    
    if (norma=="no"){
      norma<-F
      cova<-"naive"
    } else if (norma=="naive"){
      norma<-T
      cova<-"naive"
    } else if (norma=="np"){
      norma<-T
      cova<-"np"
    }
    
    
    ##### Default bandwidth if type = "zheng" or type = "pala" or type ="sicm"
    
    if (type=="zheng" & cch=="default"){
      cch<-1.06*n^(-1/(4+k))
    }
    
    if ((type=="pala" || type=="sicm") & cch=="default"){
      cch<-1.06*n^(-1/5)
    }
    
    ##### Default bandwidth if the statistic is normalized with the nonparametric
    ##### covariance estimator
    
    if ((cova=="np" || boot=="smooth") & hv=="default"){
      hv=1.06*n^(-1/(4+k))
    }
    
    ##### Default number of betas in the unit hypersphere
    
    if (nbeta=="default"){
      nbeta=20*floor(sqrt(k))
    }
    
    ##### Default preferred alternative is the NLS estimator if type = "pala"
    
    if (type=="pala" & prod(direct=="default")){
      if (inherits(eq,"lm")){
        if (names(coefficients(eq))[1]!="(Intercept)"){
          direct<-as.numeric(coefficients(eq))
        } else {
          direct<-as.numeric(coefficients(eq))[-1]
        }
      } else if (inherits(eq,"nls")) {
        direct<-rep(0,k)
      }
    }
    
    ##### Default preferred direction in the unit hypersphere if type = "sicm"
    
    if (type=="sicm" & k>1){
      direct<-rep(0,k)
    }
    
    ##### Unit hypersphere if type = "pala" or if type = "sicm"
    
    if (type=="pala"){
      hyper<-sapply(1:nbeta, function(e) r_beta(k))
    }
    
    if (type=="sicm"){
      if (k==1){
        hyper<-1
      } else {
        hyper<-sapply(1:nbeta, function(e) r_beta_direct(direct))
      }
    }
    
    ##### Default preferred alternative if type = "pala"
    
    if (type=="pala" & alphan=="default"){
      alphan<-log(n)*n^(-3/2)
    }
    
    
    ##### If type = "icm" or if type = "esca" then the rejection rule is 
    ##### automatically based on the bootstrap
    
    if (type== "icm" || type=="esca"){
      
      rejection<-"bootstrap"
      
    }
    
    ##### If rejection is based on asymptotic theory, the statistic is automatically
    ##### normalized
    
    if (rejection=="asymptotics" & norma==F){
      
      norma<-T
      cova<-"naive"
      
    }
    
    SStat <- switch(type,"icm"=icm(uhat=uhat,x=x_sd,norma=norma,cova=cova,ker=ker,knorm=knorm,hv=hv),
                    "zheng"=zheng(uhat=uhat,x=x_sd,norma=norma,cova=cova,ker=ker,knorm=knorm,
                                  cch=cch,hv=hv),
                    "esca"=esca(uhat=uhat,W_mat=W_mat,norma=norma,cova=cov,ker=ker,knorm=knorm,hv=hv),
                    "pala"=pala(uhat=uhat,x=x_sd,hyper=hyper,norma=norma,cova=cova,ker=ker,knorm=knorm,cch=cch,
                                hv=hv,direct=direct,
                                alphan=alphan),
                    "sicm"=sicm(uhat=uhat,x=x_sd,hyper=hyper,norma=norma,cova=cova,ker=ker,knorm=knorm,cch=cch,
                                hv=hv))
    
    if (rejection=="bootstrap"){
      
      if (inherits(eq,"nls")){
        x<-data.frame(x)
        names(x)<-x_names
      }
      
      if (para==FALSE){
        
        Sboot<-v_boot(uhat=uhat,x=x,x_sd=x_sd,eq=eq,W_mat=W_mat,hyper=hyper,
                      norma=norma,cova=cova,type=type,boot=boot,nboot=nboot,
                      ker=ker,knorm=knorm,cch=cch,hv=hv,
                      direct=direct,alphan=alphan)
        
      } else if (para==TRUE){
        
        if (inherits(eq,"lm")){
          fit<-eq$fitted.values
        } else if (inherits(eq,"nls")){
          fit<-as.numeric(fitted(eq))
        }
        
        if (boot=="smooth"){
          remov2<-switch(type,"icm"=F,"zheng"=T,"esca"=F,"pala"=T,"sicm"=T)
        }
        uhatb<-switch(boot,"wild"=uhat,"smooth"=sqrt(est_cova(uhat=uhat,x=x_sd,i=1:n,ker=ker,
                                                              knorm=knorm,
                                                              remov=remov2,h=hv)))
        
        ncores<-parallel::detectCores()-1
        cluster<-parallel::makeCluster(mc <- getOption("cl.cores", ncores))
        doParallel::registerDoParallel(cluster)
        
        `%dopar%`<-foreach::`%dopar%`
        Sboot<-rep(0,nboot)
        Sboot<-foreach::foreach(is=1:nboot,.combine=c)%dopar%{s_boot(uhat=uhatb,w=w_boot(n),x=x,x_sd=x_sd,fit=fit,
                                                          eq=eq,W_mat=W_mat,hyper=hyper,type=type,norma=norma,
                                                          cova=cova,ker=ker,knorm=knorm,
                                                          cch=cch,hv=hv,direct=direct,alphan=alphan)}
        parallel::stopCluster(cluster)
      }
      
      PPval <- mean(Sboot>SStat)
      
    } else if (rejection=="asymptotics"){
      if (type=="zheng" || type=="pala" || type=="sicm"){
        PPval <- 1-pnorm(abs(SStat))
      } else if (type=="icm" || type=="esca"){
        PPval <-1-pchisq(SStat,df=1)
      }
      
    }
    
    final<-list(SStat,PPval,type,rejection,norma,cova,boot,nboot,ker,knorm,cch,hv,nbeta,direct,alphan)
    names(final)<-c("stat","pval","type","rejection","norma","cova","boot","nboot","ker","knorm","cch","hv","nbeta","direct","alphan")
    class(final)<-"STNP"
    
    
    return(final)
    
  } else {
    
    warning("'eq' is neither of class lm or of class nls")
    
  }
  
}




###########################################################################
########################### Utility functions #############################

##### Print function

print.STNP<-function(x, ...){
  
  name<-switch(x$type,"icm"="Bierens (1982) integrated conditional moment test","zheng"="Zheng (1996) test",
               "esca"="Escanciano (2006) test","pala"="Lavergne and Patilea (2008) test",
               "sicm"="Lavergne and Patilea (2012) smooth integrated conditional moment test")
  
  
  
  if (x$rejection=="asymptotics"){
    
    statf <- c("Normalized test statistic : ",round(x$stat,digits=5))
    pvalf<-c("Asymptotic p-value : ",round(x$pval,digits=5))
    
  } else if (x$rejection=="bootstrap"){
    if (x$norma==F){
      statf <- c("Test statistic : ",round(x$stat,digits=5))
    } else if (x$norma==T){
      statf <- c("Normalized test statistic : ",round(x$stat,digits=5))
    }  
    pvalf<-c("Bootstrap p-value : ",round(x$pval,digits=5))
    
  }
  
  cat("\n ",name,"\n\n ",statf,"\n ",pvalf,"\n ")
  
}


##### Summary function

summary.STNP<-function(object, ...){
  
  name<-switch(object$type,"icm"="Bierens (1982) integrated conditional moment test","zheng"="Zheng (1996) test",
               "esca"="Escanciano (2006) test","pala"="Lavergne and Patilea (2008) test",
               "sicm"="Lavergne and Patilea (2012) smooth integrated conditional moment test")
  
  cat("\n ", name)
  
  if (object$norma==F){
    
    statf <- c("Test statistic : ",round(object$stat,digits=5))
    
    cat("\n\n ", statf)
    
  } else if (object$norma==T){
    
    statf <- c("Normalized test statistic : ",round(object$stat,digits=5))
    cat("\n\n ", statf)
    
    if (object$norma==T){
      
      if (object$cova=="np"){
        
        hvp<-c("Nonparametric conditional variance estimator bandwidth : ",round(object$hv,digits=5))
        cat("\n ","Conditional variance estimator for normalization : Nonparametric","\n ",hvp)
        
      } else if (object$cova=="naive") {
        
        cat("\n ","Conditional variance estimator for normalization : Naive")
        
      }
      
    }
    
  }
  
  if (object$rejection=="asymptotics"){
    
    pvalf<-c("Asymptotic p-value : ",round(object$pval,digits=5))
    
  } else if (object$rejection=="bootstrap"){
    
    pvalf<-c("Bootstrap p-value : ",round(object$pval,digits=5))
    
  }
  
  cat("\n\n ",pvalf)
  
  if (object$rejection=="bootstrap"){
    
    bootp<-switch(object$boot,"wild"="Bootstrap type :  Wild bootstrap","smooth"="Bootstrap type :  Smooth conditional moments bootstrap")
    nbootp<-c("Number of bootstrap samples : ",object$nboot)
    cat("\n ",bootp,"\n ",nbootp)
    
    if (object$boot=="smooth"){
      
      hvpp<-c("Bootstrap conditional variance estimator bandwidth : ",round(object$hv,digits=5))
      cat("\n ",hvpp)
      
    }
    
  }
  
  kerp<-switch(object$ker,"normal"="Central matrix kernel function :  Normal p.d.f",
               "triangle"="Central matrix kernel function :  Triangular p.d.f",
               "logistic"="Central matrix kernel function :  Logistic p.d.f",
               "sinc"="Central matrix kernel function :  Sine Cardinal function")
  
  knormp<-switch(object$knorm,"sd"="Kernel standardization :  Standard deviation using the Kernel as a density = 1",
                 "sq"="Kernel standardization :  Integral of squared Kernel function = 1")
  
  cat("\n\n ", kerp, "\n ", knormp)
  
  if (object$type!= "esca" & object$type!= "icm"){
    
    cchp<-c("Kernel bandwidth : ",round(object$cch,digits=5))
    cat("\n ",cchp)
    
  }
  
  if (object$type=="pala" || object$type=="sicm"){
    
    nbetap<-c("Number of directions in the unit hypersphere : ", object$nbeta)
    cat("\n\n ",nbetap)
    
  }
  
  if (object$type=="pala"){
    
    directp<-c("Initial direction for beta = ",object$direct)
    cat("\n ",directp)
    
  }
  
  if (object$type=="pala"){
    
    alphanp<-c("Weight given to the favored alternative : ",object$alphan)
    cat("\n ",alphanp)
    
  }
  
  if (object$type=="sicm" & sum(object$direct!=0)>0){
    
    directp<-c("Initial direction for beta : ",round(object$direct,digits=5))
    cat("\n ",directp)
    
  }
  
  cat(" \n ")
  
}



