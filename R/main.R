#'Standard Algorithm for Bayesian Network Discovery
#'@param pvalue a vector of p-values obtained from large scale statistical hypothesis testing
#'@param net a n by n network configuration, n is the length of pvalue
#'@param iter number of iterations. The default is 5000
#'@param nburns number of burn-in. The default is 2000
#'@param piall a vector of selections of pi0. The default vector is 0.75, 0.8, 0.85, 0.9. The selections of pi0 should be placed in sequence, from smaller to larger.
#'@param rhoall a vector of selections of rho0 and rho1. The default vector is 0.5, 1, 5, 10, 15. The selections of rho0 and rho1 should be placed in sequence, from smaller to larger.
#'@details This generic function fits a Bayesian Nonparametric Mixture Model for gene selection incorporating network information (Zhao et al., 2014):
#' \itemize{
#' \item  r_i| g_i, \strong{theta} ~ N(mu_{g_i}, sigma_{g_i}),
#' \item	g_i | z_i=k, \strong{q}_k ~ Discrete(\strong{a}_k, \strong{q}_k),
#' \item	\strong{theta} ~ G_{0k}, for g in \strong{a}_k,
#' \item	\strong{q}_k ~ Dirichlet(tau_k \strong{1}_{L_k}/L_k),
#' \item	\strong{theta}={\strong{theta}_g}_{g in \strong{a}_0 and \strong{a}_1} 
#' \item	\strong{theta}_g=(mu_g, sigma_g)
#' }
#' where we define 
#' \describe{
#' \item{Index}{\strong{a}_0=(-L_0+1,-L_0+2,...,0) , \strong{a}_1=(1,2,...,L_1) and the correspondent probability q_0=(q_{-L_0+1}, q_{-L_0+2}, ...,q_0), q_1=(q_1, q_2, ..., q_{L_1}), according to the defination of Discrete(\strong{a}_k, \strong{b}_k), for example, Pr(g_i={L_0+2})=q_{-L_0+2}. }
#' \item{Assumption}{We have an assumption that "selected" gene or image pixel should have larger statiscs comparing to "unselected" ones without the loss of generality. In this regard, we set the restriction mu_g<mu_{g+1} for g=-L_0+1, -L_0+2,...,L_1.}
#' }
#' For this function, The NET-DPM-1, considered as standard function is applied , and more details about the algorithm can be referred from Appnendix B.1 of Zhao et al., 2014
#'@return The trace of gi showing the evolution of the Monte Carlo Markov Chain
#'@examples #' ##Creating the network of 10X10 image
#' library(igraph)
#' library(BayesNetDiscovery)
#' g <- graph.lattice(length=10,dim=2)
#' net=as(get.adjacency(g,attr=NULL),"matrix")##this is the input of argument \code{net}
#' ##Assign the signal elements with signal intenstion as normal distribution N(1,0.2). While noise is set as N(0,0.2) 
#' newz=rep(0,100)
#' for (i in 3:7)
#' {
#'  newz[(i*10+3):(i*10+7)]=1
#' }
#' testcov<-0
#' for(i in 1:100){
#'  if(newz[i]==0){
#'    testcov[i]<-rnorm(1,mean=0,sd=0.2)
#'  
#'  }else{
#'   testcov[i]<-rnorm(1,mean=1,sd=0.2)
#'    
#'  }
#' }
#' ##The profile of the impage
#' image(matrix(testcov,10,10),col=gray(seq(0,1,length=255)))
#' ##Transform the signals into pvalue form and begin identification
#' pvalue=pnorm(-testcov)
#' total=Networks.STD(pvalue,net,iter=5000,piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(0.5,1,5,10,15))
#'@references Zhao, Y.*, Kang, J., Yu, T. A Bayesian nonparametric mixture model for gene and gene-sub network selection
#' Annals of Applied Statistics, In press: 2014. 
#'@export
Networks.STD=function(pvalue,net,iter=5000, nburns=2000, piall=c(0.75, 0.8, 0.85, 0.9),rhoall=c(0.5, 1, 5, 10, 15),trace=TRUE,show.steps=1){
  ####Preparing Stage:Defining Parameters, 
  rstat=Transfer(pvalue)###Step1 : Transferring 
  wholeindex=Kmeans(rstat)###Step2: Using Kmeans for giving the initial values
  initial_mu_var=Initial_mu_var(rstat,wholeindex) ###Step3: Genering initial mu and var based on intialindex
  pirhopair=Pi_rho_gen(piall,rhoall)####Step4.1 Generating all the possible selections of pi and rho
  znew=Generating_Zi(net,pirhopair,wholeindex)####Step4.2 Generating Zi by different setting of pi and rho.
  choice=Pi_rho_selecting(znew,rstat,piall,rhoall)####Step 4.3 Using likelyhood to select best set of pi and rho
  parameter=Parameter_Define(initial_mu_var,pirhopair,choice)####Defining parameters
  #parameter$pi0=0.85
  #parameter$pi1=0.15
  #parameter$rho0=0.5
  #parameter$rho1=1
  #return(parameter)
  ####Iteration
  model=Model(wholeindex,net,rstat,parameter,initial_mu_var)#####Defining the values which would be used in the later iterations
  
  ####iteration begin~
  sample.save = matrix(NA,ncol=length(pvalue),nrow=iter-nburns)
  for (jj in 1:iter)
  {
    if(trace==TRUE){
    if(jj%%show.steps==0){
    cat("Iter:",jj,"\n")
    flush.console()}}
   for (num1 in 1:length(rstat)){
     model=Step1_1_Update_gi_zi(model,num1)
     model=Step1_2_Check(model=model,num1)
     model=Step_2_Update_mu_var(model,num1)
     model=Step_3_switch_the_label(model,num1)
   }
   if(jj>nburns)
     sample.save[jj-nburns,]=model$wholeindex
  }
  
  if (trace==TRUE){
       show=c()
      
                        for (kk in 1:nrow(sample.save))
                        {
                         if(kk%%show.steps==0){show=rbind(show,sample.save[kk,])}
                        }
          print(show)
  }
  invisible(sample.save)
  
}



#'Fast Algorithm for Bayesian Network Discovery
#'@param pvalue a vector of p-values obtained from large scale statistical hypothesis testing
#'@param net a n by n network configuration, n is the length of pvalue
#'@param iter number of iterations. The default is 5000
#'@param nburn number of burn-in. The default is 2000
#'@param initials the way you choose to get the inital posteritor samples of parameters produced by DPM fitting. Generating initals by function Mclust() when you set "mclust" or by DPdensity() when you set "DPdensity". It is recommended to choose "Mclust" when the dimension of data is large, and to choose "DPdensity" when the dimenstion is small.
#'@param v number of iterations set for DPM fitting by "DPdensity". v is only valid when you choose initals as "DPdensity"
#'@param piall a vector of selections of pi0. The default vector is 0.75, 0.8, 0.85, 0.9. The selections of pi0 should be placed in sequence, from smaller to larger.
#'@param rhoall a vector of selections of rho0 and rho1. The default vector is 1, 2, 5, 10, 15. The selections of rho0 and rho1 should be placed in sequence, from smaller to larger.
#'@details This generic function fits a Bayesian Nonparametric Mixture Model for gene selection incorporating network information (Zhao et al., 2014):
#' \itemize{
#' \item  r_i| g_i, \strong{theta} ~ N(mu_{g_i}, sigma_{g_i}),
#' \item  g_i | z_i=k, \strong{q}_k ~ Discrete(\strong{a}_k, \strong{q}_k),
#' \item	\strong{theta} ~ G_{0k}, for g in \strong{a}_k,
#' \item	\strong{q}_k ~ Dirichlet(tau_k \strong{1}_{L_k}/L_k),
#' \item	\strong{theta}={\strong{theta}_g}_{g in \strong{a}_0 and \strong{a}_1} 
#' \item	\strong{theta}_g=(mu_g, sigma_g)
#' }
#' where we define 
#' \describe{
#' \item{Index}{\strong{a}_0=(-L_0+1,-L_0+2,...,0) , \strong{a}_1=(1,2,...,L_1) and the correspondent probability q_0=(q_{-L_0+1}, q_{-L_0+2}, ...,q_0), q_1=(q_1, q_2, ..., q_{L_1}), according to the defination of Discrete(\strong{a}_k, \strong{b}_k), for example, Pr(g_i={L_0+2})=q_{-L_0+2}. }
#' \item{Assumption}{We have an assumption that "selected" gene or image pixel should have larger statiscs comparing to "unselected" ones without the loss of generality. In this regard, we set the restriction mu_g<mu_{g+1} for g=-L_0+1, -L_0+2,...,L_1.}
#' }
#' For this function, The NET-DPM-3, considered as Fast function is applied , and more details about the algorithm can be referred from Appnendix B.3 of Zhao et al., 2014
#'@return a list cantaining Bayesian Statistics information of the distribution of zi
#'\describe{
#'\item{trace}{a n by (iter-nburns) matrix (n is the length of elements and iter-nburns is the length of the saved chain), showing the evolution of the chain for each element}
#'\item{mean}{the mean of the distribution for each element}
#'\item{median}{the median of the distribution for each element}
#'\item{var}{the variance of the distribution for each element}
#'\item{quantile}{the quantiles of the distribution for each element}
#'}
#'@examples ####Example1. For a 10X10 image with 5X5 signal for example
#' ##Creating the network of 10X10 image
#' library(igraph)
#' library(BayesNetDiscovery)
#' g <- graph.lattice(length=10,dim=2)
#' net=as(get.adjacency(g,attr=NULL),"matrix")##this is the input of argument \code{net}
#' ##Assign the signal elements with signal intenstion as normal distribution N(1,0.2). While noise is set as N(0,0.2) 
#' newz=rep(0,100)
#' for (i in 3:7)
#' {
#'  newz[(i*10+3):(i*10+7)]=1
#' }
#' testcov<-0
#' for(i in 1:100){
#'  if(newz[i]==0){
#'    testcov[i]<-rnorm(1,mean=0,sd=0.2)
#'  }else{
#'   testcov[i]<-rnorm(1,mean=1,sd=0.2)
#'  }
#' }
#' ##The profile of the impage
#' image(matrix(testcov,10,10),col=gray(seq(0,1,length=255)))
#' ##Transform the signals into pvalue form and begin identification
#' pvalue=pnorm(-testcov)
#' total2=Networks.Fast(pvalue,net,iter=5000,initials="mclust",piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(0.5,1,5,10,15))
#' 
#' 
#' 
#' 
#' ####Example.2 Gene Network discovery
#' ##Generating Scale free Gene Network
#' library(igraph)
#' library(BayesNetDiscovery)
#' g <- barabasi.game(50, power=1, zero.appeal=1.5,directed = F)
#' tkplot(g, layout=layout.kamada.kawai)
#' net=as(get.adjacency(g,attr=NULL),"matrix")
#' ##Random assign selected genes and make the signal intension as gaussian mixture
#' newz=rep(c(1,0,0,1,0),10)
#' Simnorm=function(n){
#' weight = c(0.4, 0.6)
#'mu = c(5,4)
#'sigma = c(1,0.5)
#'z = sample(c(1,2),size=n, prob=weight,replace=TRUE)
#'r = rnorm(n,mean=mu[z],sd=sigma[z])
#'return(r)
#'}
#'testcov<-0
#'for(i in 1:50){
#'  if(newz[i]==0){
#'    testcov[i]<-rnorm(1,mean=0,sd=1)
#'  }else{
#'   testcov[i]<-Simnorm(1)
#'  }
#'}
#'pvalue=pnorm(-testcov)
#'total1=Networks.Fast(pvalue,net,iter=5000,v=20,initials="DPdensity",piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(1, 2, 5, 10, 15))
#'@export
Networks.Fast=function(pvalue,net,iter=5000,nburns=2000,initials=c("mclust","DPdensity"),v=20,piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(1, 2, 5, 10, 15),trace=TRUE,show.steps=10)
{
  ####Preparing Stage:Defining Parameters, 
  print("NOW_Transferring p-values to r statistics")
  rstat=Transfer(pvalue)###Step1 : Transferring 
  print("NOW_Getting initials by Kmeans")
  wholeindex=Kmeans(rstat)###Step2: Using Kmeans for giving the initial values
  
  #initial_mu_var=Initial_mu_var(rstat,wholeindex) ###Step3: Genering initial mu and var based on intialindex
  
  pirhopair=Pi_rho_gen(piall,rhoall)####Step4.1 Generating all the possible selections of pi and rho
  print("NOW_Generating_Zi for likelyhood comparison")
  znew=Generating_Zi(net,pirhopair,wholeindex)####Step4.2 Generating Zi by different setting of pi and rho.
  print("NOW_Comparing the likelyhood to select the best set")
  choice=Pi_rho_selecting(znew,rstat,piall,rhoall)####Step 4.3 Using likelyhood to select best set of pi and rho
  #parameter=Parameter_Define(initial_mu_var,pirhopair,choice)####Defining parameters
  
  if(length(choice)!=1){choice=sample(choice,1,prob=c(rep(1/length(choice),length(choice))))}###in case there are multiple choices
  
  print("Interation Begins~")
  if(initials=="mclust")
  {
      mclust=Mclust(rstat)
      if (length(mclust$parameter$mean)==1) {print("warning: the input is not appropriate for mclust since only one cluster was detected by the function Mclust" )}
      if (length(mclust$parameter$mean)==1) break
      hodcmclust=HODCMclust(mclust)        
      total=Iteration3_Mclust(iter,wholeindex,hodcmclust,net,pirhopair,choice,rstat,show.steps,trace)
  }else if (initials=="DPdensity"){
    dpdensitycluster=DPdensitycluster(v,rstat)
     total=Iteration3_DPdensity(iter,wholeindex,dpdensitycluster,net,pirhopair,choice,rstat,v,show.steps,trace)
  }
  
  
  networks.fast=list()
  networks.fast$trace=total[(nburns+1):iter,]
  networks.fast$parameter$mean=sapply(1:length(pvalue), function(kk) return(mean(total[,kk])))
  networks.fast$parameter$median=sapply(1:length(pvalue), function(kk) return(median(total[,kk])))
  networks.fast$parameter$var=sapply(1:length(pvalue), function(kk) return(networks.fast$parameter$mean[kk]*(1-networks.fast$parameter$mean[kk])))
  networks.fast$parameter$quantile=sapply(1:length(pvalue), function(kk) return(quantile(total[,kk])))
  
  if (trace==FALSE){
    show(networks.fast$parameter)       
  }else if (trace==TRUE){
    if (show.steps=="All") {
      show(networks.fast)
    }else{steps=c()
          for (kk in 1:nrow(networks.fast$trace))
          {
            if(kk%%show.steps==0){steps=rbind(steps,total[kk,])}
          }
          showdata=list()
          showdata$steps=steps
          show(showdata)
          show(networks.fast$parameter)
          
    }
  }
  
  
  invisible(networks.fast) 
  
  
  
  
}


#'Hierachical ordered density clustering (HODC) Algorithm with input generated by Mclust
#'@param pvalue a vector of p-values obtained from large scale statistical hypothesis testing
#'@details Without the information of networking, we can have an approximation to the marginal density by DPM model fitting on \strong{r}. Suppose the number of finite mixture normals is equal to L_0+L_1, which means the number of classes we have, we apply HODC algorithm in partitioning the $L_0$ and $L_1$ components into two classes,
#'For this function, the input is generated by Mclust
#'@return a list of HODC algorithm returned parameters. 
#'\describe{
#'\item{mean}{the mean of each of two clusters}
#'\item{variance}{the variance of each of two clusters}
#'\item{pro}{the probability of each of two clusters}
#'\item{classificaiton}{The classification corresponding to each cluster}
#'}
#'@examples rstat=c(rnorm(50,mean=1),rnorm(50,mean=2),rnorm(100,mean=4),rnorm(100,mean=8))
#'pvalue=pnorm(-rstat)
#'mclustHODC=MclustHODC(pvalue)
#'@export
MclustHODC=function(pvalue)
{
  rstat=Transfer(pvalue)
  mclust=Mclust(rstat)

    while (length(mclust$parameter$variance$sigmasq)!=length(unique(mclust$parameter$variance$sigmasq))){
      for (i in 1:length(index)){
        for (j in 1:length(index)){
          if(i!=j){if(mclust$parameter$variance$sigmasq[i]==mclust$parameter$variance$sigmasq[j]){mclust$parameter$variance$sigmasq[j]=mclust$parameter$variance$sigmasq[j]+0.0001}}
        }
      }
    }
    
    while (length(mclust$parameter$pro)!=length(unique(mclust$parameter$pro))){
      for (i in 1:length(index)){
        for (j in 1:length(index)){
          if(i!=j){if(mclust$parameter$pro[i]==mclust$parameter$pro[j]){mclust$parameter$pro[j]=mclust$parameter$pro[j]+0.0001}}
          
        }
      }
    }
    
    if (length(mclust$parameter$mean)==1) {warning("warning: the input is not appropriate for mclust since only one cluster was detected by the function Mclust" )}
    ###Step1 find the min distance
    if (length(mclust$parameter$mean)==2) {
      hodcmclust=list()
      hodcmclust$mean=unique(mclust$parameter$mean)
      hodcmclust$pro=unique(mclust$parameter$pro)
      hodcmclust$variance=unique(mclust$parameter$variance$sigmasq)
    }else{
      repeat{
        distance_all=0
        mclust$parameter$mean=unique(mclust$parameter$mean)
        mclust$parameter$pro=unique(mclust$parameter$pro)
        mclust$parameter$variance$sigmasq=unique(mclust$parameter$variance$sigmasq)
        for (i in 1:(length(mclust$parameter$mean)-1))
        {
          distance=Inte_Distance(i,mclust)
          distance_all=c(distance_all,distance)
        }
        lmin=which(distance_all[-1]==min(distance_all[-1]))
        
        for (l in 1:(length(mclust$parameter$mean)-1))
        {
          if (l<lmin){mclust$parameter$mean[l]=mclust$parameter$mean[l]
                      mclust$parameter$pro[l]=mclust$parameter$pro[l]}
          
          else if (l==lmin){mclust$parameter$mean[l]=mclust$parameter$mean[l]*mclust$parameter$pro[l]/(mclust$parameter$pro[l]+mclust$parameter$pro[l+1])+mclust$parameter$mean[l+1]*mclust$parameter$pro[l+1]/(mclust$parameter$pro[l]+mclust$parameter$pro[l+1])
                            mclust$parameter$pro[l]=mclust$parameter$pro[l]+ mclust$parameter$pro[l+1]
                            k=lmin
                            repeat{
                              k=k+1                           
                              if (length(mclust$classification[which(mclust$classification==k)])!=0){mclust$classification[which(mclust$classification==k)]=lmin
                                                                                                     break}
                            }
                            
                            if (mclust$parameter$variance$modelName!="E"){mclust$parameter$variance$sigmasq[l]=var(rstat[which(mclust$classification==l)])}}
          else if (l>lmin){mclust$parameter$mean[l]=mclust$parameter$mean[l+1]
                           mclust$parameter$variance$sigmasq[l]=mclust$parameter$variance$sigmasq[l+1]
                           mclust$parameter$pro[l]=mclust$parameter$pro[l+1]
          }
          
        }
        
        if (lmin==(length(mclust$parameter$mean)-1)){mclust$parameter$mean=mclust$parameter$mean[-length(mclust$parameter$mean)]
                                                     mclust$parameter$pro=mclust$parameter$pro[-length(mclust$parameter$pro)]
                                                     if (mclust$parameter$variance$modelName!="E"){mclust$parameter$variance$sigmasq=mclust$parameter$variance$sigmasq[-length(mclust$parameter$variance$sigmasq)] }}
        if (length(unique(mclust$parameter$mean))==2) break
      }
      hodcmclust=list()
      index=sort(unique(mclust$classification))
      hodcmclust$mean[1]=mean(rstat[which(mclust$classification==index[1])])
      hodcmclust$mean[2]=mean(rstat[which(mclust$classification==index[2])])
      hodcmclust$variance[1]=var(rstat[which(mclust$classification==index[1])])
      hodcmclust$variance[2]=var(rstat[which(mclust$classification==index[2])])
      hodcmclust$pro[1]=length(rstat[which(mclust$classification==index[1])])/length(rstat)
      hodcmclust$pro[2]=length(rstat[which(mclust$classification==index[2])])/length(rstat)
      hodcmclust$classificaiton=mclust$classification}
    return(hodcmclust) 
}

#'Hierachical ordered density clustering (HODC) Algorithm with input generated by DPdensity
#'@param pvalue a vector of p-values obtained from large scale statistical hypothesis testing
#'@param v number of iterations set for DPM fitting by "DPdensity"
#'@details Without the information of networking, we can have an approximation to the marginal density by DPM model fitting on \strong{r}. Suppose the number of finite mixture normals is equal to L_0+L_1, which means the number of classes we have, we apply HODC algorithm in partitioning the $L_0$ and $L_1$ components into two classes.
#'For this function, the input is generated by Mclust
#'@return a list of HODC algorithm returned parameters.
#'\describe{
#'\item{mean}{the means of each of two cluster for every DPM fitting by "DPdensity"
#'\describe{
#'\item{mu0}{the means of the cluster with smaller mean}
#'\item{mu1}{the means of the cluster with larger mean}
#'}
#'}
#'\item{variance}{the variance of each of two cluster for every DPM fitting by "DPdensity"
#'\describe{
#'\item{var0}{the variances of the cluster with smaller mean}
#'\item{var1}{the variances of the cluster with larger mean}
#'}
#'}
#'\item{probability}{the probability of each of two cluster for every DPM fitting by "DPdensity"
#'\describe{
#'\item{pro0}{the probabilities of the cluster with smaller mean}
#'\item{pro1}{the probabilities of the cluster with larger mean}
#'}
#'}
#'\item{classification}{The classification corresponding to each cluster for every DPM fitting by "DPdensity"}
#'}
#'@examples rstat=c(rnorm(50,mean=1),rnorm(50,mean=2),rnorm(100,mean=4),rnorm(100,mean=8))###random make the density
#'pvalue=pnorm(-rstat)###transformed into pvalue
#'dpdensityHODC=DPdensityHODC(v=5,pvalue)
#'@export
DPdensityHODC<-function(v,pvalue){
  results<-list()
  rstat=Transfer(pvalue)
  
  nburn <-2000
  nsave <-1
  nskip <-0
  ndisplay <-10
  mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
  #prior1 <- list(alpha=0.5,m1=0.2,nu1=100,psiinv1=0.1,k0=100)
  prior2 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
                 psiinv2=solve(diag(0.5,1)),
                 nu1=4,nu2=4,tau1=1,tau2=100)
  fit<-DPdensity(rstat,prior=prior2,mcmc=mcmc,status=TRUE)
  
  for(oo in 1:v){
    cat("iter: ", oo,"\n")
    flush.console()
    nburn <-0
    nsave <-1
    nskip <-0
    ndisplay <-10
    mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
    #prior1 <- list(alpha=0.5,m1=0.2,nu1=100,psiinv1=0.1,k0=100)
    prior2 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
                   psiinv2=solve(diag(0.5,1)),
                   nu1=4,nu2=4,tau1=1,tau2=100)
    fit<-DPdensity(rstat,prior=prior2,mcmc=mcmc,state=fit$state,status=FALSE)
    ##calculate cluster
    ss<-fit$state$ss
    num<-fit$state$ncluster
    prop<-table(fit$state$ss)/num
    mu<-fit$state$muclus[1:num]
    sigma<-fit$state$sigmaclus[1:num]
    index<-order(mu,decreasing=T)
    newm<-mu[index]
    news<-sigma[index] 
    newp<-prop[index]
    
    
    if(num==2 | num==1){results[[oo]]<-ss-1}else{
      #hh1<-1/2/sqrt(pi*news[1])+1/2/sqrt(pi*news[2])-2*dnorm(newm[1],newm[2],sd=sqrt(news[1]+news[2]))
      
      #hh2<-1/2/sqrt(pi*news[2])+1/2/sqrt(pi*news[3])-2*dnorm(newm[2],newm[3],sd=sqrt(news[2]+news[3]))
      
      clmatrix<-matrix(rep(0,num*num),ncol=num,nrow=num)
      for(i in 1:num){
        for(j in 1:num){
          clmatrix[i,j]<-dnorm(newm[i],newm[j],sd=sqrt(news[i]+news[j]))
          
        }
      }
      
      
      temp1<-0
      temp2<-0
      com1<-0
      com2<-0
      res1<-0
      res2<-0
      final1<-0
      final2<-0
      latent<-0
      latent[1]<-0
      final1<-c(1,-1,rep(0,num-2))
      final2<-c(0,-1,1,rep(0,num-3))
      
      for(iii in 1:(num-2)){
        dec1<-t(final1)%*%clmatrix%*%final1
        dec2<-t(final2)%*%clmatrix%*%final2
        
        if(dec1>dec2){latent[iii+1]<-1
        }else{
          latent[iii+1]<-0
        }
        if(latent[iii+1]==1){
          ddd<-max(which(latent==0))
          temp1<-rep(0,num)
          temp1[1:ddd]<-1
          com1<-temp1*newp
          com1<-scale(com1,center=F,scale=sum(com1))
          temp1<-rep(0,num)
          temp1[(ddd+1):(iii+1)]<--1
          res1<-temp1*newp
          res1<-scale(res1,center=F,scale=sum(abs(res1)))
          final1<-com1+res1
          temp2<-rep(0,num)
          temp2[(ddd+1):(iii+1)]<--1
          res2<-temp2*newp
          res2<-scale(res2,center=F,scale=sum(abs(res2)))
          com2<-rep(0,num)
          com2[iii+2]<-1
          final2<-com2+res2
        }else{
          temp1<-rep(0,num)
          temp1[1:iii]<-1
          com1<-temp1*newp
          com1<-scale(com1,center=F,scale=sum(com1))
          res1<-rep(0,num)
          res1[iii+1]<--1
          final1<-com1+res1
          final2<-rep(0,num)
          final2[iii+1]<--1
          final2[iii+2]<-1
        }
        
        
      }
      latent[num]<-1
      
      uuu<-max(which(latent==0))
      ss[ss %in% index[(uuu+1):num]]<-0
      ss[ss %in% index[1:uuu]]<-1
      results[[oo]]<-ss
    }
  }
  dpdensitycluster<-do.call(rbind,results)
  mu0=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==0)])))
  mu1=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==1)])))
  var0=sapply(1:v, function(kk) return(var(rstat[which(dpdensitycluster[kk,]==0)])))
  var1=sapply(1:v, function(kk) return(var(rstat[which(dpdensitycluster[kk,]==1)])))
  pro0=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==0)])/length(pvalue)))
  pro1=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==1)])/length(pvalue)))
  for (kk in 1:v){
    if (is.na(var0[kk])){var0[kk]=var1[kk]
    }else if (is.na(var1[kk])){var1[kk]=var0[kk]}
  }
  DPdensityHODC=list()
  DPdensityHODC$mean$mu0=mu0
  DPdensityHODC$mean$mu1=mu1
  DPdensityHODC$variance$var0=var0
  DPdensityHODC$variance$var1=var1
  DPdensityHODC$probility$pro0=pro0
  DPdensityHODC$probility$pro1=pro1
  DPdensityHODC$classification=dpdensitycluster
  return(DPdensityHODC)
}



