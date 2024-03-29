% Generated by roxygen2 (4.0.0): do not edit by hand
\name{Networks.Fast}
\alias{Networks.Fast}
\title{Bayesian Network Discovery using a Hybrid Fast Algorithm}
\description{
This function implements a hybrid fast algorithm to perform feature selection and sub-network discovery using a Bayesian nonparametric model based on Dirichlet process mixture models, finite mixture of normals model and the Ising model.  

}


\usage{
Networks.Fast(pvalue, net, iter = 5000, nburns = 2000, 
initials=c("mclust","DPdensity"), v = 20,
piall = c(0.8, 0.85, 0.9, 0.95), rhoall = c(1, 2, 5, 10, 15),
trace=TRUE,show.steps=10)
}
\arguments{
\item{pvalue}{a vector of p-values obtained from large scale statistical hypothesis testing}

\item{net}{an "n"" by "n" binary (0/1) adjacent matrix of network configurations, where n=length(pvalue)}

\item{iter}{number of iterations; the default is 5000}

\item{nburns}{number of burn-in; the default is 2000}

\item{initials}{character taking value "DPdensity" or "mclust" indiciating the function to be used to obtain initial density estimates. It is recommended to choose "Mclust" when the dimension of data is large, and to choose "DPdensity" when the dimenstion is small.}

\item{v}{number of iterations set for DPM fitting by "DPdensity". v is only valid when you choose initals as "DPdensity"}

\item{piall}{a vector of possible choices for "pi0" in an increasing order; the default value is c(0.75, 0.8, 0.85, 0.9)}

\item{rhoall}{a vector of possible choices of "rho0" and "rho1" in an increasing order; the default value is c(0.5, 1, 5, 10, 15)
}
\item{trace}{a logical variable indicating whether this run showing the results of iterations (\code{TRUE}), or not (\code{FALSE}), the default setting is \code{TRUE}}
\item{show.steps}{integer representing the frequency of the results of iterations presented, the defualt setting is 10. The setting is invalid when \code{trace=FALSE}. The setting would not affect the data saved, only for printing}
}

\value{a list of values containing the posterior samples of "z_i" and summary statistics
\describe{
\item{trace}{an length(pvalue) by (iter-nburns) matrix}
\item{parameter}{a list of summary statistics characterizing the posterior distribution of "z_i"
\describe{
\item{mean}{posterior mean for each feature}
\item{median}{posterior median  for each feature}
\item{var}{posterior variance for each feature}
\item{quantile}{posterior quantiles for each feature}}}
}
}

\details{
This function implements a Bayesian nonparametric mixture model for feature selection incorporating network information (Zhao et al., 2014):
\itemize{
\item  r_i| g_i, \strong{theta} ~ N(mu_{g_i}, sigma_{g_i}),
\item  g_i | z_i=k, \strong{q}_k ~ Discrete(\strong{a}_k, \strong{q}_k),
\item	\strong{theta} ~ G_{0k}, for g in \strong{a}_k,
\item	\strong{q}_k ~ Dirichlet(tau_k \strong{1}_{\{L_k\}}/L_k),
\item	\strong{theta}=\{{\strong{theta}_g}\}_{\{g in \strong{a}_0 and \strong{a}_1\}}
\item	\strong{theta}_g=(mu_g, sigma_g)
}
where we define
\describe{
\item{Index}{\strong{a}_0=(-L_0+1,-L_0+2,...,0) , \strong{a}_1=(1,2,...,L_1) and the correspondent probability q_0=(q_{-L_0+1}, q_{-L_0+2}, ...,q_0), q_1=(q_1, q_2, ..., q_{L_1}), according to the defination of Discrete(\strong{a}_k, \strong{b}_k), for example, Pr(g_i={L_0+2})=q_{-L_0+2}. }
\item{Assumption}{In this algorithm, we assume that "important" features should have larger statiscs comparing to "unimportant" ones without the loss of generality. In this regard, we set the restriction mu_g<mu_{g+1} for g=-L_0+1, -L_0+2,...,L_1.}
}
This function implements the NET-DPM-3 Zhat et al.(2014). Please refer to the Appnendix B.3 for more details. 
}

\references{
Zhao, Y., Kang, J., Yu, T. A Bayesian nonparametric mixture model for selecting gene and gene-sub network, Annals of Applied Statistics, In press: 2014.

Zhou Lan, Jian Kang, Tianwei Yu, Yize Zhao, BayesNetDiscovery: an R package for network identifications via Bayesian nonparametric mixture models, working paper.

}

\author{
Zhou Lan, Jian Kang, Tianwei Yu and Yize Zhao
}



\examples{
####Gene Network discovery
##Generating Scale free Gene Network
library(igraph)
library(BayesNetDiscovery)
g <- barabasi.game(50, power=1, zero.appeal=1.5,directed = FALSE)
net=as(get.adjacency(g,attr=NULL),"matrix")
##Random assign selected genes and make the signal intension as gaussian mixture
newz=rep(c(1,0,0,1,0),10)
Simnorm=function(n){
weight = c(0.4, 0.6)
mu = c(8,6)
sigma = c(1,0.5)
z = sample(c(1,2),size=n, prob=weight,replace=TRUE)
r = rnorm(n,mean=mu[z],sd=sigma[z])
return(r)
}
testcov<-0
for(i in 1:50){
 if(newz[i]==0){
   testcov[i]<-rnorm(1,mean=0,sd=1)
 }else{
  testcov[i]<-Simnorm(1)
 }
}
pvalue=pnorm(-testcov)
total1=Networks.Fast(pvalue,net,iter=10,nburns=2,
v=20,initials="DPdensity",
piall=c(0.8, 0.85, 0.9, 0.95),rhoall=c(1, 2, 5, 10, 15),
trace=FALSE)
}

