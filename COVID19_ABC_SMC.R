library(tmvtnorm)
source('COVID19_preamble.R')

# outbreak data
agedistr_cases <- 
c(416, 549, 3619, 7600, 8571, 10008, 8583, 3918, 1408)

incidence_cases <- 
c(rep(0,1), 20, 7, 27, 13, 40, 27, 40, 76, 67, 138, 125, 177, 230, 250, 414, 461, 616, 661, 844,
1286, 1484, 1893, 2339, 2609, 2484, 2598, 2553, 2454, 2203, 2143, 2009, 2879, 1744, 1483, 1267, 1108, 730,
622, 436, 292, 160, 83)

tswitch <- 24

# Length of observed outbreak (weeks)
n_obs<-length(incidence_cases)

#  Transmission between age classes
CM<-read.csv("ContactMatrixChina.csv", header=T)
CM<-as.matrix(CM)
A<-matrix(0, nrow=nrow(CM)+2, ncol=ncol(CM)+2)
A[1:nrow(CM), 1:ncol(CM)]<-CM
A[nrow(CM)+1,]<-A[nrow(CM),]
A[,ncol(CM)+1]<-A[,ncol(CM)]
A[nrow(CM)+1, ncol(CM)+1] <- A[nrow(CM), ncol(CM)]
A[nrow(A),]<-A[nrow(A)-1,]
A[,ncol(A)]<-A[,ncol(A)-1]
A[nrow(A), ncol(A)] <- A[nrow(A)-1, ncol(A)-1]

# Age classes
age<-data.frame(lower=seq(0,17*5, by=5), upper=c(seq(5, 17*5, by=5), 100)) 
n_age<-nrow(age)  

#  Adjust observed in eaxh class so that the sum is equal to sum of daily obsrved
age_inc<- sum(incidence_cases)*agedistr_cases/sum(agedistr_cases)

gamma <- 5.6  # incubation period
mu <- 2.9   # infectousness period

# Number of particles
N <- 1000

# Values for threshols 
epsilon_T <- Inf   
epsilon_A <- Inf

# Max number of iterations
G <- 100

# Number of simulations for each parameter set
n <- 1

#  Lower and upper boundaries for priors
# N0, age_sh, age_rt, beta, f_E, delta, age_bg, rho
lm.low<-c(1*10^6,  0.1, 0.1, 0, 0, 0,  0, 0.5)
lm.upp<-c(60*10^6,  20, 20, 1, 1,  1, 1, 1)

run<-1

no  <- sample(1:1000,1)

sim<-0 
while(run <10^6){ 
	res<-read.table(res.file, header=F)
	if(nrow(res)==0){
	  g<-1
	} else {
	  colnames(res)<-c("no", "g", "sim", "N0", "age_sh", "age_rt", "beta", "f_E", "delta", "age_bg", "rho", "w", "D1", "D2")
	  g<- max(res$g)
	  if( g>1) {
	    epsilon_T<-as.numeric(quantile(res$D1[res$g==g-1], prob=0.5))
	    epsilon_A<-as.numeric(quantile(res$D2[res$g==g-1], prob=0.5))
	  }
	  if (nrow(res[res$g==g,])> N) {
	    epsilon_T<-as.numeric(quantile(res$D1[res$g==g], prob=0.5))
	    epsilon_A<-as.numeric(quantile(res$D2[res$g==g], prob=0.5))
		    g<-g+1
	  }
	}	
	
	if(g<= G){
		if(g==1){
    		# Sample from prior distributions 
 			N0_star<- round(runif(1,min=lm.low[1], max=lm.upp[1]))
			age_sh_star<-runif(1, min=lm.low[2], max=lm.upp[2])
			age_rt_star<-runif(1, min=lm.low[3], max=lm.upp[3])
			beta_star<- runif(1, min=lm.low[4], max=lm.upp[4]) 
			f_E_star<- runif(1,min=lm.low[5], max=lm.upp[5]) 
			delta_star<- runif(1,min=lm.low[6], max=lm.upp[6]) 
			age_bg_star<- runif(1,min=lm.low[7], max=lm.upp[7]) 
			rho_star<- runif(1,min=lm.low[8], max=lm.upp[8]) 
			w.old<-rep(1/N,N)
		} else {
			w.old<-res$w[res$g==g-1]/sum(res$w[res$g==g-1])
			res.old<-res[res$g==g-1, c("N0", "age_sh", "age_rt", "beta", "f_E", "delta", "age_bg", "rho")]
			#  Select particle from previous generation
			p<-sample(seq(1, length(w.old)),1,prob=w.old)		
			sigma <- cov(res.old) 
			par<- rK(as.numeric(res.old[p,]),sigma)
			N0_star<-round(par[1])
			age_sh_star<-par[2]
			age_rt_star<-par[3]
			beta_star<-par[4]
			f_E_star<-par[5]
			delta_star<-par[6]
			age_bg_star<-par[7]
			rho_star<-par[8]
		}
      	#  Test if prior non zero
      	if(prior.non.zero(c(N0_star, age_sh_star, age_rt_star, beta_star, f_E_star, delta_star, age_bg_star, rho_star))) {
    			# Set number of accepted simulations to zero
    			m<-0
    			for(j in 1:n){
    				D_star<-run_model(N0_star, age_sh_star, age_rt_star, beta_star, f_E_star, delta_star, age_bg_star, rho_star)     
    				sim<-sim+1
    				# Calculate distances 
    				dist<-calc_distance(D_star[[1]], D_star[[2]])
   
    				if((dist[1] <= epsilon_T) & (dist[2] <= epsilon_A)){ # If all distances are less than their tolerances
    					m<-m+1
    				}
    			}	
    			if (m>0){
    				res.new<-c(N0_star, age_sh_star, age_rt_star, beta_star, f_E_star, delta_star, age_bg_star, rho_star)  
      			w1<-prod(sapply(1:7, function(b) dunif(res.new[b], min=lm.low[b], max=lm.upp[b])))
				if(g==1){
					w2<-1
				} else {
					w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new, mean=as.numeric(res.old[a,]), sigma=sigma, lower=lm.low, upper=lm.upp)))
				}
      			w.new <- (m/n)*w1/w2
       			cat(c(no, g, sim, res.new, w.new, dist, "\n"), file=res.file, append = TRUE, sep="\t")
      			cat(1, file=wrfl, append = F)
      			sim<-0 
      		} 
		}
	} else {
	run<-10^8}
	run<-run+1
} 

