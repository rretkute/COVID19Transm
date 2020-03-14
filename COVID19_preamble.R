
run_model<-function(N0, age_sh, age_rt, beta, f_E, delta, age_bg, rho){
	
	# Number of individuals in each class
	dist_pop<-pbeta(age$upper/max(age$upper), age_sh, age_rt)
	dist_pop<-c(dist_pop[1], sapply(2:(length(dist_pop)), function(a) dist_pop[a]-dist_pop[a-1])) + age_bg
	dist_pop<-round(N0*dist_pop/sum(dist_pop)) 
	
	# Define the initial conditions
	t0 <- 0; 
	I0 <- rep(0,n_age)
	E0<-pmax(round(f_E*dist_pop),1)
	R0 <- rep(0, n_age)
	S0<-pmax(dist_pop-I0-E0-R0,1)
	
	t<-t0; S<-S0; E<-E0; I<-I0; R<-R0
	
	cases<-rep(0,t_end)
	cases_age<-rep(0,n_age)
	
	while (t<t_end) {
		t<-t+1
		S_old<-S; E_old<-E; I_old<-I; R_old<-R; 
		# S->E  
		if(t<tswitch){
			lambda <- beta*( A %*% I_old)
		} else {
			lambda <- exp(-delta*(t-tswitch))*beta*( A %*% I_old)
		}
		n.exp<-sapply(1:n_age, function(a) rbinom(1, S_old[a], 1-exp(-lambda[a])))
		S<-pmax(S-n.exp, 0); E<-E+n.exp
		# E->I  
		n.inf<-sapply(1:n_age, function(a) rbinom(1, E_old[a],1-exp(-1/gamma)))
		E <- pmax(E-n.inf,0); I <- I+n.inf
		reported.cases<-sapply(1:n_age, function(a) rbinom(1, n.inf[a], rho))
		cases[t] <- sum(reported.cases)
		cases_age <- cases_age+reported.cases
		# I->R  
		n.rem<-sapply(1:n_age, function(a) rbinom(1, I_old[a], 1-exp(-(1/mu))))
		I <- pmax(I-n.rem, 0); R <-R+n.rem
	}
	ac<-sapply(1:length(cl.to.test), function(a) sum(cases_age[cl.to.test[[a]]]))
	return(list(cases, ac))
}

calc_distance <- function(cases, cases_age){
  sre1<-sqrt(sum((cases[1:(n_obs)] - incidence_cases[1:n_obs])^2)) 
  if(sum(cases_age)>0) {
    sre2<-sqrt(sum((age_inc_perc-cases_age)^2))
  } else { 	
    sre2<-sqrt(sum((age_inc_perc)^2)) 
  }
  return(c(sre1, sre2))
}


# Perturbation kernel 
rK <- function(mean, sigma){   
	return(rtmvnorm(1,mean=mean, sigma=sigma, lower=lm.low, upper=lm.upp)) 
}

#  Identity function
H <- function(x) as.numeric(x>0)

#  Test if prior is non zero
prior.non.zero<-function(par){
	prod(sapply(1:length(par), function(a) H(par[a]-lm.low[a])* H(lm.upp[a]-par[a])))
}
		    


