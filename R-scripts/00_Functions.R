

##########
###### Briere model (truncated) ----
##########

sink("briere_T.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(20, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()


##########
###### Briere model (with random effects) ----
##########

sink("briere_T_randeff.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(20, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Random effect for q
    sigma_q ~ dunif(0, 1000) # standard deviation of random effect (variance between unique ids)
	  tau_q <- 1 / (sigma_q * sigma_q) # convert to precision
	  
    for (j in 1:Nids){
		q[j] ~ dnorm(0, tau_q) # random q for each unique id
		}

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- (cf.q + q[unique.id[i]]) * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- (cf.q + q[unique.id[i]]) * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()

##########
###### Briere Model with gamma priors (except sigma) ----
##########

sink("briere_inf.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### Quadratic model (truncated) ----
##########

sink("quad_T.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(20, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()

##########
###### Quadratic Model with gamma priors (except sigma) ----
##########

sink("quad_inf.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()
