

##########
###### Briere model (truncated) ----
##########

sink("R-scripts/briere_T.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
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
###### Briere model Probability (truncated) ----
##########

sink("R-scripts/briereprob.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i]) * (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) < 1) + (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) > 1)
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### Briere model (with random effects) ----
##########

sink("R-scripts/briere_T_randeff.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)  # tau is precision (1 / variance)
    
    ## Random effect priors
    sigma_q ~ dunif(prior[1,4], prior[2,4])
    tau_q <- 1 / (sigma_q * sigma_q)
    
    sigma_T0 ~ dunif(prior[1,5], prior[2,5])
    tau_T0 <- 1 / (sigma_T0 * sigma_T0)
    
    sigma_Tm ~ dunif(prior[1,6], prior[2,6])
    tau_Tm <- 1 / (sigma_Tm * sigma_Tm)
    
    ## Random effects for each species-study combination (unique_id)
     
    for (j in 1:Nids) {
    q[j] ~ dnorm(0, tau_q)
    T0[j] ~ dnorm(0, tau_T0)
    Tm[j] ~ dnorm(0, tau_Tm)
    }
		
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- (cf.q + q[unique.id[i]]) * temp[i] * (temp[i] - (cf.T0 + T0[unique.id[i]])) * sqrt(((cf.Tm + Tm[unique.id[i]]) - temp[i]) * ((cf.Tm + Tm[unique.id[i]]) > temp[i])) * ((cf.T0 + T0[unique.id[i]]) < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred.pop[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])}
    
    for (j in 1:Nids) {
      for(i in 1:N.Temp.xs){
        z.trait.mu.pred.id[j,i] <- (cf.q + q[j]) * Temp.xs[i] * (Temp.xs[i] - (cf.T0 + T0[j])) * sqrt(((cf.Tm + Tm[j]) - Temp.xs[i]) * ((cf.Tm + Tm[j]) > Temp.xs[i])) * ((cf.T0 + T0[j]) < Temp.xs[i])
      }
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### Briere model Probability (truncated) (with random effects) ----
##########

sink("R-scripts/briereprob_randeff.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)  # tau is precision (1 / variance)
    
    ## Random effect priors
    sigma_q ~ dunif(prior[1,4], prior[2,4])
    tau_q <- 1 / (sigma_q * sigma_q)
    
    sigma_T0 ~ dunif(prior[1,5], prior[2,5])
    tau_T0 <- 1 / (sigma_T0 * sigma_T0)
    
    sigma_Tm ~ dunif(prior[1,6], prior[2,6])
    tau_Tm <- 1 / (sigma_Tm * sigma_Tm)
    
    ## Random effects for each species-study combination (unique_id)
     
    for (j in 1:Nids) {
    q[j] ~ dnorm(0, tau_q)
    T0[j] ~ dnorm(0, tau_T0)
    Tm[j] ~ dnorm(0, tau_Tm)
    }
		
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- (cf.q + q[unique.id[i]]) * temp[i] * (temp[i] - (cf.T0 + T0[unique.id[i]])) * sqrt(((cf.Tm + Tm[unique.id[i]]) - temp[i]) * ((cf.Tm + Tm[unique.id[i]]) > temp[i])) * ((cf.T0 + T0[unique.id[i]]) < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred.pop[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i]) * (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) < 1) + (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) > 1)
    }
  
    
    for (j in 1:Nids) {
      for(i in 1:N.Temp.xs){
          z.trait.mu.pred.id[j,i] <- (cf.q + q[j]) * Temp.xs[i] * (Temp.xs[i] - (cf.T0 + T0[j])) * sqrt(((cf.Tm + Tm[j]) - Temp.xs[i]) * ((cf.Tm + Tm[j]) > Temp.xs[i])) * ((cf.T0 + T0[j]) < Temp.xs[i]) * ((cf.q + q[j]) * Temp.xs[i] * (Temp.xs[i] - (cf.T0 + T0[j])) * sqrt(((cf.Tm + Tm[j]) - Temp.xs[i]) * ((cf.Tm + Tm[j]) > Temp.xs[i])) < 1) + ((cf.q + q[j]) * Temp.xs[i] * (Temp.xs[i] - (cf.T0 + T0[j])) * sqrt(((cf.Tm + Tm[j]) - Temp.xs[i]) * ((cf.Tm + Tm[j]) > Temp.xs[i])) > 1)
      }
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### Briere Model with gamma priors (except sigma) ----
##########

sink("R-scripts/briere_inf.txt")
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
###### Briere Model Probability with gamma priors (except sigma) ----
##########

sink("R-scripts/briereprob_inf.txt")
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
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i]) * (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) < 1) + (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) > 1)
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### Briere Model with random effect and gamma priors (except sigma) ----
##########

sink("R-scripts/briere_inf_raneff.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    
    ## Random effect priors
    sigma_q ~ dunif(0, 0.001)
    tau_q <- 1 / (sigma_q * sigma_q)
    
    sigma_T0 ~ dunif(0, 10)
    tau_T0 <- 1 / (sigma_T0 * sigma_T0)
    
    sigma_Tm ~ dunif(0, 10)
    tau_Tm <- 1 / (sigma_Tm * sigma_Tm)
    
    ## Random effects for each species-study combination (unique_id)
     
    for (j in 1:Nids) {
    q[j] ~ dnorm(0, tau_q)
    T0[j] ~ dnorm(0, tau_T0)
    Tm[j] ~ dnorm(0, tau_Tm)
    }
		
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- (cf.q + q[unique.id[i]]) * temp[i] * (temp[i] - (cf.T0 + T0[unique.id[i]])) * sqrt(((cf.Tm + Tm[unique.id[i]]) - temp[i]) * ((cf.Tm + Tm[unique.id[i]]) > temp[i])) * ((cf.T0 + T0[unique.id[i]]) < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred.pop[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])}
    
    for (j in 1:Nids) {
      for(i in 1:N.Temp.xs){
        z.trait.mu.pred.id[j,i] <- (cf.q + q[j]) * Temp.xs[i] * (Temp.xs[i] - (cf.T0 + T0[j])) * sqrt(((cf.Tm + Tm[j]) - Temp.xs[i]) * ((cf.Tm + Tm[j]) > Temp.xs[i])) * ((cf.T0 + T0[j]) < Temp.xs[i])
      }
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### Quadratic model (truncated) ----
##########

sink("R-scripts/quad_T.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
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
###### Quadratic model (with random effects) ----
##########

sink("R-scripts/quad_T_randeff.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Random effect priors
    sigma_q ~ dunif(prior[1,4], prior[2,4])
    tau_q <- 1 / (sigma_q * sigma_q)
    
    sigma_T0 ~ dunif(prior[1,5], prior[2,5])
    tau_T0 <- 1 / (sigma_T0 * sigma_T0)
    
    sigma_Tm ~ dunif(prior[1,6], prior[2,6])
    tau_Tm <- 1 / (sigma_Tm * sigma_Tm)
    
    ## Random effects for each species-study combination (unique_id)
     
    for (j in 1:Nids) {
    q[j] ~ dnorm(0, tau_q)
    T0[j] ~ dnorm(0, tau_T0)
    Tm[j] ~ dnorm(0, tau_Tm)
    }
		
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * (cf.q + q[unique.id[i]]) * (temp[i] - (cf.T0 + T0[unique.id[i]])) * (temp[i] - (cf.Tm + Tm[unique.id[i]])) * ((cf.Tm + Tm[unique.id[i]]) > temp[i]) * ((cf.T0 + T0[unique.id[i]]) < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred.pop[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    for (j in 1:Nids) {
      for(i in 1:N.Temp.xs){
        z.trait.mu.pred.id[j,i] <- -1 * (cf.q + q[j]) * (Temp.xs[i] - (cf.T0 + T0[j])) * (Temp.xs[i] - (cf.Tm + Tm[j])) * ((cf.Tm + Tm[j]) > Temp.xs[i]) * ((cf.T0 + T0[j]) < Temp.xs[i])
      }
    }
    } # close model
    ",fill=T)
sink()



##########
###### Quadratic Model with gamma priors (except sigma) ----
##########

sink("R-scripts/quad_inf.txt")
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



##########
###### Quadratic model Probability (truncated) ----
##########
## For proportion and probability traits: maximum trait value set to 1

sink("R-scripts/quadprob.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### Quadratic model Probability (with random effects) ----
##########

sink("R-scripts/quadprob_randeff.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(prior[1,1], prior[2,1])
    cf.T0 ~ dunif(prior[1,2], prior[2,2])
    cf.Tm ~ dunif(prior[1,3], prior[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Random effect priors
    sigma_q ~ dunif(prior[1,4], prior[2,4])
    tau_q <- 1 / (sigma_q * sigma_q)
    
    sigma_T0 ~ dunif(prior[1,5], prior[2,5])
    tau_T0 <- 1 / (sigma_T0 * sigma_T0)
    
    sigma_Tm ~ dunif(prior[1,6], prior[2,6])
    tau_Tm <- 1 / (sigma_Tm * sigma_Tm)
    
    ## Random effects for each species-study combination (unique_id)
     
    for (j in 1:Nids) {
    q[j] ~ dnorm(0, tau_q)
    T0[j] ~ dnorm(0, tau_T0)
    Tm[j] ~ dnorm(0, tau_Tm)
    }
		
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * (cf.q + q[unique.id[i]]) * (temp[i] - (cf.T0 + T0[unique.id[i]])) * (temp[i] - (cf.Tm + Tm[unique.id[i]])) * ((cf.Tm + Tm[unique.id[i]]) > temp[i]) * ((cf.T0 + T0[unique.id[i]]) < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred.pop[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
    }
    
    for (j in 1:Nids) {
      for(i in 1:N.Temp.xs){
        z.trait.mu.pred.id[j,i] <- -1 * (cf.q + q[j]) * (Temp.xs[i] - (cf.T0 + T0[j])) * (Temp.xs[i] - (cf.Tm + Tm[j])) * ((cf.Tm + Tm[j]) > Temp.xs[i]) * ((cf.T0 + T0[j]) < Temp.xs[i]) * (-1 * (cf.q + q[j]) * (Temp.xs[i] - (cf.T0 + T0[j])) * (Temp.xs[i] - (cf.Tm + Tm[j])) < 1) + (-1 * (cf.q + q[j]) * (Temp.xs[i] - (cf.T0 + T0[j])) * (Temp.xs[i] - (cf.Tm + Tm[j])) > 1)
      }
    }
    } # close model
    ",fill=T)
sink()



##########
###### Quadratic Model with gamma priors (except sigma) ----
##########

sink("R-scripts/quadprob_inf.txt")
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
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### Quadratic Model with random effect and gamma priors (except sigma) ----
##########

sink("R-scripts/quad_inf_raneff.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    
    ## Random effect priors
    sigma_q ~ dunif(0, 0.01)
    tau_q <- 1 / (sigma_q * sigma_q)
    
    sigma_T0 ~ dunif(0, 10)
    tau_T0 <- 1 / (sigma_T0 * sigma_T0)
    
    sigma_Tm ~ dunif(0, 10)
    tau_Tm <- 1 / (sigma_Tm * sigma_Tm)
    
    ## Random effects for each species-study combination (unique_id)
     
    for (j in 1:Nids) {
    q[j] ~ dnorm(0, tau_q)
    T0[j] ~ dnorm(0, tau_T0)
    Tm[j] ~ dnorm(0, tau_Tm)
    }
		
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * (cf.q + q[unique.id[i]]) * (temp[i] - (cf.T0 + T0[unique.id[i]])) * (temp[i] - (cf.Tm + Tm[unique.id[i]])) * ((cf.Tm + Tm[unique.id[i]]) > temp[i]) * ((cf.T0 + T0[unique.id[i]]) < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred.pop[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    for (j in 1:Nids) {
      for(i in 1:N.Temp.xs){
        z.trait.mu.pred.id[j,i] <- -1 * (cf.q + q[j]) * (Temp.xs[i] - (cf.T0 + T0[j])) * (Temp.xs[i] - (cf.Tm + Tm[j])) * ((cf.Tm + Tm[j]) > Temp.xs[i]) * ((cf.T0 + T0[j]) < Temp.xs[i])
      }
    }
    
    } # close model
    ",fill=T)
sink()
