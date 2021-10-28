library(rjags)

moths <- data.frame(site = rep(1:7, each = 2),
                    morph = rep(1:2, times = 7),
                    distance = rep(c(0, 7.2, 24.1, 30.2, 36.4, 41.5, 51.2), each = 2),
                    placed = rep(c(56, 80, 52, 60, 60, 84, 92), each = 2),
                    removed = c(13, 14, 28, 20, 18, 22, 9, 16, 16, 23, 20, 40, 24, 39))

datlist<-list(y=moths$removed,
              x=moths$distance,
              m=moths$morph-1,
              #s=moths$site,
              n=moths$placed,
              n_row=nrow(moths))

mscript<-"model{
  for(i in 1:n_row){
    y[i] ~ dbin(p[i], n[i])
    logit(p[i]) <- beta1 + beta2 * x[i] + beta3*m[i] + beta4*x[i]*m[i]
  }
  
  beta1 ~ dnorm(0, 1/10^2)
  beta2 ~ dnorm(0, 1/5^2)
  beta3 ~ dnorm(0, 1/5^2)
  beta4 ~ dnorm(0, 1/5^2)
  
  beta2_morph1 <- beta2
  beta2_morph2 <- beta2 + beta4
  
}"

jm <- jags.model(file=textConnection(mscript),
                 data=datlist,
                 n.chains=3)

jm_coda <- coda.samples(jm,
                        variable.names=c("p",
                                         "beta1",
                                         "beta2",
                                         "beta3",
                                         "beta4",
                                         "beta2_morph1",
                                         "beta2_morph2"),
                        n.iter=3000)

traplot(jm_coda,
        parms = c("p",
                  "beta1",
                  "beta2",
                  "beta3",
                  "beta4",
                  "beta2_morph1",
                  "beta2_morph2"))

mcmcplot(jm_coda,
         parms=c("beta1",
                 "beta2",
                 "beta3",
                 "beta4"))

jm_coda <- coda.samples(jm,
                        variable.names=c("p",
                                         "beta1",
                                         "beta2",
                                         "beta3",
                                         "beta4",
                                         "beta2_morph1",
                                         "beta2_morph2"),
                        n.iter=3000*40,
                        thin=40)

mcmcplot(jm_coda,
         parms=c("beta1",
                 "beta2",
                 "beta3",
                 "beta4"))
# Rhat, or Gelman diagnostic
gelman.diag(jm_coda, multivariate=FALSE)

caterplot(jm_coda, 
          parms = c("beta1",
                    "beta2",
                    "beta3",
                    "beta4",
                    "beta2_morph1",
                    "beta2_morph2"),
          reorder = FALSE)

caterplot(jm_coda, 
          parms = c("beta2_morph1",
                    "beta2_morph2"),
          reorder = FALSE)

summary(jm_coda)

mbrm<-brm(removed | trials(placed) ~ morph*distance, family=binomial,
          data=mutate(moths, morph=as.factor(morph)),
          prior=prior(normal(0,10), class="Intercept")+
                        prior(normal(0,5), class="b"),
          iter=5e3, warmup=5e3/4)

conditional_effects(mbrm)
summary(mbrm)

mbrm2<-brm(removed | trials(placed) ~ morph*distance + (1|site), family=binomial,
          data=mutate(moths, morph=as.factor(morph)),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sd"),
          iter=5e3, warmup=5e3/4, control=list(adapt_delta=0.9))

conditional_effects(mbrm2)
