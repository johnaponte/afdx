library(rjags)
library(afdx)
#' Estimate the AF for the malaria_df1 dataset using different categories
estimate_af <- function(ncat) {
  model_af = afdx::get_latent_model()
  model_aft = substring(model_af, 114) # removing the part that defines the data
  cutpoints = c(0,
                with(
                    malaria_df1,quantile(
                      density[density > 0], 
                      probs = seq(0, 1, 1/ncat))))
  K = length(cutpoints)-1
  cat = dplyr::mutate(malaria_df1, cats = cut(density, breaks = cutpoints) )
  df = data.frame(table(cat$cats,cat$fever))
  n=df$Freq[1:K]
  n[is.na(n)] <- 0
  m=df$Freq[K+1:K*2]
  m[is.na(m)] <- 0
  data_af = list(n=n,m=m,K=K,Sn=sum(n),Sm=sum(m))
  requiredVars = c('lambda','sens','spec')
  jags_af = rjags::jags.model(
    textConnection(model_aft), 
    data=data_af, 
    inits = list(.RNG.name = "base::Wichmann-Hill"), 
    n.adapt=1000)
  jagssamples_af = coda.samples(jags_af,
                                    variable.names = requiredVars,
                                    n.iter = 10000, 
                                    n.burnin=2000,
                                    n.thinning = 5)
  
  stats_af <- summary(jagssamples_af)
  # Return the lambda row
  data.frame(t(cbind(stats_af[[1]], stats_af[[2]])["lambda",]))
}

library(plyr)
af_cats <- ldply(5:25, estimate_af)
af_cats$ncats = 5:25

library(ggplot2)
ggplot(
  af_cats,
  aes(x=ncats, y = Mean, ymin = X2.5., ymax = X97.5.)
) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  scale_y_continuous("Estimated Attributable fraction") +
  scale_x_continuous("Number of categories of positive values") +
  ggtitle("Dependency of the AF estimation with the number of categories")


