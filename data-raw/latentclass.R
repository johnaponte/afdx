library(afdx)
library(coda)
library(rjags)
library(ggmcmc)

cutoffs <- c(0,1,100,200,400,800,1600,3200,6400,12800,25600,51200, 102400, 204800)
data <- 
  malaria_df1 %>%
  mutate(k = cut(density,c(cutoffs,Inf), include.lowest =T, labels = cutoffs)) %>%
  group_by(k,fever) %>%
  tally() %>%
  mutate(category = ifelse(fever ==1,"n (fever)","m (no fever)")) %>%
  select(-fever) %>%
  pivot_wider(names_from = "category", values_from = "n", values_fill = list(n = 0)) %>%
  rename(`k (category lower limit)` = k)

model <- get_latent_model()

# compile the model
af_latent <-
  jags.model(
    textConnection(get_latent_model()),
    data = list(n = data$`n (fever)`,
                m = data$`m (no fever)`),
    n.chains = 4,
    n.adapt = 1000,
    inits = list(
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1111),
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2222),
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 3333),
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 4444)
    )
  )

# Simulate the posterior
latent_sim <-
  coda.samples(
    model = af_latent,
    variable.names = c('lambda','sens','spec','ppv','npv'),
    n.thinning = 5,
    n.iter =  10000 )

latent_sum <-  summary(latent_sim)
latent_eff <-  effectiveSize(latent_sim)


saveRDS(latent_sum, file = "inst/vignette_data/latent_sum.RDS", compress = "gzip")
saveRDS(latent_eff, file = "inst/vignette_data/latent_eff.RDS", compress = "gzip")