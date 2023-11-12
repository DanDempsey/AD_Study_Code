##### Benchmark Test using Results from UPG package vignette
##### Daniel Dempsey

### Load in inference code and libraries
setwd( 'Code/Main_Software' )
source( 'Polychotomous_VS_Implementation.R' )
dyn.load( 'KS_RNG.so' )
library( tictoc ) # For timing runs 
library( UPG ) # Multinomial regression package

### Load in data and prep response variable
data( 'program' )
y <- factor( program$program )

### MCMC settings
MCMC_len <- 5000
privar <- 4
gam_start <- rep( TRUE, 4 )
burn_in <- 2500
baseline <- levels( y )[ 1 ]

### HH fit (no VS)
set.seed( 111 )
tic()
program_fit <- bayes_logistic( y, program[ , 3:5 ], MCMC_iterations = MCMC_len, 
                               prior_variance = privar, gamma_start = gam_start, 
                               gamma_update_length = MCMC_len + 1, burn_in = burn_in,
                               thin = 1, baseline = baseline, scaling = FALSE )
toc() # 60.53 sec elapsed

### HH fit (with VS)
set.seed( 222 )
tic()
program_fit_VS <- bayes_logistic( y, program[ , 3:5 ], MCMC_iterations = MCMC_len, 
                                  prior_variance = privar, gamma_start = gam_start,
                                  burn_in = burn_in, thin = 1, baseline = baseline )
toc() # 67.495 sec elapsed

### UPG fit
set.seed( 333 )
tic()
upg_fit <- UPG( y = y, X = program[-1], model = 'mnl', A0 = privar,
                B0 = privar, verbose = FALSE, baseline = baseline,
                draw = MCMC_len, burnin = burn_in, gamma.boost = FALSE, 
                delta.boost = FALSE )
toc() # 23.88 sec elapsed

### Compare results
keep <- (burn_in + 1):MCMC_len
summary( program_fit )
summary( program_fit_VS )
summary( upg_fit )

#### Examine how intercept changes
# female and ses removed
set.seed( 444 )
upg_fit2 <- UPG( y = y, X = program[-c(1, 3:4)], model = 'mnl', A0 = privar,
                B0 = privar, verbose = FALSE, baseline = baseline,
                draw = MCMC_len, burnin = burn_in, gamma.boost = FALSE, 
                delta.boost = FALSE )

summary( upg_fit2 )

bb <- program_fit_VS$beta$general[keep, 1]
plot( density( bb ) )


gg <- program_fit_VS$gamma$general[ keep, 'ses' ]
bb_with_ses <- program_fit_VS$beta$general[keep, 1][ gg == TRUE ]
bb_without_ses <- program_fit_VS$beta$general[keep, 1][ gg == FALSE ]

plot( density( bb_with_ses ) )
mean( bb_with_ses )

plot( density( bb_without_ses ) )
mean( bb_without_ses )

