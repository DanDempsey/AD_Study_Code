##### Simulation to test Bayesian logistic model inference code
##### Daniel Dempsey

### Load in inference code and libraries
source( 'Code/Main_Software/Polychotomous_VS_Implementation.R' )
dyn.load( 'Code/Main_Software/KS_RNG.so' )
setwd( 'Output/Simulation_Output/' )
library( tictoc ) # For timing runs 
library( parallel ) # For parallel processing
library( corrplot ) # Plotting correlation matrices

### Convience wrapper for looping over MCMC algorithm
run_fit <- function( reps, dat, ... ) {
  bayes_logistic( X = dat, ... )
}

### Create function that performs a simulation experiment
sim_function <- function( reps, loops = 10, folder, N, beta_sim, n_cores, ... ) {
 
  tic()
  
  res <- vector( 'list', loops )
  for ( i in 1:loops ) {
    
    # Simulate covariates, including intercept
    G <- ncol( beta_sim )
    nvars <- nrow( beta_sim ) - 1
    X <- matrix( c( rep( 1, N ), rnorm( nvars*N ) ), nrow = N )
    
    # Simulate linear predictor
    eta <- cbind( V0 = 0, as.data.frame( lapply( beta_sim, function ( b ) { X%*%b } ) ) )
    
    # Simulate response variable
    multiprob <- pmultinom( eta )
    y <- apply( multiprob, 1, rmultinom2 )
    
    # Run model
    res[[i]] <- mclapply( reps, run_fit, y = y, dat = X[ , -1 ], 
                          baseline = 'V0', ..., mc.cores = n_cores )
    cat( paste0( 'Iteration ', i, ' completed!\n' ) )
    
  }
  
  cat( 'Saving results and creating visualisations...\n' )
  
  # Save results
  save( res, file = paste0(folder, '/res.Rdata') )
  
  ### Visualisations
  # Plots gamma results
  at_ind <- rev(seq( nrow(beta_sim) - 1 ))
  for ( i in 1:length(res) ) {
    
    V1_sel <- sel_fun( res[[i]], 1 ) 
    V2_sel <- sel_fun( res[[i]], 2 )
    
    colnames( V1_sel ) <- identify_par( beta_sim$V1[ -1 ] )
    colnames( V2_sel ) <- identify_par( beta_sim$V2[ -1 ] )
    
    # Group 1
    pdf( paste0(folder, '/gammaV1_', i, '.pdf'), height = 10, width = 10 )
    boxplot( V1_sel, border = col_fun( V1_sel ), main = 'Group 1', 
             ylim = c( 0, 1 ), horizontal = TRUE, las = 2, at = at_ind,
             xlab = 'Proportion Selected' )
    grid()
    abline( v = 0.5, lty = 2 )
    dev.off()
    
    # Group 2
    pdf( paste0(folder, '/gammaV2_', i, '.pdf'), height = 10, width = 10 )
    boxplot( V2_sel, border = col_fun( V2_sel ), main = 'Group 2', 
             ylim = c( 0, 1 ), horizontal = TRUE, las = 2, at = at_ind,
             xlab = 'Proportion Selected' )
    grid()
    abline( v = 0.5, lty = 2 )
    dev.off()
    
  }
  
  cat( 'Complete!\n' )
  
  toc()
  
}  

### Visualisation Functions
sel_fun <- function( dat, ind ) {
  
  keep <- seq( dat[[1]]$burn_in + 1, dat[[1]]$MCMC_iterations, dat[[1]]$thin )
  res <- lapply( dat, function(x) { apply( x$gamma[[ ind ]][ keep, ], 2, mean ) } )
  t( do.call( 'cbind', res ) )[ , -1 ]
  
}

col_fun <- function( x ) {
  
  ifelse( colnames( x ) == '', 'orange', 'dodgerblue' )
  
}

identify_par <- function( x ) {
  
  x <- abs( x )
  y <- character( length( x ) )
  y[ which( x == 0 ) ] <- ''
  y[ which( x == 0.5 ) ] <- '0.5'
  y[ which( x == 0.8 ) ] <- '0.8'
  y[ which( x == 1 ) ] <- '1'
  y
  
}

### Run simulations
set.seed( 7 )
reps <- 1:6
num_cores <- length( reps )
its <- 50000
bi <- ( its / 2 ) - 1
privar <- 4

# Large
l_sim <- as.data.frame( matrix( c( 0, -0.8, 0.5, 0, 0, 0, 0,
                                   0, 0, 0, 0.8, -0.5, 0, 0 ), ncol = 2 ) )
sim_function( reps = reps, folder = 'Large', N = 1000, beta_sim = l_sim, 
              n_cores = num_cores, prior_variance = privar, MCMC_iterations = its, 
              burn_in = bi ) # 40667.088 sec elapsed

# Mid
m_sim <- as.data.frame( matrix( c( 0, -1, 1, -0.8, 0.8, 0, 0, rep( 0, 10 ),
                                   0, 0, 0, 0.8, -0.8, 1, -1, rep( 0, 10 ) ), ncol = 2 ) )
sim_function( reps = reps, folder = 'Mid', N = 500, beta_sim = m_sim, 
              n_cores = num_cores, prior_variance = privar, MCMC_iterations = its, 
              burn_in = bi ) # 22257.665 sec elapsed

# Small
s_sim <- as.data.frame( matrix( c( 0, -1, 1, -0.8, 0.8, 0, 0, rep( 0, 17 ),
                                   0, 0, 0, 0.8, -0.8, 1, -1, rep( 0, 17 ) ), ncol = 2 ) )
sim_function( reps = reps, folder = 'Small', N = 94, beta_sim = s_sim, 
              n_cores = num_cores, prior_variance = privar, MCMC_iterations = its, 
              burn_in = bi ) # 6739.451 sec elapsed

