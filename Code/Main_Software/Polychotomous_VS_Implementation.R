##### Multinomial Logistic Regression with Variable Selection using Latent Variables
##### Method taken from Holmes & Held (2006) Appendices A3, A4 and A5
##### Daniel Dempsey

### Load library function for efficient sparse matrix handling
library( Matrix ) # For sparse matrix representation and functionality
library( matrixStats ) # For logSumExp function
library( ggplot2 ) # For plotting ridgeplots
library( ggridges ) # Ditto ^

### Function that updates the latent variable z (independent of lambda)
### Updates via truncated logistic distribution
z_update <- function( y, Xb, category ) {
  
  # Compute one-dimensional identifier
  yy <- ifelse( y == category, 1, 0 )
  
  # Find the quantile that truncates the distribution
  trunc <- plogis( 0, Xb )
  
  # (Inversion method) sample from a uniform distribution
  unif_lower <- ifelse( yy == 1, trunc, 0 )
  unif_upper <- ifelse( yy == 1, 1, trunc )
  u <- runif( length( yy ), unif_lower, unif_upper )
  
  # Transform to a truncated logistic variable
  qlogis( u, Xb )
  
}

### Function that updates the latent variable lambda
### Non-standard distribution; method taken from appendix of Holmes and Held (2006)
lambda_update <- function( z, Xb ) {
  
  # Initialise
  N <- length(z)
  r_all <- abs( z - Xb )
  lambda <- numeric( N )
  
  for ( i in 1:N ) {
    
    r <- r_all[i]
    ok <- FALSE
    
    # Rejection sampling
    while ( !ok ) {
      
      x <- rnorm( 1 )^2
      y <- 1 + ( x - sqrt(x*(4*r + x)) )/(2*r)
      
      u1 <- runif( 1 )
      lambda[i] <- ifelse( u1 <= 1/(1 + y), r/y, r*y )
      
      u2 <- runif( 1 )
      ok <- ifelse( lambda[i] > 4/3, rightmost(u2, lambda[i]), leftmost(u2, lambda[i]) )
      
    }
    
  }
  
  lambda
  
}

### rightmost function for lambda acceptance
rightmost <- function( u, l ) {
  
  z <- 1
  x <- exp( -l/2 )
  j <- 0
  
  while ( TRUE ) {
    
    j <- j + 1
    z <- z - ( (j + 1)^2 )*x^( (j + 1)^2 - 1 )
    if ( z > u ) { return( TRUE ) }
    
    j <- j + 1
    z <- z + ( (j + 1)^2 )*x^( (j + 1)^2 - 1 )
    if ( z < u ) { return( FALSE ) }
    
  }
  
}

### leftmost function for lambda acceptance
leftmost <- function( u, l ) {
  
  H <- ( log(2)/2 ) + ( 2.5*log(pi) ) - ( 2.5*log(l) ) - ( pi^2/(2*l) ) + ( l/2 )
  lu <- log( u )
  z <- 1
  x <- exp( -pi^2 / (2*l) )
  k <- l / pi^2
  j <- 0
  
  while ( TRUE ) {
    
    j <- j + 1
    z <- z - k*x^( j^2 - 1 )
    if ( (H + log(z)) > lu ) { return( TRUE ) }
    
    j <- j + 1
    #z <- z + ( (j + 1)^2 ) * x^( (j + 1)^2 - 1 ) # I think this is a mistake in H&H 2006
    z <- z + k*x^( j^2 - 1 )
    if ( (H + log(z)) < lu ) { return( FALSE ) }
    
  }
  
}

### Function that updates the regression parameters beta
### Gaussian Distribution
beta_update <- function( x, z, lambda, v_inv ) {
  
  # Compute necessary values
  bpars <- beta_parameters( x, z, lambda, v_inv )
  u <- rnorm( length( bpars$mean ) )
  L <- chol( bpars$variance )
  
  # Update beta
  bpars$mean + crossprod( L, u )
  
}

### Function that calculates beta's parameters
beta_parameters <- function( x, z, lambda, v_inv ) {
  
  V <- chol2inv( chol( crossprod( x, x/lambda ) + v_inv ) )
  B <- V %*% crossprod( x, z/lambda )
  list( mean = B, variance = V )
  
}

### Function that updates the covariate indicatos gamma
### Metropolis-Hastings algorithm
gamma_update <- function( gamma, x, x_full, v_inv, v_inv_full, z, 
                          lambda, gamma_prior, working_env ) {
  
  # Select component of gamma to propose update
  gamma_proposal <- gamma
  change_ind <- sample( 2:length( gamma ), 1 )
  gamma_proposal[ change_ind ] <- !gamma[ change_ind ]
  
  # Create proposal parameters
  x_proposal <- x_full[ , gamma_proposal, drop = FALSE ]
  v_inv_proposal <- v_inv_full[ gamma_proposal, gamma_proposal, drop = FALSE ]
  beta_pars <- beta_parameters( x, z, lambda, v_inv )
  beta_pars_proposal <- beta_parameters( x_proposal, z, lambda, v_inv_proposal )
  gamma_prior_eval <- sum( dbinom( gamma, 1, gamma_prior, log = TRUE ) )
  gamma_prior_eval_proposal <- sum( dbinom( gamma_proposal, 1, gamma_prior, log = TRUE ) )
  
  # Compute M-H acceptance probability
  accept <- gamma_acceptance( beta_pars_proposal, v_inv_proposal, gamma_prior_eval_proposal ) - 
    gamma_acceptance( beta_pars, v_inv, gamma_prior_eval )
  
  # Accept or reject and return result
  if ( log( runif( 1 ) ) < accept ) {
    gamma <- gamma_proposal
    working_env$x <- x_proposal
    working_env$v_inv <- v_inv_proposal
  }
  
  gamma
  
}

### Function that computes the components of the acceptance ratio
gamma_acceptance <- function( beta_pars, v_inv, gamma_prior ) {
  
  # Compute all necessary values
  B <- beta_pars$mean
  V <- beta_pars$variance
  L <- chol( V )
  V_inv <- chol2inv( L )
  l_inv <- chol( v_inv )
  detV <- sum( log( diag( L ) ) )
  detv_inv <- sum( log( diag( l_inv ) ) )
  kernel <- crossprod( B, V_inv %*% B ) / 2
  
  # Compute acceptance component
  detV + detv_inv + kernel + gamma_prior
  
}

### Function that updates the prior probability of variable inclusion (currently unused)
### Beta Distribution
p_update <- function( g ) {
  
  rbeta( 1, 1 + sum( g ), 1 + sum( 1 - g ) )
  
}

### Function that computes multinomial probabilities
pmultinom <- function( q ) {
  
  lse <- apply( q, 1, logSumExp )
  exp( q - lse )
  
}

### Function that generates 1 sample multinomial draws
rmultinom2 <- function( w ) {
  
  ind <- as.logical( rmultinom( 1, 1, w ) )
  c( 'V0', 'V1', 'V2' )[ ind ]
  
}


### Function that calculates the log-posterior
log_posterior_eval <- function( Xb_terms, y, categories, beta, 
                                gamma, prior_sd, gamma_prior ) {
  
  preds <- pmultinom( Xb_terms )
  colnames( preds ) <- categories
  ll <- log( sapply( 1:length( y ), function( x ) { preds[ x, y[ x ] ] } ) )
    
  sd <- lapply( gamma, function( x ) { rep( prior_sd, sum( x ) ) } )
  lbp <- Map( dnorm, x = beta, sd = sd, MoreArgs = list( log = TRUE ) )
  lgp <- Map( dbinom, x = gamma, prob = gamma_prior, MoreArgs = list( size = 1, log = TRUE ) )
  
  sum( ll, sapply( lbp, sum ), sapply( lgp, sum ) )
  
}

### Function that computes Monte Carlo Posterior Predictive Distribution
ppd_eval <- function( X_new, MCMC_res, scaling = TRUE ) {
  
  if ( scaling ) { X <- cbind( 1, apply( X_new, 2, scale ) ) }
  else { X <- cbind( 1, X_new ) }
  nvars <- 1:length( MCMC_res$beta[[ 1 ]][ 1, ] )
  
  Xb_list <- lapply( MCMC_res$beta, tcrossprod, x = X )
  N <- nrow( Xb_list[[1]] )
  G <- length( Xb_list ) + 1
  K <- nrow( MCMC_res$beta[[1]] )
  Xb_mats <- rep( list( matrix( 0, ncol = G, nrow = N ) ), K )
  
  for ( k in 1:K ) {
    
    col_ext <- do.call( 'cbind', lapply( Xb_list, '[', i = 1:N, j = k ) )
    Xb_mats[[ k ]] <- cbind( col_ext, 0 )
    colnames( Xb_mats[[ k ]] )[ G ] <- MCMC_res$baseline
    
  }
  
  Reduce( '+', lapply( Xb_mats, pmultinom ) ) / K
  
}

### Full MCMC function
bayes_logistic <- function( y, X, beta_start, gamma_start, gamma_update_length = 1,
                            hyperprior = FALSE, prior_variance = 10, 
                            indicator_prior = 0.5, scaling = TRUE, baseline,
                            MCMC_iterations = 100000, burn_in, thin = 10,
                            progress_report = 500 ) {
  
  ### Initialisation
  working_env <- environment()
  P <- ncol( X ) + 1
  if ( is.matrix( X ) ) { X <- data.frame( X ) }
  if ( missing( burn_in ) ) { burn_in <- MCMC_iterations / 2 }
  
  # Find number of groups
  if ( !is.factor( y ) ) { y <- factor( y ) }
  categories <- levels( y )
  G <- length( unique( y ) ) - 1
  N <- length( y )
  
  # Set baseline
  if ( !missing( baseline ) ) {
    base_ind <- which( categories == baseline )
    categories <- c( categories[ -base_ind ], categories[ base_ind ] )
  }
  baseline <- categories[ G + 1 ]
  
  # Initialise objects where MCMC iterations will be stored
  all_beta <- lapply( 1:G, function( x ) { matrix( 0, nrow = MCMC_iterations, ncol = P ) } )
  all_gamma <- lapply( 1:G, function( x ) { matrix( F, nrow = MCMC_iterations, ncol = P ) } )
  
  nms <- c( 'Intercept', colnames( X ) )
  all_beta <- lapply( all_beta, function( x ) { colnames( x ) <- nms; x } )
  all_gamma <- lapply( all_gamma, function( x ) { colnames( x ) <- nms; x } )
  
  names( all_beta ) <- names( all_gamma ) <- categories[ 1:G ]
  
  if ( missing( gamma_start ) ) { 
    gamma_start <- lapply( 1:G, function( x ) { c( TRUE, rep( FALSE, P - 1 ) ) } ) 
  }
  else {
    gamma_start <- rep( list( gamma_start ), G )
  }
  
  current_gamma <- gamma_start
  all_gamma <- Map( x = all_gamma, y = gamma_start, function( x, y ) { x[ 1, ] <- y; x } )
  
  if ( missing( beta_start ) ) { 
    beta_start <- Map( x = 1:G, y = gamma_start, function( x, y ) { rep( 0, sum( y ) ) } )
  }
  else {
    beta_start <- rep( list( beta_start ), G )
  }
  current_beta <- beta_start
  
  # Set hyperprior if desired
  if ( hyperprior ) {
    all_p <- matrix( 0, nrow = MCMC_iterations, ncol = G )
    all_p[ 1, ] <- sapply( all_gamma, function( x ) { p_update( x[ 1, ] ) } )
    gamma_prior <- lapply( all_p[ 1, ], function( x ) { c( 1, rep( x, P - 1 ) ) } )
  }
  else {
    if ( length( indicator_prior ) == 1 ) { indicator_prior <- rep( indicator_prior, P - 1 ) }
    gamma_prior <- rep( list( c( 1, indicator_prior ) ), G )
  }
  
  # Set beta prior variance
  v_inv_full <- diag( 1 / prior_variance, P )
  prior_sd <- sqrt( prior_variance )
  
  # Scale design matrix (if desired) and add an intercept
  if ( scaling ) { x_full <- cbind( 1, apply( X, 2, scale ) ) }
  else { x_full <- cbind( 1, as.matrix( X ) ) }
  
  # Initialise beta
  all_beta <- Map( function( x, y, z ) { x[ 1, y[ 1, ] ] <- z; x }, 
                   x = all_beta, y = all_gamma, z = beta_start )
  
  # Initialise offset and latent variables
  offset_terms <- matrix( 0, nrow = N, ncol = G + 1 )
  z <- lambda <- matrix( 1, nrow = N, ncol = G )
  for ( j in 1:G ) {
    
    offset_terms[ , j ] <- x_full[ , gamma_start[[ j ]], drop = FALSE ] %*% beta_start[[ j ]]
    z[ , j ] <- z_update( y, 0, categories[ j ] )
    
  }
  
  log_posterior <- numeric( length = MCMC_iterations )
  log_posterior[ 1 ] <- log_posterior_eval( Xb = offset_terms, y = y, 
                                            categories = categories, 
                                            beta = current_beta,
                                            gamma = current_gamma, 
                                            prior_sd = prior_sd, 
                                            gamma_prior = gamma_prior )
  
  ### Main MCMC Loop
  for ( i in 2:MCMC_iterations ) {
    
    for ( j in 1:G ) {
      
      # Filter dataset and metrics by gamma
      g <- all_gamma[[ j ]][ i - 1, ]
      x <- x_full[ , g, drop = FALSE ]
      Xb <- offset_terms[ , j ]
      v_inv <- v_inv_full[ g, g, drop = FALSE ]
      
      # Compute offset
      C <- apply( offset_terms[ , -j ], 1, logSumExp )
      
      # Update gamma
      if ( i%%gamma_update_length == 0 ) {
        all_gamma[[ j ]][ i, ] <- g <- gamma_update( g, x, x_full, v_inv, v_inv_full, 
                                                     z[ , j ] + C, lambda[ , j ], 
                                                     gamma_prior[[ j ]], working_env ) 
      }
      else {
        all_gamma[[ j ]][ i, ] <- g
      }
      current_gamma[[ j ]] <- g
      
      # Update gamma hyper-parameter (if desired)
      if ( hyperprior ) {
        all_p[ i, j ] <- p_update( all_gamma[[ j ]][ i, -1 ] )
        gamma_prior[[ j ]] <- c( 1, rep( all_p[ i, j ], P - 1 ) )
      }
      
      # Update beta
      all_beta[[ j ]][ i, g ] <- beta_update( x, z[ , j ] + C, lambda[ , j ], v_inv )
      Xb <- offset_terms[ , j ] <- x %*% all_beta[[ j ]][ i, g ]
      current_beta[[ j ]] <- all_beta[[ j ]][ i, g ]
      
      # Update latent variables
      z[ , j ] <- z_update( y, Xb - C, categories[ j ] )
      #lambda[ , j ] <- lambda_update( z[ , j ] + C, Xb ) # R Implementation
      l <- lambda[ , j ]
      lambda[ , j ] <- .C( "KS_RNG", as.integer(N), as.double(z[ , j ] + C), as.double(Xb), as.double(l) )[[4]] # C Implementation
      
    }
    
    # Update log posterior
    log_posterior[ i ] <- log_posterior_eval( Xb = offset_terms, y = y, 
                                              categories = categories, 
                                              beta = current_beta, 
                                              gamma = current_gamma, 
                                              prior_sd = prior_sd, 
                                              gamma_prior = gamma_prior )
    
    # User progress report
    if ( i%%progress_report == 0 ) { cat( paste0('Iteration ', i, ' completed.\n') ) }
    
  }
  
  # Remove burn-in
  keep <- seq( burn_in, MCMC_iterations, thin )
  keep_fun <- function( x ) { x[keep, ] }
  beta_retain <- lapply( all_beta, keep_fun )
  gamma_retain <- lapply( all_gamma, keep_fun )
  lp_retain <- log_posterior[ keep ]
  
  # Return results
  res <- list( beta = beta_retain, gamma = gamma_retain, log_posterior = lp_retain, 
               baseline = baseline, data = X, response = y, MCMC_iterations = MCMC_iterations, 
               burn_in = burn_in, thin = thin )
  if ( hyperprior ) { res$p <- all_p }
  
  class( res ) <- 'bayes_multinomial'
  res
  
}

### Utility function that returns posterior sample filtered by gamma
gamma_support <- function( b, g ) {
  
  Map( function( x, y ) { x[ y ] }, x = as.data.frame( b ), y = as.data.frame( g ) )

}

### Print method for the model output
#print.bayes_multinomial <- function( x ) {
#  
#  cat( 'Bayesian Multinomial Object.\n' )
#
#}

### Summary method for the model output
summary.bayes_multinomial <- function( x ) {
  
  G <- length( x$beta )
  sel_props <- lapply( x$gamma, function( y ) { apply( y, 2, mean ) } )
  beta_list <- Map( gamma_support, b = x$beta, g = x$gamma )
  
  bm <- lapply( beta_list, function( y ) { sapply( y, mean ) } )
  bsd <- lapply( beta_list, function( y ) { sapply( y, sd ) } )
  bq1 <- lapply( beta_list, function( y ) { sapply( y, quantile, probs = 0.025 ) } )
  bq2 <- lapply( beta_list, function( y ) { sapply( y, quantile, probs = 0.975 ) } )
  
  Map( function( x, y, z1, z2, p ) { round( data.frame( mean = x, sd = y, Q2_5 = z1, 
                                                        Q97_5 = z2, prop = p ), 2 ) }, 
       bm, bsd, bq1, bq2, sel_props )
  
}

### Plot method for the model output (ridgeline)
plot.bayes_multinomial <- function( x, which_vars, which_groups, var_order, ... ) {
  
  G <- length( x$beta )
  beta_list <- Map( gamma_support, b = x$beta, g = x$gamma )
  
  if ( !missing( which_vars ) ) {
    beta_list <- lapply( beta_list, '[', which_vars )
  }
  
  if ( !missing( which_groups ) ) {
    beta_list <- beta_list[ which_groups ]
  }
  
  varlens <- lapply( beta_list, function( x ) { sapply( x, length ) } )
  varvec <- lapply( varlens, function( x ) { rep( names( x ), x ) } )
  vargroupvec <- Map( function( x, y ) { data.frame( var = x, group = y ) }, x = varvec, y = names( varvec ) )
  beta_vecs <- lapply( beta_list, do.call, what = 'c' )
  graph_dat <- do.call( 'rbind', Map( function( x, y ) { cbind( beta = x, y ) }, x = beta_vecs, y = vargroupvec ) )
  graph_dat$group <- factor( graph_dat$group, levels = rev( unique( graph_dat$group ) ) )
  
  if( !missing(var_order) ) {
    if( is.numeric( var_order ) ) { var_order <- names( beta_list[[1]] )[var_order] }
    graph_dat$var <- factor( graph_dat$var, levels = var_order )
  }
  
  ggplot( graph_dat, aes( x = beta, y = var, fill = group ), show.legend = FALSE ) + 
    geom_density_ridges( ... ) + scale_fill_discrete( name = NULL ) +  ylab( '' ) + xlab( 'Regression Coefficient' ) +
    scale_fill_manual( values = c("#56B4E9", '#E69F00') ) + geom_vline( xintercept = 0 )
  
}

