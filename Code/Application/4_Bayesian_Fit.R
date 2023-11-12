##### Analysis of AD data using polychotomous H&H (2006) model with variable selection 
##### Daniel Dempsey

### Load code and libraries
source( 'Code/Main_Software/Polychotomous_VS_Implementation.R' )
source( 'Code/Application/3_Output_Function.R' )
dyn.load( 'Code/Main_Software/KS_RNG.so' )
library( readr ) # For reading in csv files
library( dplyr ) # For data manipulation tools
library( parallel ) # For parallelising over runs
library( tictoc ) # For timing the code
library( corrplot ) # For plotting correlation matrices
library( ggplot2 ) # Extra plotting functionality
library( extrafont ) # So the Arial font is available 
library( tidyr ) # Provides extra data frame manipulation tools
RNGkind("L'Ecuyer-CMRG") # Sets a RNG stream that carries across parallel computing, so that code is reproducible

### Colour pallette for plots, mostly taken from: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
mycol <- c( "#FFBF00", "#56B4E9", "#009E73", "#F0E442",  
            "#0072B2", "#D55E00", "#CC79A7", "#999999", "#FFFFFF" )
scales::show_col( mycol )

### Read in data
complete_dat <- read_csv( 'Data/analysis_dat.csv' )
drop_cols <- c( 'D8_Memory Effector', 'D8_CD4+', 'CLA- Effector', 'D0_CD3', # Removed as excluded from analysis
                'TS_IL-4', 'TS_IL-17A', # Removed as they are missing for a lot of patients
                'D0_Effector T cells', 'CLA+ Effector' )  # Removed as they are highly correlated with another variable in the same group
#TS_inds <- colnames(complete_dat)[ which( substr(colnames(complete_dat), 1, 2) == 'TS' ) ]
#drop_cols <- unique( c( drop_cols, TS_inds ) )
final_dat <- select( complete_dat, -all_of(drop_cols) )
dim( final_dat ) # 93 x 59

### Create a function that runs model and creates output from result
setwd( 'Output/AD_Analysis_Output/' )
ad_analysis <- function( dat, cohort_keep, imsup_keep = FALSE, n_runs = 6, n_cores = 6, mits = 100000, ... ) {
  
  ## Set specified directory
  wd <- paste( cohort_keep, collapse = '_' )
  if ( imsup_keep ) { wd <- paste0( 'Imsup_', wd ) }
  original_wd <- getwd()
  setwd( wd )
  
  ## Select dataset
  if ( missing( cohort_keep ) ) { dat_cohort <- dat }
  else {
    col_inds <- c( 1:2, which( substr( colnames( dat ), 1, 2 ) %in% cohort_keep ) )
    dat_cohort <- select( dat, all_of( col_inds ) )
  }
  
  ## Clean column names
  all_names <- colnames( dat_cohort )
  #change_name_inds <- which( ( substr( all_names, 1, 2 ) %in% c('D0', 'D8') ) & !( substr( all_names, 4, 13 ) == 'Memory CD4' ))
  change_name_inds <- which( substr( all_names, 1, 2 ) %in% c('D0', 'D8') )
  colnames( dat_cohort )[change_name_inds] <- substr( all_names[change_name_inds], 4, nchar(all_names[change_name_inds]) )
  
  ## Remove missing values
  miss_inds <- apply( dat_cohort, 1, function( x ) { any( is.na( x ) ) } )
  dat_final <- filter( dat_cohort, !miss_inds )
  
  X <- select( dat_final, -all_of(c('ID', 'Group')) )
  y <- dat_final$Group
  
  ## Correlation matrix of covariates
  cor_mat <- cor( X )
  pdf( 'Correlation_Matrix.pdf', height = 12, width = 12 )
  corrplot( cor_mat )
  dev.off()
  
  # Violin Plots
  group_labs <- c('AD S.aureus', 'AD Control', 'H Control')
  group_cols <- mycol[ c(1, 3, 2) ]
  use_cols <- group_cols[ match( y, group_labs ) ]
  
  dir.create( 'Violin' )
  
  #pdf( 'Violin.pdf' )
  for ( i in 1:ncol( X )  ) {
    
    nm <- colnames( X )[i]
    vplot <- ggplot(X, aes(x = factor(y, levels = group_labs), 
                           y = X[[i]], col = use_cols)) +
      geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
      geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
      theme(legend.position = "none", text = element_text('TT Arial')) +
      ggtitle( nm ) + ylab( '' ) + xlab( '' )
    ggsave( vplot, filename = nm, device = cairo_pdf, path = 'Violin' )
    
  }
  #dev.off()
  
  ## Create list of starting values for the MCMC algorithm
  p <- ncol( X )
  gamma_start_list <- vector( 'list', n_runs )
  gamma_start_list[[ 1 ]] <- c( TRUE, rep( FALSE, p ) ) # All FALSE
  gamma_start_list[[ 2 ]] <- rep( TRUE, p + 1 ) # All TRUE
  
  for ( i in 3:n_runs ) { # Random starts 
    gamma_start_list[[ i ]] <- c( TRUE, as.logical( rbinom( p, 1, 0.5 ) ) )
  }
  
  beta_start_list <- vector( 'list', n_runs )
  beta_start_list[[ 1 ]] <- 0 # All zero
  beta_start_list[[ 2 ]] <- rep( 0, p + 1 ) # All zero
  
  for ( i in 3:n_runs ) { # Random starts
    beta_start_list[[ i ]] <- rnorm( sum( gamma_start_list[[ i ]] ), sd = 3 )
  }
  
  gamma_prior <- 0.5
  
  ### Add immunology data if requested by user
  if ( imsup_keep ) {
    
    im_dat <- select( dat, ID, Immunosuppressed ) %>% filter( ID %in% final_dat$ID )
    dat_final <- merge( dat_final, im_dat, by = 'ID' )
    X <- select( dat_final, -all_of(c('ID', 'Group')) )
    y <- dat_final$Group
    
    ## Violin plots with immunosuppression
    highlights <- ifelse( X$Immunosuppressed == 1, 1, 0.1 )
    
    dir.create( 'Violin_Immunosuppression' )
    
    #pdf( 'Violin_Immunosuppression.pdf' )
    for ( i in 1:(ncol( X )-1)  ) {
      
      nm <- colnames( X )[i]
      vplot <- ggplot(X, aes(x = factor(y, levels = group_labs), 
                             y = X[[i]], col = use_cols)) +
        geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
        geom_jitter(position = position_jitter(seed = 1, width = 0.2), alpha = highlights) +
        theme(legend.position = "none", text = element_text('TT Arial')) +
        ggtitle( nm ) + ylab( '' ) + xlab( '' )
      ggsave( vplot, filename = nm, device = cairo_pdf, path = 'Violin_Immunosuppression' )
      
    }
    #dev.off()
    
    gamma_start_list <- lapply( gamma_start_list, c, TRUE )
    beta_start_list <- lapply( beta_start_list, c, 0 )
    
    gamma_prior <- c( rep( gamma_prior, p ), 1 ) 
    
  }
  
  ## Save dataset
  write_csv( dat_final, file = 'dat.csv' )
  
  #### Run MCMC inference
  tic( )
  ad_fit <- mcMap( bayes_logistic, gamma_start = gamma_start_list, 
                   beta_start = beta_start_list, 
                   MoreArgs = list( y = y, X = X, MCMC_iterations = mits,
                                    prior_variance = 4, indicator_prior = gamma_prior,
                                    baseline = 'AD Control', ... ),
                   mc.cores = n_cores )
  toc( )
  
  output_fun( ad_fit )
  
  setwd( original_wd )
  
  return( NULL )
  
}

### Fit data
# following code produces warnings; this is just tidyverse tibble quibbles
set.seed( 93 )
ad_analysis( final_dat, c('D0', 'D8') ) # 457.564 sec elapsed

set.seed( 69 )
ad_analysis( final_dat, c('D0', 'CL', 'TS') ) # 523.913 sec elapsed

set.seed( 931 )
ad_analysis( final_dat, c('D0', 'D8'), imsup_keep = TRUE ) # 519.207 sec elapsed

set.seed( 691 )
ad_analysis( final_dat, c('D0', 'CL', 'TS'), imsup_keep = TRUE ) # 577.642 sec elapsed

