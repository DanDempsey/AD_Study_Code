##### Function that creates Output from MCMC fits
##### Daniel Dempsey

output_fun <- function( X ) {
  
  covar_names <- colnames( X[[1]]$beta[[1]] )
  group_names <- names( X[[1]]$beta )
  mycol <- c( "#FFBF00", "#56B4E9", "#009E73", "#F0E442",  
              "#0072B2", "#D55E00", "#CC79A7", "#999999", "#FFFFFF" )
  
  ### Variable Selection, Gamma
  #for ( i in 1:length( X ) ) {
  #  
  #  for ( j in 1:length( group_names ) ) {
  #    
  #    gamma_fit <- as.matrix( X[[i]]$gamma[[j]] )
  #    sel_tab <- apply( gamma_fit, 2, mean )[-1]
  #    sel_tab_frame <- data.frame( Var = names( sel_tab ), Num = sel_tab )
  #    sel_tab_ord <- order( sel_tab )
  #    sel_tab_frame$Var <- factor( sel_tab_frame$Var, levels = sel_tab_frame$Var[sel_tab_ord] )
  #    
  #    # Barplot
  #    var_ind_bar <- ggplot( sel_tab_frame, aes(x = Num, y = Var) ) + geom_bar( stat = 'identity', fill = mycol[2] ) + 
  #      xlim( 0, 1 ) + xlab( '' ) + ylab( '' ) + ggtitle( paste0( group_names[j], ', Run ', i) ) + 
  #      theme( text = element_text(size = 18, family = 'TT Arial') )
  #    
  #    #pdf( paste0( 'Gamma_Barplots/Run_', i, '_', group_names[j], '.pdf' ), width = 12, height = 12 )
  #    ggsave( var_ind_bar, filename = paste0( 'Gamma_Barplots/Run_', i, '_', group_names[j], '.pdf' ), device = cairo_pdf )
  #    #dev.off()
  #  
  #  }
  #  
  #}
  
  ### Individual boxplots of gamma for each group
  gamma_boxplot <- function( ind, filename ) {
    
    fit <- lapply( X, function( x ) { apply( x$gamma[[ ind ]], 2, mean ) } )
    fit_sums <- do.call( 'rbind', fit )[ , -1 ] %>% as.data.frame
    fit_sums_long <- pivot_longer( fit_sums, cols = colnames(fit_sums) )
    
    fit_quantiles <- sapply( fit_sums, quantile, probs = 0.5 )
    ords <- order( fit_quantiles )
    fit_sums_long$name <- factor( fit_sums_long$name, levels = fit_sums_long$name[ords] )
    
    bplot <- ggplot( fit_sums_long, aes( x = value, y = name ) ) + geom_vline( xintercept = 0.5, linetype = 'dashed' ) + 
      geom_boxplot( fill = mycol[2] ) + xlab( '' ) + ylab( '' ) + xlim( 0, 1 ) + 
      theme( text = element_text(size = 18, family = 'TT Arial') ) +
      ggtitle( paste0( 'All Runs Selection Proportion of the ', names( X[[1]]$gamma )[ ind ], ' Group' ) )
    
    ggsave( bplot, filename = filename, device = cairo_pdf, height = 12, width = 12, units = 'in' )
    ords
    
  }
  
  #pdf( 'All_Boxplot_Active.pdf', height = 12, width = 12 )
  g1_ords <- gamma_boxplot( 1, 'All_Boxplot_Active.pdf' ) 
  #dev.off()
  
  #pdf( 'All_Boxplot_Healthy.pdf', height = 12, width = 12 )
  g2_ords <- gamma_boxplot( 2, 'All_Boxplot_Healthy.pdf' ) 
  #dev.off()
  
  ### Both Groups Together
  dat_construct <- function( i ) {
  
    fit <- lapply( X, function( x ) { apply( x$gamma[[i]], 2, mean ) } )
    fit_sums <- do.call( 'rbind', fit )[ , -1 ] %>% as.data.frame
    fit_sums_long <- pivot_longer( fit_sums, cols = colnames(fit_sums) )
    
    fit_quantiles <- sapply( fit_sums, quantile, probs = 0.5 )
    ords <- order( fit_quantiles )
    fit_sums_long$name <- factor( fit_sums_long$name, levels = fit_sums_long$name[ords] )
    
    fit_sums_long$group <- factor( names( X[[1]]$gamma )[i], levels = names( X[[1]]$gamma ) )
    fit_sums_long
    
  }
  
  group_dats <- do.call( 'rbind', lapply( 1:2, dat_construct ) )
  
  bplot <- ggplot( group_dats, aes( x = value, y = name, fill = group ) ) + 
    geom_vline( xintercept = 0.5, linetype = 'dashed' ) + 
    geom_boxplot( ) + xlab( '' ) + ylab( '' ) + xlim( 0, 1 ) + ggtitle( 'All Runs Selection Proportion' ) + 
    scale_fill_manual( name = '', values = mycol[1:2] ) + 
    theme( legend.position = 'bottom', text = element_text(size = 18, family = 'TT Arial') )
  
  #pdf( 'All_Selection_Boxplots.pdf', height = 12, width = 12 )
  ggsave( bplot, filename = 'All_Selection_Boxplots.pdf', device = cairo_pdf, height = 12, width = 12, units = 'in' )
  #dev.off()
  
  ### Slope Parameters, Beta
  #density_plot <- function( i ) {
  #  
  #  vars <- colnames( X[[1]]$beta[[1]] )
  #  beta_draws_1 <- lapply( X, function( x ) { x$beta[[i]] %>% as.data.frame } )
  #  gamma_draws <- do.call( 'rbind', lapply( X, function( x ) { x$gamma[[i]] %>% as.data.frame } ) ) 
  #  beta_draws <-do.call( 'rbind', Map( function(x, y) { x$Run <- y; x }, x = beta_draws_1, y = 1:6 ) )
  #  
  #  beta_long <- pivot_longer( beta_draws, cols = colnames(beta_draws)[-ncol(beta_draws)] )
  #  gamma_long <- pivot_longer( gamma_draws, cols = colnames( gamma_draws ) )
  #  
  #  beta_long_retain <- filter( beta_long, gamma_long$value )
  #  beta_long_retain$Run <- as.factor(beta_long_retain$Run)
  #  
  #  pdf( paste0('Beta_Density_Plot', '_', group_names[i], '.pdf') )
  #  for ( j in 1:length(vars) ) {
  #    
  #    var_dat <- filter( beta_long_retain, name == vars[j] )
  #    dplot <- ggplot( var_dat, aes(x = value, colour = Run) ) + geom_density() + 
  #      ylab( 'Density' ) + xlab( vars[j] ) + 
  #      ggtitle( paste0(vars[j], ' for the ', group_names[i], ' Group') ) +
  #      theme( text = element_text(size = 18, family = 'TT Arial') )
  #    
  #    print( dplot )  
  #  
  #  }
  #  dev.off()
  #  
  #}
  
  #density_plot( 1 )
  #density_plot( 2 )
  
  ### Beta Ridgeplots
  for ( i in 1:length( X ) ) {
    
    rplot <- plot( X[[ i ]], var_order = c(0, g1_ords) + 1 ) + coord_cartesian( xlim = c( -7, 7 ) ) + 
      theme( text = element_text(size = 18, family = 'TT Arial') )
    
    ggsave( rplot, filename = paste0('Beta_Ridgeplots/All_Run_', i, '.pdf'), 
            device = cairo_pdf, height = 12, width = 12, units = 'in' )
    
    for ( j in 1:2 ) {
      
      rplot <- plot( X[[ i ]], which_groups = j, var_order = c(0, g1_ords) + 1 ) + 
        theme( legend.position = "none", text = element_text(size = 18, family = 'TT Arial') ) + 
        coord_cartesian( xlim = c( -7, 7 ) ) +
        ggtitle( paste0( group_names[j], ' Group' ) )
      ggsave( rplot, filename = paste0('Beta_Ridgeplots/Run_', i, '_Group_', j, '.pdf'), device = cairo_pdf, 
              height = 12, width = 12, units = 'in' )
      
    }
    
  }
  
  ### Log Posterior Traceplot
  for ( i in 1:length( X ) ) {
    
    pdf( paste0( 'Log_Posterior_Traceplots/Run_', i, '.pdf' ) )
    plot( X[[i]]$log_posterior, type = 'l', main = paste0( 'Run ', i ), ylab = 'Log Posterior' )
    dev.off()
    
  }
  
  ### PPD
  ppd <- lapply( X, ppd_eval, X_new = X[[1]]$data )
  for ( i in 1:length( ppd ) ) {
    write_csv( as.data.frame( ppd[[ i ]] ), paste0( 'PPD/ppd_Run_', i, '.csv' ) )
  }
  
  actual <- X[[1]]$response
  
  for ( i in 1:length( ppd ) ) {
    
    categories <- colnames( ppd[[ i ]] )
    preds <- categories[ apply( ppd[[ i ]], 1, which.max ) ]
    
    write.csv( table( preds, actual ), paste0( 'PPD/ct_Run_', i, '.csv' ) )
    
  }
  
  save( X, file = 'ad_fit.R' )
  
}

