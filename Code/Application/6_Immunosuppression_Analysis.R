##### Analysis specifically of Immunosuppression data 
##### Daniel Dempsey

library( readr )
library( dplyr )
library( ggplot2 )
library( bayesQR )

### Colour pallette for plots, taken from: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
mycol <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442",  
            "#0072B2", "#D55E00", "#CC79A7", "#999999", "#FFFFFF" )
scales::show_col( mycol )
                     
### Read in data
inactive <- read_csv( 'Data/analysis_dat.csv' ) %>% 
  filter( Group == 'Inactive' ) %>% select( -ID, -Group )

### Quantile regression
immuno_prior <- prior( Immunosuppressed ~ ., data = inactive, V0 = 10*diag(ncol(inactive)) )
immuno_fit <- bayesQR( Immunosuppressed ~ ., data = inactive, prior = immuno_prior, 
                       quantile = 0.5, ndraw = 1e5, keep = 10, seed = 43 )
sum_fit <- summary( immuno_fit, burnin = 5000 )[[1]]$betadraw
sum_fit[which( sign( sum_fit[, 2] ) == sign( sum_fit[, 3] ) ), ]

### Plot overlapping densities
inactive$Immunosuppressed <- ifelse( inactive$Immunosuppressed == 1, 'Yes', 'No' )
inactive$Immunosuppressed <- factor( inactive$Immunosuppressed )
cn <- colnames( inactive )

pdf( 'Output/Extra_Visualisations/Immunosuppressed_Densities.pdf' )
for( i in 1:(ncol(inactive)-1)) {
  
  dens_plot <- ggplot(inactive, aes(x = inactive[[i]], fill = Immunosuppressed)) +
    geom_density(alpha = 0.8, color = NA) + 
    scale_fill_manual(values = mycol[1:2]) +
    ylab('') + xlab(cn[i])
  print( dens_plot )
  
}
dev.off()

