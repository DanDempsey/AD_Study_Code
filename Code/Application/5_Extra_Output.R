##### Extra output requested by the Immunologists
##### Daniel Dempsey

### Read in library
library( readxl ) # Reading in csv files
library( dplyr ) # Data frame manipulation tools
library( ggplot2 ) # For creating violin plots
library( readr ) # For reading csv filesc( 'AD S.aureus', 'AD Control', 'H Control' )
library( AMR ) # For performing G-tests
library( extrafont ) # So the Arial font is available 
#font_import()
source( 'Code/Main_Software/Polychotomous_VS_Implementation.R' )
source( 'Code/Application/3_Output_Function.R' )
dyn.load( 'Code/Main_Software/KS_RNG.so' )
#setwd( 'Output/Extra_Visualisations/' )

### Colour pallette for plots, mostly taken from: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
mycol <- c( "#FFBF00", "#56B4E9", "#009E73", "#F0E442",  
            "#0072B2", "#D55E00", "#CC79A7", "#999999", "#FFFFFF" )
scales::show_col(mycol)

### Group vector
groups <- c( 'AD S.aureus', 'AD Control', 'H Control' )

### CLA+Th17 and D0 Th17 violin plots
# Plotting function
vio_fun <- function( x, y, z, ... ) {
  
  dat <- data.frame( group = y, vals = x )
  
  # Draw violin plot
  group_cols <- mycol[ c(6, 2, 3) ]
  use_cols <- group_cols[ match( dat$group, groups ) ]
  
  vplot <- ggplot(dat, aes(x = group, y = vals, col = use_cols)) +
    geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
    geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
    theme(legend.position = "none", text = element_text('TT Arial', size = 18)) +
    ggtitle( z ) + ylab( '' ) + xlab( '' )
 
  ggsave( plot = vplot, filename = ..., device = cairo_pdf, path = "Output/Extra_Visualisations/" )
   
}

# Load in data
d8_th17 <- read_xlsx( 'Data/Daniel - Data Set - Nov 2022.xlsx', range = 'K6:K79', sheet = 3, col_names = FALSE ) %>% unlist
d0_th17 <- read_xlsx( 'Data/Daniel - Data Set - Nov 2022.xlsx', range = 'L5:L113', sheet = 1, col_names = FALSE ) %>% unlist
d8_group_raw <- read_xlsx( 'Data/Daniel - Data Set - Nov 2022.xlsx', range = 'A6:A79', sheet = 3, col_names = FALSE ) %>% unlist
d0_group_raw <- read_xlsx( 'Data/Daniel - Data Set - Nov 2022.xlsx', range = 'A5:A113', sheet = 1, col_names = FALSE ) %>% unlist

# Fix groups
d8_group <- factor( groups[ as.numeric( substr(d8_group_raw, 7, 7) ) ], levels = groups )
d0_group <- factor( groups[ as.numeric( substr(d0_group_raw, 7, 7) ) ], levels = groups )

# Remove outlier from day 0 th17
out_ind <- which.max( d0_th17 )
d0_th17_zoom <- d0_th17[-out_ind]
d0_group_zoom <- d0_group[-out_ind]

# Create plots
vio_fun( x = d8_th17, y = d8_group, z = 'CLA+ Th17', 'D8_Th17.pdf' )
vio_fun( x = d0_th17, y = d0_group, z = 'Th17', 'D0_Th17.pdf' )
vio_fun( x = d0_th17_zoom, y = d0_group_zoom, z = 'Th17 (AD Control Outlier Excluded)', 'D0_Th17_nooutlier.pdf' )

### CD8 plot
# Read in Data
complete_dat <- read_csv( 'Data/all_dat.csv' )
complete_dat$Group <- factor( complete_dat$Group, levels = groups )

# Complete data
drop_ind <- which.max( complete_dat$D0_CD8 )
xx <- complete_dat[-drop_ind, ]

group_cols <- mycol[ c(6, 2, 3) ]
use_cols <- group_cols[ match( xx$Group, groups ) ]

vplot <- ggplot(xx, aes(x = Group, y = D0_CD8, col = use_cols)) +
  geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none", text = element_text(family = 'TT Arial', size = 18)) +
  ggtitle( 'CD8 - Outlier Excluded' ) + ylab( '' ) + xlab( '' )
ggsave( plot = vplot, filename = 'Output/Extra_Visualisations/CD8_zoomed.pdf', device = cairo_pdf )

### Ridgeplots with less features
ridge_plot <- function( X, nums, file_name, ind = 1, remove_vars = NULL, xrange = c(-7, 7), 
                        use_groups = groups[-2], leg = TRUE ) {
  fit <- lapply( X, function( x ) { apply( x$gamma[[ ind ]], 2, mean ) } )
  fit_sums <- do.call( 'rbind', fit )[ , -1 ]
  fit_quantiles <- apply( fit_sums, 2, quantile, probs = 0.5 )
  use_cols <- names( fit_quantiles )[order( fit_quantiles, decreasing = TRUE )[nums]]
  im_ind <- which( use_cols %in% remove_vars )
  if ( length(im_ind) > 0 ) {
    use_cols <- use_cols[-im_ind]
  }
  rplot <- plot( X[[1]], use_cols, var_order = rev(use_cols), which_groups = use_groups ) + 
    xlim( xrange ) + theme( text = element_text(family = 'TT Arial', size = 18), legend.title = element_blank() )
  if( !leg ) {
    rplot <- rplot + theme( legend.position = 'none' ) 
  }
  ggsave( plot = rplot, filename = file_name, device = cairo_pdf )
}

load( 'Output/AD_Analysis_Output/D0_D8/ad_fit.R' )
ridge_plot( X, 1:7, 'Output/Extra_Visualisations/Small_Ridge/D0_D8_nonImsup.pdf' )
ridge_plot( X, 1:7, use_groups = groups[1], leg = F, 'Output/Extra_Visualisations/Small_Ridge/D0_D8_nonImsup_ActiveOnly.pdf' )
load( 'Output/AD_Analysis_Output/D0_CL_TS/ad_fit.R' )
ridge_plot( X, 1:11, 'Output/Extra_Visualisations/Small_Ridge/D0_CL_TS_nonImsup.pdf' )
ridge_plot( X, 1:11, use_groups = groups[1], leg = F, 'Output/Extra_Visualisations/Small_Ridge/D0_CL_TS_nonImsup_ActiveOnly.pdf' )
load( 'Output/AD_Analysis_Output/Imsup_D0_D8/ad_fit.R' )
ridge_plot( X, 1:10, remove_vars = 'Immunosuppressed', 'Output/Extra_Visualisations/Small_Ridge/D0_D8_Imsup.pdf' )
ridge_plot( X, 1:10, use_groups = groups[1], leg = F, remove_vars = 'Immunosuppressed', 
            'Output/Extra_Visualisations/Small_Ridge/D0_D8_Imsup_ActiveOnly.pdf' )
load( 'Output/AD_Analysis_Output/Imsup_D0_CL_TS/ad_fit.R' )
ridge_plot( X, 1:9, remove_vars = 'Immunosuppressed', 'Output/Extra_Visualisations/Small_Ridge/D0_CL_TS_Imsup.pdf' )
ridge_plot( X, 1:9, use_groups = groups[1], leg = F, remove_vars = 'Immunosuppressed', 
            'Output/Extra_Visualisations/Small_Ridge/D0_CL_TS_Imsup_ActiveOnly.pdf' )

### Table 4 p-values
HC <- read_excel( 'Data/Swab table results - for Daniel 6 11 23.xlsx', sheet = 'Healthy Controls', 
                  range = 'A6:B16', col_names = FALSE ) %>% na.omit
HC[[3]] <- 35 - HC[[2]]

AD <- read_excel( 'Data/Swab table results - for Daniel 6 11 23.xlsx', sheet = 'AD Control', 
                  range = 'A4:B21', col_names = FALSE ) %>% na.omit
AD[[3]] <- 46 - AD[[2]]

SA <- read_excel( 'Data/Swab table results - for Daniel 6 11 23.xlsx', sheet = 'AD S.aureus', 
                  range = 'A4:B19', col_names = FALSE ) %>% na.omit
SA[[3]] <- 12 - SA[[2]]

tab_fun <- function( ... ) {
  do.call( 'rbind', list(...) )
}

p_vals <- numeric( nrow(SA) )
for( i in 1:4 ) {
  xtab <- tab_fun( AD[i, -1], SA[i, -1] )
  p_vals[i] <- fisher.test( xtab )$p.value
}

for( i in 5:nrow(SA) ) {
  xtab <- tab_fun( HC[i-4, -1], AD[i, -1], SA[i, -1] )
  p_vals[i] <- g.test( xtab )$p.value
}

sevs <- HC[1:4, 1] %>% unlist
groups <- c( 'Lesional skin', 'Non lesional skin', 'Nose' )
cc_vec <- paste0( rep(groups, each = 4), ' : ', sevs )
data.frame( 'Colony Counts' = cc_vec, 'p-vals' = p_vals ) %>% write_csv( 'Output/Misc/Tab4.csv' )

