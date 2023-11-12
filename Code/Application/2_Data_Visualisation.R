##### Visualise Data
##### Daniel Dempsey

### Read in library
library( readr ) # Reading in csv files
library( dplyr ) # Data frame manipulation tools
library( corrplot ) # For plotting correlation matrices
library( ggplot2 ) # Visualisation tools
library( extrafont ) # So the Arial font is available 

### Colour pallette for plots, mostly taken from: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
mycol <- c( "#FFBF00", "#56B4E9", "#009E73", "#F0E442",  
            "#0072B2", "#D55E00", "#CC79A7", "#999999", "#FFFFFF" )
scales::show_col( mycol )

### Read in Data
complete_dat <- read_csv( 'Data/analysis_dat.csv' )
complete_dat$Group <- factor( complete_dat$Group, levels = c( 'AD S.aureus', 'AD Control', 'H Control' ) )

### Set colours
group_cols <- mycol[ c(1, 3, 2) ]
use_cols <- group_cols[ match( complete_dat$Group, c( 'H Control', 'AD S.aureus', 'AD Control' ) ) ]

### Clean column names
all_names <- colnames( complete_dat )
#change_name_inds <- which( ( substr( all_names, 1, 2 ) %in% c('D0', 'D8') ) & !( substr( all_names, 4, 13 ) == 'Memory CD4' ) )
change_name_inds <- which( substr( all_names, 1, 2 ) %in% c('D0', 'D8') )
colnames(complete_dat)[change_name_inds] <- substr( all_names[change_name_inds], 4, nchar(all_names)[change_name_inds] )

### Barplot of cohort sizes
cohort_sizes <- table( complete_dat$Group ) %>% as.data.frame

#pdf( 'Output/Data_Visualisation/Response_Group_Sizes.pdf' )
response_bar <- ggplot( cohort_sizes, aes(x = Var1, y = Freq) ) + geom_bar( stat = 'identity', fill = group_cols ) + 
  ylab( '' ) + xlab( '' ) + ggtitle( 'Patient Group Sizes' ) + theme( axis.text.x = element_text(face = 'italic'),
                                                                      text = element_text(family = 'TT Arial') )
ggsave( response_bar, filename = 'Output/Data_Visualisation/Response_Group_Sizes.pdf', device = cairo_pdf )
#dev.off()

### Violin plots of each variate for each group
num_dat <- select( complete_dat, -'Immunosuppressed' )

dir.create( 'Output/Data_Visualisation/Violin' )

#pdf( 'Output/Data_Visualisation/Violin.pdf' )
for ( i in 3:ncol( num_dat )  ) {
  
  nm <- colnames( num_dat )[i]
  vplot <- ggplot(num_dat, aes(x = num_dat$Group, y = num_dat[[i]], col = use_cols)) +
    geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
    geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
    theme(legend.position = "none", axis.text.x = element_text(face = 'italic'),
          text = element_text(family = 'TT Arial')) +
    ggtitle( nm ) + ylab( '' ) + xlab( '' )
  ggsave( vplot, filename = nm, device = cairo_pdf, path = 'Output/Data_Visualisation/Violin' )
  
}
#dev.off()

### Correlation matrix of numerical covariates
cor_mat <- cor( select( num_dat, -c( 'ID', 'Group' ) ), use = 'complete.obs' )

pdf( 'Output/Data_Visualisation/Correlation_Matrix.pdf', height = 12, width = 12 )
corrplot( cor_mat, tl.cex = 0.8 )
dev.off()

# Highlight correlations larger than 0.9
inds <- which( abs(cor_mat) > 0.95, arr.ind = TRUE )
inds[ which( inds[, 1] != inds[, 2] ), ]

### Immunosuppression
table( complete_dat$Group, complete_dat$Immunosuppressed )

### Violin plot again, but immunosuppresants highlighted
highlights <- ifelse( complete_dat$Immunosuppressed == 1, 1, 0.1 )

dir.create( 'Output/Data_Visualisation/Violin_Immunosuppression' )

#pdf( 'Output/Data_Visualisation/Violin_Immunosuppression.pdf' )
for ( i in 3:ncol( num_dat )  ) {
  
  nm <- colnames( num_dat )[i]
  vplot <- ggplot(num_dat, aes(x = num_dat$Group, y = num_dat[[i]], col = use_cols)) +
    geom_violin(alpha = 0.3, size = 0.1, fill = 'grey', col = 'grey') +
    geom_jitter(position = position_jitter(seed = 1, width = 0.2), alpha = highlights) +
    theme(legend.position = "none", axis.text.x = element_text(face = 'italic'),
          text = element_text(family = 'TT Arial')) +
    ggtitle( nm ) + ylab( '' ) + xlab( '' )
  ggsave( vplot, filename = nm, device = cairo_pdf, path = 'Output/Data_Visualisation/Violin_Immunosuppression' )
  
}
#dev.off()

