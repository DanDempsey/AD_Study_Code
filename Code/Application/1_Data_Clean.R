##### Read in and clean data
##### Daniel Dempsey

### Load in libraries and set working directory
library( readxl ) # Reading in excel files
library( readr ) # Writing csv files
library( dplyr ) # Data frame manipulation tools
library( tidyr ) # More data frame manipulation tools
library( visdat ) # Tools for inspecting missing values
setwd( 'Data' )

### Read in day 0 data
xl_name <- 'Daniel - Data Set - Nov 2022.xlsx'
day0_raw <- read_xlsx( xl_name, range = 'A4:N113', sheet = 1 )
colnames( day0_raw ) <- c( 'Group', 'ID', paste0('D0_', colnames(day0_raw)[-c(1:2)]) )

### Read in day 8 non-CLA data
day8nc_raw <- read_xlsx( xl_name, range = 'A5:O101', sheet = 2, col_names = FALSE )
day8nc_header <- read_xlsx( xl_name, range = 'C3:O3', sheet = 2, col_names = FALSE )
day8nc_header2 <- c( 'Group', 'ID', paste0('D8_', unlist( day8nc_header )) )
colnames( day8nc_raw ) <- day8nc_header2

### Read in day 8 CLA data
day8c_raw <- read_xlsx( xl_name, range = 'A6:AA79', sheet = 3, col_names = FALSE )
day8c_header <- read_xlsx( xl_name, range = 'C4:AA4', sheet = 3, col_names = FALSE )
day8c_header2 <- c( 'Group', 'ID', unlist( day8c_header ) )
colnames( day8c_raw ) <- day8c_header2 

### Check that CD4 is the same in day8 datasets
common_ids <- intersect( day8c_raw$ID, day8nc_raw$ID ) 
all( filter( day8nc_raw, ID %in% common_ids )$`D8_CD4+ - Not Corrected` == 
       filter( day8c_raw, ID %in% common_ids )$`CD4+ - Uncorrected` )

# TRUE; remove that column from CLA data
day8c_raw2 <- select( day8c_raw, -'CD4+ - Uncorrected' )

### Fix the ID so that the repeat sample marker is removed
ID_fix <- function( x ) { x$ID <- substr(x$ID, 1, 6); x }
day0_idfix <- ID_fix( day0_raw )
day8nc_idfix <- ID_fix( day8nc_raw )
day8c_idfix <- ID_fix( day8c_raw2 )

# Merge these data
day08_all <- merge( day0_idfix, day8nc_idfix, by = c( 'ID', 'Group'), all = TRUE ) %>%
  merge( y = day8c_idfix, by = c( 'ID', 'Group'), all = TRUE )

# Redo Group assignments
gids <- day08_all$Group
groups <- c( Active = 'AD S.aureus', Inactive = 'AD Control', Healthy = 'H Control' )
day08_all$Group <- groups[ as.numeric( substr(gids, 7, 7) ) ]

### Read in tape strip data
tapestrip_filename <- 'Tape strip tables for Daniel including IL4 & IL17 Oct 2022.xlsx'

measure_header <- read_xlsx( tapestrip_filename, range = 'C7:BH7', sheet = 'Results', col_names = FALSE ) %>% 
  unlist

needed_columns_dat <- which( measure_header %in% c( 'Total NMF', 'pg/µg protein' ) )
needed_columns_name <- which( measure_header %in% c( 'Total NMF', 'pg/mL' ) )

ts_header_raw <- read_xlsx( tapestrip_filename, range = 'C6:BF6', sheet = 'Results', col_names = FALSE ) %>%
  select( all_of( needed_columns_name ) ) %>% unlist
ts_header <- paste0( 'TS_', c( 'NMF', ts_header_raw[ 2:length( ts_header_raw ) ] ) )

ts_data <- read_xlsx( tapestrip_filename, range = 'C8:BH115', sheet = 'Results', 
                      col_names = FALSE ) %>% select( all_of( needed_columns_dat ) )
colnames( ts_data ) <- ts_header

ts_cohort_data <- read_xlsx( tapestrip_filename, range = 'A8:B115', sheet = 'Results', 
                             col_names = FALSE )
colnames( ts_cohort_data ) <- c( 'ID', 'Group' ) 
ts_cohort_data$ID <- gsub( '-', '', ts_cohort_data$ID )
ts_cohort_data$Group <- groups[ ts_cohort_data$Group ]

tapestrip <- cbind( ts_cohort_data, ts_data )

# Combine Day 0, Day 8 and Tapestrip data
all_dat <- merge( day08_all, tapestrip, by = c( 'ID', 'Group' ), all = TRUE )
nrow( all_dat ) == length( unique(all_dat$ID) ) # Check to make sure there are no conflicts

# Fix the column names
gsub_alt <- function( x, pattern, replacement ) {
  gsub( pattern, replacement, x )
}

colnames( all_dat ) <- colnames( all_dat ) %>% gsub_alt( 'Mem ', 'Memory ' ) %>%
  gsub_alt( 'Mem\\+', 'Memory\\+' ) %>% gsub_alt( 'Tregs', 'Treg' ) %>%
  gsub_alt( 'IL4\\+IL13\\+Th2', 'IL4\\+IL13\\+' ) %>% gsub_alt( 'IL4\\+ Th2', 'IL4\\+') %>%
  gsub_alt( 'IL13\\+ Th2', 'IL13\\+' )

# Change D0 Memory CD4 specifically
d0memcd4_ind <- which( colnames( all_dat ) == 'D0_Memory CD4' )
colnames( all_dat )[d0memcd4_ind] <- 'D0_CD4+CD45RO+'

# Write this file
write_csv( all_dat, 'all_dat.csv' )

# Delete the rows where Day 0 and Day 8 is missing
TS_CL_inds <- which( substr(colnames(all_dat), 1, 2) %in% c('TS', 'CL') )
D0_D8_miss <- apply( select( all_dat, -all_of(TS_CL_inds) ), 1, function( x ) { any( is.na( x ) ) } )
complete_dat <- filter( all_dat, !D0_D8_miss )
dim( complete_dat )

# Make sure patients missing non-CLA data also missing CLA
D8_inds <- which( substr(colnames(all_dat), 1, 2) == 'D8' )
CL_inds <- which( substr(colnames(all_dat), 1, 2) == 'CL' )

D8_miss <- all_dat[apply( select( all_dat, all_of(D8_inds) ), 1, function( x ) { any( is.na( x ) ) } ), ]$ID
CL_miss <- all_dat[apply( select( all_dat, all_of(CL_inds) ), 1, function( x ) { any( is.na( x ) ) } ), ]$ID

all( D8_miss %in% CL_miss ) # okay

# Inspect the remaining missing values
pdf( '../Output/Data_Visualisation/Missing_Values_Raw_Dataset.pdf', width = 12 )
vis_miss( complete_dat )
dev.off()

#any_miss <- apply( complete_dat, 1, function( x ) { any( is.na( x ) ) } )
#table( complete_dat$Group, any_miss ) # 3 missing in Active, 15 in Inactive, 5 in Healthy
#complete_dat$ID[ any_miss ] # Patients that are missing data 

# Take dataset with only complete observations
#complete_dat <- filter( all_dat, !any_miss )
#nrow( all_dat ) # 118 rows total
#nrow( complete_dat ) # 71 complete rows 

# No value should be negative
neg_inds <- which( complete_dat < 0, arr.ind = TRUE ) # four negative values
colnames( complete_dat )[neg_inds[, 2]]
complete_dat[ neg_inds ] <- 0 # change to zero

# Remove the following patient as they have Netherdon's disease; not truly part of study group
final_dat <- filter( complete_dat, !(ID == 'VAD093') )

# Change some column names
nc <- grep( 'not corrected', tolower(colnames(final_dat)) ) 
nbc <- grep( 'not background corrected', tolower(colnames(final_dat)) )
num_char <- nchar( colnames(final_dat) )
colnames(final_dat)[nc] <- substr( colnames(final_dat)[nc], 1, num_char[nc] - 16 )
colnames(final_dat)[nbc] <- substr( colnames(final_dat)[nbc], 1, num_char[nbc] - 27 )

### Read in immuno-suppression Data
imsup <- read_xlsx( 'Immunosuppressed patients - 93.xlsx' ) %>% select( all_of(1:3) )
imsup$Immunosuppressed <- ifelse( imsup$Immunosuppressed == 'Y', 1, 0 )
imsup$ID <- sub( '-', '', imsup$ID ) # Remove the dashes from IDs
imsup$Group <- groups[imsup$Group]
all( imsup$ID == final_dat$ID ) # Check to make sure coherence between imsup data and immunology data

# Merge
analysis_dat <- merge( final_dat, imsup, by = c('ID', 'Group') )

# Export dataset
write_csv( analysis_dat, 'analysis_dat.csv' )

