### Script for running ubms dyn occupancy model on cpu cluster ###

setwd("xxx")

## add packages ##
library(ubms)
library(unmarked)
library(reshape2)

## load data ##
spp_data <- read.csv("data/combined_trait_data.csv", header = TRUE) # load species trait info
inc_spp_list <- as.character(spp_data[spp_data$major_group == "bee" & spp_data$No_gens_DATA_NOT_RELIABLE_FOR_HOVS == "Obligate Single", "Species"]) # create a list of species to model (single flight period bee species)
completed_spp_list <- read.csv("species_complete_spatiotemporal.csv", header = TRUE) # A list of species that have been run. Needed as some species failed previously
inc_spp_list <- inc_spp_list[!inc_spp_list%in%gsub(".rdata", "", completed_spp_list$species_done)] # ensure species that were already modelled are not modelled again. 

## set up index for cluster run ##
index <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
index

## select the species of interest for the given run
k <- inc_spp_list[index]

# load the species and predictor data (all in one rdata object)
load(paste("spp_data_2km/", k, "_formattedOccData_britain_2km.rdata", sep = "")) # load the data for the species in question
taxa_name <- k # record the species name
occDetdata <- spp_formattedOccData_britain_2km$occDetdata # split out the occupancy data
spp_vis <- spp_formattedOccData_britain_2km$spp_vis # extract the species visit info
occDetdata <- merge(occDetdata, spp_vis[,c("visit", taxa_name)]) # combine the above into a single file needed for ubms
names(occDetdata)[names(occDetdata) == taxa_name] <- "focal" # name the species of interest focal (essentially pick out the species of the current run)

# extract min and max year
min_year <- min(occDetdata$year)
max_year <- max(occDetdata$year)

# create a site look up
site_lookup <- unique(data.frame(occDetdata$site, occ_mod_site = as.numeric(as.factor(occDetdata$site))))
occDetdata$site <- as.numeric(as.factor(occDetdata$site))

# extract the max number of visits to a given site within a given year
site_visit_table <- unique(occDetdata[,c("visit", "site", "year")])
site_visits <- as.data.frame(table(paste(site_visit_table$site, site_visit_table$year, sep = "-")))
site_visits$Var1 <- as.character(site_visits$Var1)
max_site_year_visit <- max(site_visits$Freq)

# add visit number data
site_visit_table$site_year <- paste(site_visit_table$site, site_visit_table$year, sep = "-")
site_visit_table$visit_num = 1
site_visits_above_1 <- as.character(site_visits[site_visits$Freq > 1, "Var1"])

# loop through the site visit combos and update the visit number to be informative (i.e. which visit is it in the given year)
for (i in site_visits_above_1){
  site_visit_table[site_visit_table$site_year == i, "visit_num"] <- 1:site_visits[site_visits$Var1 == i,"Freq"]
}

# add visit number info back in the the main data frame
occDetdata_and_vis_num <- merge(occDetdata, site_visit_table[,c("visit", "visit_num")])
occDetdata_and_vis_num$year_vis <- paste0(occDetdata_and_vis_num$year, "_", occDetdata_and_vis_num$visit_num) # add a column to fill in the matrices below
occDetdata_and_vis_num$focal <- as.numeric(occDetdata_and_vis_num$focal) # convert this to binary rather than T/F

# convert to categorical list length
occDetdata_and_vis_num[occDetdata_and_vis_num$L %in% 2:3, "L"] <- 2
occDetdata_and_vis_num[occDetdata_and_vis_num$L > 3, "L"] <- 3  
occDetdata_and_vis_num$L <- as.character(occDetdata_and_vis_num$L)

## create a matrix to be filled with the presence/inferred absence, site and observation covariates data ##
# nrow = number of sites
# ncol = number of years multiplied by the maximum num of visits per year
# fill it with NAs
matrix_to_fill <- matrix(data = NA, nrow = length(unique(occDetdata$site)), ncol = max_site_year_visit*length(unique(occDetdata$year)))
row.names(matrix_to_fill) <- sort(unique(occDetdata$site)) # set the row names to be the site names
colnames(matrix_to_fill) <- paste0(rep(min_year:max_year, each = max_site_year_visit), "_", rep(1:max_site_year_visit, length(unique(occDetdata$year)))) # add column names which is the visit number within the given year

## create a matrix for each spatiotemporal terms ##
# create y - presence/absence
y_matrix <- matrix_to_fill
for(i in 1:nrow(occDetdata_and_vis_num)){
  print(paste(i, "of", nrow(occDetdata_and_vis_num)))
  y_matrix[rownames(y_matrix) == as.character(occDetdata_and_vis_num[i,"site"]), colnames(y_matrix) == as.character(occDetdata_and_vis_num[i,"year_vis"])] <- occDetdata_and_vis_num[i,"focal"]
}

colnames(y_matrix) <- paste0("det", "_", rep(min_year:max_year, each = max_site_year_visit), "_", rep(1:max_site_year_visit, length(unique(occDetdata$year)))) # add column names which is the site and the visit number to said site

# List length (observation model covariate)
LL_matrix <- matrix_to_fill
for(i in 1:nrow(occDetdata_and_vis_num)){
  print(paste(i, "of", nrow(occDetdata_and_vis_num)))
  LL_matrix[rownames(LL_matrix) == as.character(occDetdata_and_vis_num[i,"site"]), colnames(LL_matrix) == as.character(occDetdata_and_vis_num[i,"year_vis"])] <- occDetdata_and_vis_num[i,"L"]
}

colnames(LL_matrix) <- paste0("ll", "_", rep(min_year:max_year, each = max_site_year_visit), "_", rep(1:max_site_year_visit, length(unique(occDetdata$year)))) # add column names which is the site and the visit number to said site

# the following are annual-site covariates (not visit/observation level)
# local site mean annual temperature (state model covariate - persistence)

# need to update this to use the following:
load(paste("data/species_fp_summary_CHESS_data_2000/", k, "_tas_mean_1990_2015_anomaly.rdata", sep = "")) # load climate data 
clim_anom_data <- merge(final_clim_data, site_lookup, by.x = "gridref", by.y = "occDetdata.site", all.x = TRUE)
clim_anom_data <- clim_anom_data[!is.na(clim_anom_data$occ_mod_site),] # drop rows with NA in the anom_tas_mean column

spatiotemporal_temp <- matrix(data = NA, nrow = nrow(y_matrix), ncol = length(unique(occDetdata$year)))
rownames(spatiotemporal_temp) <- rownames(y_matrix)
colnames(spatiotemporal_temp) <- paste0("spatiotemporal_temp_", min_year:max_year) # add column names which is the visit number within the given year

for(i in 1:nrow(clim_anom_data)){
  print(paste(i, "of", nrow(clim_anom_data)))
  spatiotemporal_temp[rownames(spatiotemporal_temp) == as.character(clim_anom_data[i,"occ_mod_site"]), colnames(spatiotemporal_temp) == paste0("spatiotemporal_temp_", as.character(clim_anom_data[i,"year"]))] <- clim_anom_data[i,"tas_mean"]
}

spatiotemporal_temp <- scale(spatiotemporal_temp) # scale this variable

## create a matrix for the site parameters ##

# local mean temperature (state model covariate - persistence) - spatial term (not temporal)
load(paste("data/species_fp_summary_CHESS_data_2000/", k, "_MEAN_tas_mean_1990_2015.rdata", sep = "")) # load climate data
clim_mean_data <- merge(summ_data, site_lookup, by.x = "gridref", by.y = "occDetdata.site", all.x = TRUE)
clim_mean_data <- clim_mean_data[!is.na(clim_mean_data$occ_mod_site),] # drop rows with NA in the anom_tas_mean column
temp_mean <- matrix(data = NA, nrow = nrow(y_matrix), ncol = 1)
rownames(temp_mean) <- rownames(y_matrix)

for(i in 1:nrow(clim_mean_data)){
  print(paste(i, "of", nrow(clim_mean_data)))
  temp_mean[rownames(temp_mean) == as.character(clim_mean_data[i,"occ_mod_site"])] <- clim_mean_data[i,"tas_mean"]
}

colnames(temp_mean) <- "temp_mean" # add column names which is the site and the visit number to said site
temp_mean <- scale(temp_mean) # scale this variable

## combine the data into a single data frame and convert into a unmarkedMultiFrame ##

# add year identifiers
years <- as.character(min_year:max_year)
years <- matrix(years, nrow(matrix_to_fill), length(unique(occDetdata$year)), byrow=TRUE)

# convert matrices to dataframes
y_data <- as.data.frame(y_matrix)
LL_data <- as.data.frame(LL_matrix)

spatiotemporal_temp_data <- as.data.frame(spatiotemporal_temp)
temp_mean_data <- as.data.frame(temp_mean)

## create unmarked multiframe
umf_bee <- unmarkedMultFrame(y=y_matrix,
                             siteCovs = temp_mean_data,
                             yearlySiteCovs = list(year = years,
                                                   spatiotemporal_temp = spatiotemporal_temp_data),
                             obsCovs = list(LL = LL_matrix),
                             numPrimary = length(unique(occDetdata$year)))


## Run and save the dynamic occupancy model
fit_bee_4 <- stan_colext(psiformula = ~temp_mean, gammaformula = ~year-1, epsilonformula = ~spatiotemporal_temp, pformula = ~year + LL, umf_bee, chains=3, iter=3000)

save(fit_bee_4, file = paste0("output/dyn_occ_2km_spatiotemporal/", k, ".rdata"))

quit(save = 'no', runLast = FALSE)