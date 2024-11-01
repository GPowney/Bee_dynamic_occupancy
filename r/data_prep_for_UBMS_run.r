###### Data prep for dynamic occupancy model ######

#### Details ####
# Aim is to create a dataset to ready to be fed into a dynamic occupancy. First step is cleaning the species data,
# This will be records >= 1990, 2km scale, GB only. Next we need to convert the data into the format required by the 
# dynamic occupany model, to do this we run the formatOccData funciton from sparta. We then need to add climate 
# metrics to the occData dataframe. 

# add packages
require(BRCmap)
require(sparta)

## load data ##
taxa_data <- read.table("../raw_records_data.txt", header = TRUE) # load raw occurrence records

## drop non bee species ##
full_spp_list <- read.csv("../Aculeate_Species_Family_Group.csv", header = TRUE) # load species/group info table

taxa_run <- "Bee" # which broad taxonomic group are we running?
scale_run <- 2000 # set the spatial scale of the analysis
spp_list <- full_spp_list[full_spp_list$Popular_group_Name == taxa_run,] # select species of interest
taxa_data <- taxa_data[tolower(taxa_data$BWARS_SPECIES) %in% tolower(spp_list$FullName),] # filter records to match species of interest

## clean the data ##
# check precision of grid refs - Note all are 1 x 1km maximum
taxa_data$precision <- det_gr_precision(taxa_data$TO_GRIDREF)

# take all 1km prec or lower prec - save them to 1km - Note max prec is 1000 for the bee dataset
if(!all(taxa_data$precision == 1000)){
  taxa_data <- taxa_data[taxa_data$precision <= 1000, ] # monad or finer
  taxa_data$TO_GRIDREF <- reformat_gr(taxa_data$TO_GRIDREF, prec_out = 1000)
}

taxa_data <- unique(taxa_data) # remove duplicates 

# take all 1km prec or lower prec - reformat them to 2 x 2km grid square scale
taxa_data$TO_GRIDREF_Xkm <- reformat_gr(taxa_data$TO_GRIDREF, prec_out = scale_run)

# check dates are in the correct format
if(is.factor(taxa_data$TO_STARTDATE)){
  taxa_data$TO_STARTDATE <- as.Date(taxa_data$TO_STARTDATE, format = "%d/%m/%Y")
}

if(is.factor(taxa_data$TO_ENDDATE)){
  taxa_data$TO_ENDDATE <- as.Date(taxa_data$TO_ENDDATE, format = "%d/%m/%Y")
}

# subset to data with day precision
taxa_data <- taxa_data[taxa_data$TO_STARTDATE == taxa_data$TO_ENDDATE,]

# create single date column - ensure only day precision records
taxa_data$time_period <- as.Date(taxa_data$TO_STARTDATE, format = "%d/%m/%Y")

# add year column 
taxa_data$year <- format(taxa_data$time_period,"%Y") # add a year column

# subset to years of interest
taxa_data <- taxa_data[taxa_data$year>=1990 & taxa_data$year<=2015,]

## Drop non GB records ##
cn_info <- read.csv("../../../sq1km_country_id_border_dropped.csv", header = TRUE) # add country cell id info
taxa_data <- taxa_data[taxa_data$TO_GRIDREF %in% cn_info[cn_info$ENGLAND == 1 | cn_info$WALES == 1 | cn_info$SCOTLAND == 1, "SQ1_SQUARE"],] # exclude non GB records

# drop duplicate rows
taxa_data <- unique(taxa_data)

# drop unwanted columns
taxa_data <- taxa_data[c("BWARS_SPECIES", "TO_GRIDREF_Xkm", "time_period", "year")]
names(taxa_data) <- gsub("TO_GRIDREF_Xkm", "TO_GRIDREF", names(taxa_data)) # change column name

# save data for risk-of-bias assessment #
save(taxa_data, file = paste0("..\\bee_taxa_data_", scale_run, ".rdata"))

## Create the object required for the occupancy model
bee_formattedOccData_britain_2km <- formatOccData(taxa = taxa_data$BWARS_SPECIES,  
                                                  site = taxa_data$TO_GRIDREF,
                                                  survey = taxa_data$time_period)

bee_formattedOccData_britain_2km$spp.list <- names(bee_formattedOccData_britain_2km$spp_vis)[2:length(bee_formattedOccData_britain_2km$spp_vis)] # First column is "visit", not a species

## Add site and year specific climate data ##
spp_data <- read.csv("../combined_trait_data.csv", header = TRUE) # load trait data that's needed to extract species with a single flight period only
inc_spp_list <- as.character(spp_data[spp_data$major_group == "bee" & spp_data$No_gens_DATA_NOT_RELIABLE_FOR_HOVS == "Obligate Single", "Species"]) # extract bee species with a single flight period only

# loop through species create a bee_formattedOccData_britain_2km.rdata for each species #
for (k in inc_spp_list){
  
  load(paste("data/species_fp_summary_CHESS_data_2000/", k, "_tas_mean_1990_2015_anomaly.rdata", sep="")) # load species specific climate data (final_clim_data)
  names(final_clim_data) <- gsub("gridref", "site", names(final_clim_data)) # change column name
  names(final_clim_data) <- gsub("year", "TP", names(final_clim_data)) # change column name
  
  spp_formattedOccData_britain_2km <- NULL # clear previous version
  spp_formattedOccData_britain_2km <- bee_formattedOccData_britain_2km # create the species-specific climate/occupancy dataset
  
  # drop rows from both spp_vis and occDetData without climate data, then merge climate data into occDetData #
  occData <- merge(spp_formattedOccData_britain_2km$occDetdata, spp_formattedOccData_britain_2km$spp_vis) # add species info back into the occurrence data
  occData <- merge(occData, final_clim_data) # merge climate data with reformatted occurence data
  spp_formattedOccData_britain_2km$occDetdata <- occData[,c("visit", "site", "L", "TP","tas_mean", "anom_tas_mean")] # take columns of interest for the occupancy data
  spp_formattedOccData_britain_2km$spp_vis <- occData[,c("visit", spp_formattedOccData_britain_2km$spp.list)] # create the spp_vis data
  names(spp_formattedOccData_britain_2km$occDetdata) <- gsub("TP", "year", names(spp_formattedOccData_britain_2km$occDetdata)) # change column name
   
  # save the resulting visitData table for use on the cpu cluster
  save(spp_formattedOccData_britain_2km, file = paste("../", k, "_formattedOccData_britain_2km.rdata", sep = "")) # Save the formattedOccData object for the cluster run
  
}
