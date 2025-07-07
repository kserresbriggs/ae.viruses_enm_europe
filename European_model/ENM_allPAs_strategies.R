# Geospatial libraries
library(ncdf4)
library(exactextractr)
library(raster)
library(terra)
library(sf)
library(sp)

# Data manipulation libraries
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readxl)
library(openxlsx)
library(overlapping)

# Statistical modeling libraries
library(gmodels)
library(devtools)
library(USE)
library(dismo)
library(gbm)
library(blockCV)
library(ecospat)

# Plotting libraries
library(ggplot2)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(hrbrthemes)
library(scales)

# Parallel processing libraries
library(furrr)
library(future.apply)

setwd("/Users/kylaserres/Dropbox/PhD/Research projects/ENM-arboviral-projections/Scripts_and_Data/github/")

################################################################################
#                    I. LOADING ALL DATA
################################################################################

# 1. Europe NUTS3 shapefile preparation 
########################################

nutsM = crop(shapefile("./Shapefiles/MOOD_EUROPE_NUTS3_WillyW/moodpolygonsjanmasks21/VectornetMAPforMOODjan21.shp"), extent(-10,33,34.5,72)) # modified NUT3 shapefile
# Only keep EU areas
nutsM = subset(nutsM, !CountryISO%in%c("MA","DZ","TN","MT","TR","CI","MD","UA","RU","FO","IS","GE","BY","IS","GL","FO","CY","SJ"))
nutsM_sf = st_as_sf(nutsM)
# load correspondences between NUTS3 areas and NUTSM areas 
correspondences = shapefile("./Shapefiles/MOOD_EUROPE_NUTS3_WillyW/moodpolygonsjanmasks21/vnMOODdatamapcodejoinpoly.shp")@data[,c("DATLOCODE","MAPLOCODE")]

# get nutsM data
nutsM_IDs = as.data.frame(matrix(nrow = dim(nutsM@data)[1], ncol = 3))
colnames(nutsM_IDs) = c("country","region_name", "nutM_code"); nutsM_IDs$nutM_code = nutsM@data[,"LocationCo"]
nutsM_IDs$country = nutsM@data[,"CountryNam"]; nutsM_IDs$region_name = nutsM@data[,"LocationNa"]

# 2. Load all presence records in Europe
########################################
Aedes_presences = read.csv("European_model/Occurrence_data/ENM_nuts_presences.csv")


################################################################################
#                    II. PSEUDO-ABSENCE SAMPLING METHODS
################################################################################

# 1. Group presence data by years
##########################################

# In order to not give more importance to same areas where presences have been 
# recorded accross multiple years, we will only keep the year range so that the 
# environmental variables can be averaged across that range

# create the funtion
years_by_nuts = function(df) {
  df %>%
    group_by(nutM_code) %>%
    summarize(years = list(sort(unique(year))), .groups = "drop") %>%
    rowwise() %>%
    mutate(
      groups = list({
        y = years
        print(paste("Processing ID:", nutM_code))
        print(paste("All years:", paste(y, collapse = ", ")))
        
        if (length(y) == 1) {
          print("Single year")
          list(y)
        } else {
          breaks = c(0, which(diff(y) >= 10), length(y))
          print(paste("Year gap positions:", paste(breaks, collapse = ", ")))
          lapply(seq_along(breaks[-1]), function(i) y[(breaks[i] + 1):breaks[i + 1]])
        }
      })
    ) %>%
    unnest(groups) %>%
    mutate(
      n_years = lengths(groups),
      year = ifelse(n_years == 1, groups, NA_integer_),
      year_start = ifelse(n_years > 1, map_int(groups, min), NA_integer_),
      year_end = ifelse(n_years > 1, map_int(groups, max), NA_integer_)
    ) %>%
    select(nutM_code, year, year_start, year_end) %>%
    {
      .
    }
}


# apply function to the Aedes presences
ae.presences_agg = years_by_nuts(Aedes_presences)
ae.presences_agg$response = 1 # add response column

# 2. Create pseudo-absence shapefile
####################################

# For each sampling strategy, the adjacent polygons to presence sites are removed

# find adjacent polygons to each disease presence location
ae.pres_geoms = which(nutsM_sf$LocationCo %in% ae.presences_agg$nutM_code)
ae.pres_geoms = nutsM_sf[ae.pres_geoms, ]
ae.adjacent = st_intersects(ae.pres_geoms$geometry, nutsM_sf$geometry)
ae.adj = unlist(ae.adjacent)
ae.adj_geoms = nutsM_sf[ae.adj, ]

# plot to visualize (control check)
plotting = FALSE
if (plotting==TRUE) {
  plot(contour, lwd=0.4, border="gray30", col=NA)
  plot(nutsM, col= "gray95", border="gray50", lwd=0.1, add=T) 
  plot(ae.adj_geoms$geometry, col= "gray80",border="gray50", lwd=0.1, add=T)
  plot(ae.pres_geoms$geometry, col= "#D73027",border="gray50", lwd=0.1, add=T)
}


# removing adjacent polygons from the possible pseudo-absence areas
absences = setdiff(nutsM_IDs$nutM_code, ae.presences_agg) # select all NUTM areas that are not in presences
Aedes_abs = data.frame(
  nutM_code = absences, 
  year = NA
)

remove_rows = which(Aedes_abs$nutM_code %in% ae.adj_geoms$LocationCo) 
Aedes_abs = Aedes_abs[(-remove_rows), ] # remove adjacent polygons from the possible pseudo-absences
NA_ids = c("EL304","EL307","EL413","EL421","EL422","EL624","ES531","ES533",
           "ES630","ES640","GG","GI","IM","JE","MEG125354","MEG125368","UKM65") # these IDs correspond to small islands in Europe, 
# unfortunately the landcover data is unavailable for these areas (remove from possible absences)
Aedes_abs = Aedes_abs[!Aedes_abs$nutM_code %in% NA_ids, ] # This is the final sp objet that includes NUTS areas where PAs can be sampled from


# For all sampling strategies, the pseudo-absences must also be associated to a certain
# year as we are training the model temporally. The sampled year is weighted on
# the presence year frequency probabilities

# presence outbreak dates (year of infection)
# remove duplicated presences
Aedes_presences.unique = Aedes_presences[!duplicated(Aedes_presences[c("year", "nutM_code")]), ]
# get the years of presences
presence_years = Aedes_presences.unique$year
year_counts = table(presence_years)  # Count occurrences of each year
year_probs = year_counts / sum(year_counts)  # Normalize to get probabilities
years = as.integer(names(year_probs))
max_year = max(years)

# year ranges 
n_total = nrow(ae.presences_agg)
has_range = !is.na(ae.presences_agg$year_start) & !is.na(ae.presences_agg$year_end)
n_ranges = sum(has_range)
p_ranges = n_ranges / n_total

# 4. Get empirical distribution of range lengths from presence data
range_lengths = ae.presences_agg$year_end[has_range] - ae.presences_agg$year_start[has_range] + 1
range_length_table = table(range_lengths)
range_length_probs = range_length_table / sum(range_length_table)
possible_lengths = as.integer(names(range_length_probs))




# 3. All pseudo-absence sampling methods
########################################

## 3A. Random PA sampling strategy ##
# Pseudo-absence sampling strategies to be evaluated within the model training loop
random_sampling = quote({
  indices1 = seq_len(nrow(Aedes_abs))
  selected_absences = sample(indices1, ratios[[r]] * length(ae.presences_agg$nutM_code), replace = F) # set ratio (1:1, 1:3, 1:6)
  absences = Aedes_abs[selected_absences, ]
  
  n_abs = nrow(absences)
  n_range_abs = round(n_abs * p_ranges)
  range_indices = sample(seq_len(n_abs), size = n_range_abs, replace = FALSE)
  single_indices = setdiff(seq_len(n_abs), range_indices)
  
  absences$year = NA_integer_
  absences$year_start = NA_integer_
  absences$year_end = NA_integer_
  
  # Assign ranges
  absences$year_start[range_indices] = sample(
    years, size = n_range_abs, replace = TRUE, prob = year_probs
  )
  sampled_lengths = sample(
    possible_lengths, size = n_range_abs, replace = TRUE, prob = range_length_probs
  )
  absences$year_end[range_indices] = pmin(absences$year_start[range_indices] + sampled_lengths - 1, max_year)
  
  # Assign single years
  absences$year[single_indices] = sample(
    years, size = length(single_indices), replace = TRUE, prob = year_probs
  )
  
  # Add response
  absences$response = 0
  
  # Combine with presences
  data = rbind(ae.presences_agg, absences)
})


## 3B. Density-weighted geographic PA sampling strategy ##
# the Aedes Albopictus ecological suitability raster is used to 
# sample PAs in geographic regions where there the vector is abundant

# Load GeoTIFF aedes albopictus ecological suitability (M.Kraemer)
geo_density_sampling = quote({
  aedes = raster("./Rasters/albopictus_suitability_Kraemer2015.tif")
  aedes_cropped = crop(aedes, nutsM) # crop to Europe
  
  # extract suitability values per nutM area
  aedes_nutsM = exact_extract(aedes_cropped, nutsM, fun ="sum")
  nutsM_data = nutsM_IDs
  nutsM_data$aedes_suitability = aedes_nutsM
  
  # Compute weighted probabilities
  Aedes_abs$aedes_suitability = NA
  Aedes_abs$probabilities = NA
  for (a in seq_len(nrow(Aedes_abs))) {
    index = which(nutsM_data$nutM_code == Aedes_abs[a, "nutM_code"])
    Aedes_abs[a, "aedes_suitability"] = nutsM_data[index, "aedes_suitability"]
  }
  Aedes_abs$probabilities = Aedes_abs$aedes_suitability / sum(Aedes_abs$aedes_suitability)
  
  # the following code is to be in introduced in large loop to sample different PA
  # for each replicate.
  indices2 = seq_len(nrow(Aedes_abs))
  selected_absences = sample(indices2, ratios[[r]] * length(ae.presences_agg$nutM_code), replace = F, prob = Aedes_abs$probabilities) # set ratio (1:1, 1:3, 1:6)
  absences = Aedes_abs[selected_absences, ]
  
  n_abs = nrow(absences)
  n_range_abs = round(n_abs * p_ranges)
  range_indices = sample(seq_len(n_abs), size = n_range_abs, replace = FALSE)
  single_indices = setdiff(seq_len(n_abs), range_indices)
  
  absences$year = NA_integer_
  absences$year_start = NA_integer_
  absences$year_end = NA_integer_
  
  # Assign ranges
  absences$year_start[range_indices] = sample(
    years, size = n_range_abs, replace = TRUE, prob = year_probs
  )
  sampled_lengths = sample(
    possible_lengths, size = n_range_abs, replace = TRUE, prob = range_length_probs
  )
  absences$year_end[range_indices] = pmin(absences$year_start[range_indices] + sampled_lengths - 1, max_year)
  
  # Assign single years
  absences$year[single_indices] = sample(
    years, size = length(single_indices), replace = TRUE, prob = year_probs
  )
  
  # Add response
  absences$response = 0
  absences = absences[, !(colnames(absences) %in% c("aedes_suitability", "probabilities"))]
  Aedes_abs = Aedes_abs[, !(colnames(Aedes_abs) %in% c("aedes_suitability", "probabilities"))]
  
  # Combine with presences
  data = rbind(ae.presences_agg, absences)  
})

## 3C. Surveillance raster PA sampling strategy ##
# the surveillance suitability raster is used to 
# sample PAs in geographic regions where there arbovirus surveillance is robust and where presences have not been reported

# Load GeoTIFF surveillance raster Lim et al 2025
surveillance_sampling = quote({
  rast_surv = raster("./Rasters/Lim_et_al_2025_data/outputs/Rasters/Surveillance_map_wmean.tif")
  rast_surv_cropped = crop(rast_surv, nutsM) # crop to Europe
  
  # extract suitability values per nutM area
  surv_nutsM = exact_extract(rast_surv_cropped, nutsM, fun="mean")
  surv.nutsM_data = nutsM_IDs
  surv.nutsM_data$surv = surv_nutsM 
  surv.nutsM_data$surv[is.na(surv.nutsM_data$surv)] = 0
  
  # Compute weighted probabilities
  Aedes_abs$surv = NA
  Aedes_abs$Probabilities = NA
  for (m in seq_len(nrow(Aedes_abs))) {
    index = which(surv.nutsM_data$nutM_code == Aedes_abs[m, "nutM_code"])
    Aedes_abs[m, "surv"] = surv.nutsM_data[index, "surv"]
  }
  Aedes_abs$Probabilities = Aedes_abs$surv / sum(Aedes_abs$surv)
  
  # the following code is to be in introduced in large loop to sample different PA
  # for each replicate.
  indices3 = seq_len(nrow(Aedes_abs))
  selected_absences = sample(indices3, ratios[[r]] * length(ae.presences_agg$nutM_code), replace = F, prob = Aedes_abs$Probabilities) # set ratio (1:1, 1:3, 1:6)
  absences = Aedes_abs[selected_absences, ]
  
  n_abs = nrow(absences)
  n_range_abs = round(n_abs * p_ranges)
  range_indices = sample(seq_len(n_abs), size = n_range_abs, replace = FALSE)
  single_indices = setdiff(seq_len(n_abs), range_indices)
  
  absences$year = NA_integer_
  absences$year_start = NA_integer_
  absences$year_end = NA_integer_
  
  # Assign ranges
  absences$year_start[range_indices] = sample(
    years, size = n_range_abs, replace = TRUE, prob = year_probs
  )
  sampled_lengths = sample(
    possible_lengths, size = n_range_abs, replace = TRUE, prob = range_length_probs
  )
  absences$year_end[range_indices] = pmin(absences$year_start[range_indices] + sampled_lengths - 1, max_year)
  
  # Assign single years
  absences$year[single_indices] = sample(
    years, size = length(single_indices), replace = TRUE, prob = year_probs
  )
  
  # Add response
  absences$response = 0
  absences = absences[, !(colnames(absences) %in% c("surv", "Probabilities"))]
  Aedes_abs = Aedes_abs[, !(colnames(Aedes_abs) %in% c("surv", "Probabilities"))]
  
  data = rbind(ae.presences_agg, absences)
})


## 3D. Environmental space PA sampling strategy ##
ES_setup = quote({
  # Load the reference raster (e.g., one of the ISIMIP climatic layers)
  reference = brick(paste0("./Rasters/Environmental_rasters/ISIMIP3_past/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"))
  reference = reference[[1]] # select the first layer
  reference = crop(reference, nutsM, snap="out") # crop to europe
  reference = mask(reference, nutsM)
  values(reference) = NA # Reference raster becomes NULL raster
  
  # duplicate reference raster
  rast_remove = reference
  
  # Use exact_extract to get cell indices for all presence geometries
  # All cells that are within a presence polygon are assigned a value of 1
  presence_indices = exact_extract(reference, ae.pres_geoms, include_cell = TRUE)
  indices = list()
  for (i in 1:length(presence_indices)){Idx = presence_indices[[i]][[2]];indices[[i]] = Idx}
  indices = unlist(indices, use.names = FALSE)
  # set presence cells of null raster to 1
  reference[indices] = 1
  
  # Use exact_extract to get cell indices for all adjacent polygons to presence geometries
  # All cells that are within an adjacent polygon are assigned a value of 1
  # this raster will be used to remove the adjacent polygons
  
  NA_ids = c("EL304","EL307","EL413","EL421","EL422","EL624","ES531","ES533",
             "ES630","ES640","GG","GI","IM","JE","MEG125354","MEG125368","UKM65") # these IDs correspond to small islands in Europe, 
  # unfortunately the landcover data is unavailable for these areas (remove from possible absences)
  na_geoms = which(nutsM_sf$LocationCo %in% NA_ids)
  na_geoms = nutsM_sf[na_geoms, ]
  Ae.adj_geoms = rbind(ae.adj_geoms, na_geoms) # add the small islands to the adjacent polygons
  
  adj_indices = exact_extract(reference, Ae.adj_geoms, include_cell = TRUE)
  Adj_indices = list()
  for (i in 1:length(adj_indices)){
    idx = adj_indices[[i]][[2]]; Adj_indices[[i]] = idx}
  Adj_indices = unlist(Adj_indices, use.names = FALSE)
  # set adjacent cells of second null raster to 2
  rast_remove[Adj_indices] = 2
  
  # Find cells where both presence and adjacent cells overlap (1 + 2)
  # remove the overlap, the final raster will only have cells corresponding to 
  # adjacent polygons
  overlap_mask = reference == 1 & rast_remove == 2
  rast_remove[overlap_mask] = NA
  
  # Extract data from the reference presence raster
  ref_data = as.data.frame(matrix(nrow=ncell(reference), ncol=4)); colnames(ref_data) = c("indices", "longitude", "latitude", "response")
  ref_data$indices = 1:ncell(reference) # get presence cell indices
  coords = xyFromCell(reference, 1:ncell(reference))
  ref_data$longitude = coords[, 1]; ref_data$latitude = coords[, 2] # get presence cell coordinates
  ref_data$response = getValues(reference) # presence cells (value = 1)
  
  myPres = ref_data[!is.na(ref_data$response), ]  # presences (value=1)
  myPres = st_as_sf(myPres, coords = c("longitude", "latitude")) # convert to sf object
  
  # Control check 
  plot = F
  if (plot==T){
    dev.new()
    plot(nutsM, col="gray95", border="gray50", lwd=0.1) 
    plot(rast_remove, col = "gray50", add=T, legend=F, axes=F)
    plot(reference, col = "tomato", add=T, legend=F, axes=F)
  }
  
  
  
  ## Load climatic rasters for the environmental space pseudo-absence sampling ##
  # load raster layers, subset over the past 20 years and by season 
  year_start = (1999 - 1901) * 12 + 1; year_end = (2019 - 1901 + 1) * 12
  temp = brick(paste0("./Rasters/Environmental_rasters/ISIMIP3_past/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"))
  temp = subset(temp, year_start:year_end)
  
  winter_indices = sort(c(12, 1, 2)) + rep((1999:2019 - 1999) * 12, each = 3); winter_indices = winter_indices[-c(1, 2, length(winter_indices))] # December (previous year), January, February
  spring_indices = sort(c(3,4,5) + rep((2000:2019-1999)*12, each= 3))
  summer_indices = sort(c(6,7,8) + rep((2000:2019-1999)*12, each= 3))
  fall_indices = sort(c(9,10,11) + rep((2000:2019-1999)*12, each= 3))
  
  tepmm_winter = mean(temp[[winter_indices]]) - 273.15
  tepmm_spring = mean(temp[[spring_indices]]) - 273.15
  tepmm_summer = mean(temp[[summer_indices]]) - 273.15
  tepmm_fall = mean(temp[[fall_indices]]) - 273.15
  
  precp = brick(paste0("./Rasters/Environmental_rasters/ISIMIP3_past/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_pr_1901_2021_monmean.nc"))
  precp = subset(precp, year_start:year_end)
  
  precp_winter = mean(precp[[winter_indices]]) * 60 * 60 * 24
  precp_spring = mean(precp[[spring_indices]]) * 60 * 60 * 24
  precp_summer = mean(precp[[summer_indices]]) * 60 * 60 * 24
  precp_fall = mean(precp[[fall_indices]]) * 60 * 60 * 24
  
  relh = brick(paste0("./Rasters/Environmental_rasters/ISIMIP3_past/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_hurs_1901_2021_monmean.nc"))
  relh = subset(relh, year_start:year_end)
  
  relh_winter = mean(relh[[winter_indices]])
  relh_spring = mean(relh[[spring_indices]]) 
  relh_summer = mean(relh[[summer_indices]])
  relh_fall = mean(relh[[fall_indices]])
  
  year_start_temp = 2000 - 1900;  year_end_temp = 2024 - 1900
  crops_temp = brick("./Rasters/Environmental_rasters/Landuse/croplands_LUH2-GCB2024_1901_2024.nc"); crops_temp = mean(crops_temp[[year_start_temp:year_end_temp]])
  past_temp = brick("./Rasters/Environmental_rasters/Landuse/pastures_LUH2-GCB2024_1901_2024.nc"); past_temp = mean(past_temp[[year_start_temp:year_end_temp]])
  ranges_temp = brick("./Rasters/Environmental_rasters/Landuse/rangeland_LUH2-GCB2024_1901_2024.nc"); ranges_temp = mean(ranges_temp[[year_start_temp:year_end_temp]])
  urban_temp = brick("./Rasters/Environmental_rasters/Landuse/urbanAreas_LUH2-GCB2024_1901_2024.nc"); urban_temp = mean(urban_temp[[year_start_temp:year_end_temp]])
  primForest_temp = brick("./Rasters/Environmental_rasters/Landuse/primaryForest_LUH2-GCB2024_1901_2024.nc"); primForest_temp = mean(primForest_temp[[year_start_temp:year_end_temp]])
  primNonForest_temp = brick("./Rasters/Environmental_rasters/Landuse/primaryNonF_LUH2-GCB2024_1901_2024.nc"); primNonForest_temp = mean(primNonForest_temp[[year_start_temp:year_end_temp]])
  secForest_temp = brick("./Rasters/Environmental_rasters/Landuse/secondaryForest_LUH2-GCB2024_1901_2024.nc"); secForest_temp = mean(secForest_temp[[year_start_temp:year_end_temp]])
  secNonForest_temp = brick("./Rasters/Environmental_rasters/Landuse/secondaryNonF_LUH2-GCB2024_1901_2024.nc"); secNonForest_temp = mean(secNonForest_temp[[year_start_temp:year_end_temp]])
  
  pop = brick("./Rasters/Environmental_rasters/population_histsoc_30arcmin_annual_1850_2021.nc", varname="total-population")
  year_start_temp = 2000 - 1900;  year_end_temp = 2019 - 1900
  pop_temp = mean(pop[[year_start_temp:year_end_temp]]); pop_temp = log10(pop_temp + 1)
  
  
  # group rasters
  envVariables = list()
  envVariables[[1]] = tepmm_winter
  envVariables[[2]] = tepmm_spring
  envVariables[[3]] =  tepmm_summer
  envVariables[[4]] =  tepmm_fall
  envVariables[[5]] =  precp_winter
  envVariables[[6]] =  precp_spring
  envVariables[[7]] =  precp_summer
  envVariables[[8]] =  precp_fall
  envVariables[[9]] =  relh_winter
  envVariables[[10]] =  relh_spring
  envVariables[[11]] =  relh_summer
  envVariables[[12]] =  relh_fall
  envVariables[[13]] = crops_temp
  envVariables[[14]] = past_temp
  envVariables[[15]] = ranges_temp
  envVariables[[16]] = urban_temp
  envVariables[[17]] = primForest_temp
  envVariables[[18]] = primNonForest_temp
  envVariables[[19]] = secForest_temp
  envVariables[[20]] = secNonForest_temp
  envVariables[[21]] = pop_temp
  
  # match extent to europe ans stack
  for (i in 1:length(envVariables)) envVariables[[i]] = crop(envVariables[[i]], nutsM, snap="out")
  for (i in 1:length(envVariables)) envVariables[[i]] = mask(envVariables[[i]], nutsM)
  variable_stack = stack(envVariables)
  variable_stack[rast_remove] = NA # remove cells corresponding to adjacent polygons from the rasters and thus the environmental space
  
  ## 3. Run the PCA
  PCA = USE::rastPCA(env.rast = variable_stack, nPC = 2, stand = TRUE)
  pca_result = PCA$pca ; pca_summary = summary(pca_result); print(pca_summary) 
  PCA_stack = c(PCA$PCs$PC1, PCA$PCs$PC2)
  
  
  ## 4. Find the optimal grid resolution for pseudo-absence sampling
  PCAstack_df = as.data.frame(PCA_stack, xy = T, na.rm = T)
  PCAstack_sp = st_as_sf(PCAstack_df, coords = c("PC1", "PC2"))
  opt_res = optimRes(sdf = PCAstack_sp, grid.res = seq(1, 15), showOpt = T)
  
  # 5. rasterize the Europe polygon shapefile to associate raster cells to a polygon ID
  rast = brick(paste0("./Rasters/Environmental_rasters/ISIMIP3_past/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"))
  rast = rast[[1]] # select the first layer to make the polygon file the same resolution and extent
  rast = crop(rast, nutsM, snap="out")
  rast = mask(rast, nutsM)
  
  nutsM$LocationCo_numeric = as.numeric(as.factor(nutsM$LocationCo)) # Convert polygon IDs to numeric indices
  rasterized_nutsM = rasterize(nutsM, rast, field="LocationCo_numeric")
  
  cell_ids = which(!is.na(values(rasterized_nutsM)))  # cell numbers with data
  polygon_indices = rasterized_nutsM[cell_ids]        # the numeric values (polygon IDs)
  lookup_table = unique(as.data.frame(nutsM)[, c("LocationCo_numeric", "LocationCo")])
  
})

# sampling strategy for Environmental space sampling
ES_sampling = quote({
  ps.abs = USE::paSampling(variable_stack,
                           pres = myPres,
                           thres = 0.75,
                           H = NULL,
                           grid.res = opt_res$Opt_res, 
                           prev = ratios[[r]], 
                           sub.ts = FALSE,
                           plot_proc = F)
  
  
  ABS = ps.abs[,c("x","y")]
  raster_extent = extent(variable_stack[[1]])
  
  # Create an empty raster with the same extent and resolution
  raster_cells = raster(raster_extent, nrows = nrow(variable_stack[[1]]), ncols = ncol(variable_stack[[1]]))
  
  # Initialize the raster with NAs (no data)
  raster_cells[] = NA
  
  # Set pseudo-absence points (ABS) as 1 in the raster grid
  # 'ABS' contains the coordinates of pseudo-absence points
  for(j in 1:nrow(ABS)) {
    x = ABS$x[j]
    y = ABS$y[j]
    
    # Use the raster cell coordinates to set the value to 1
    cell = cellFromXY(raster_cells, c(x, y))
    raster_cells[cell] = 1
  }
  
  df_cell.nutsM = as.data.frame(matrix(nrow = length(values(raster_cells)), ncol = 3))
  colnames(df_cell.nutsM) = c("longitude", "latitude", "PA_value")
  df_cell.nutsM$longitude = xFromCol(raster_cells, 1:ncol(raster_cells))
  df_cell.nutsM$latitude = yFromRow(raster_cells, 1:nrow(raster_cells))
  df_cell.nutsM$PA_value = getValues(raster_cells) # pseudo absence cells (value=1)
  df_cell.nutsM = na.omit(df_cell.nutsM)
  
  # Keep only the overlapping cells in both rasters (nutsM ids and pseudoabsences)
  overlap = mask(raster_cells, rasterized_nutsM)
  overlap_idx = which(!is.na(values(overlap)))
  
  # Extract the polygon IDs from rasterized_nutsM at those indices
  polygon_ids_numeric = values(rasterized_nutsM)[overlap_idx]
  
  # Match numeric IDs to NUTS3 ID codes
  matched_ids = lookup_table$LocationCo[match(polygon_ids_numeric, lookup_table$LocationCo_numeric)]
  df_cell.nutsM$polygon_ID = matched_ids 
  
  final_ps.abs_df = as.data.frame(matrix(nrow = dim(df_cell.nutsM[1]), ncol = 8))
  colnames(final_ps.abs_df) = c("disease", "nutM_code", "year", "year_start","year_end","region_name","case_numbers", "response")
  final_ps.abs_df$year = NA;final_ps.abs_df$year_start = NA; final_ps.abs_df$year_end = NA;
  final_ps.abs_df$nutM_code = df_cell.nutsM$polygon_ID;
  final_ps.abs_df$disease = "aedes transmitted virus"
  final_ps.abs_df$response = 0 ; final_ps.abs_df$case_numbers = NA
  final_ps.abs_df$region_name = nutsM_IDs$region_name[match(final_ps.abs_df$nutM_code, nutsM_IDs$nutM_code)]
  
  n_abs = nrow(final_ps.abs_df)
  n_range_abs = round(n_abs * p_ranges)
  range_indices = sample(seq_len(n_abs), size = n_range_abs, replace = FALSE)
  single_indices = setdiff(seq_len(n_abs), range_indices)
  
  
  # Assign ranges
  final_ps.abs_df$year_start[range_indices] = sample(
    years, size = n_range_abs, replace = TRUE, prob = year_probs
  )
  sampled_lengths = sample(
    possible_lengths, size = n_range_abs, replace = TRUE, prob = range_length_probs
  )
  final_ps.abs_df$year_end[range_indices] = pmin(final_ps.abs_df$year_start[range_indices] + sampled_lengths - 1, max_year)
  
  # Assign single years
  final_ps.abs_df$year[single_indices] = sample(
    years, size = length(single_indices), replace = TRUE, prob = year_probs
  )
  
  final_ps.abs_df = final_ps.abs_df[, !(colnames(final_ps.abs_df) %in% c("disease", "case_numbers"))]
  
  # Combine presences and absences
  data = rbind(ae.presences_agg, final_ps.abs_df)
})



################################################################################
#                   III. BRT MODEL TRAINING AND EVALUATION
################################################################################


# 1. Environmental data extraction
########################################

nruns = 30 # this will be used as the number of repetitions for the brt training
pa.strategy = c("random", "geo_ae.density", "surveillance","espace") # pseudo-absence sampling strategy
namesRatios = c("1.1", "1.3", "1.6") # ratio names
ratios = c(1,3,6) # corresponds to the ratio of desired sampled PAs (1:1, 1:3, 1:6)
ratios2 =c(3.7, 1.3, 0.64); # interpolated prevalence values for Environmental space sampling (to get correct ratio values)

croplands = brick("./Rasters/Environmental_rasters/Landuse/croplands_LUH2-GCB2024_1901_2024.nc")
pastures = brick("./Rasters/Environmental_rasters/Landuse/pastures_LUH2-GCB2024_1901_2024.nc")
rangelands = brick("./Rasters/Environmental_rasters/Landuse/rangeland_LUH2-GCB2024_1901_2024.nc")
urbanAreas = brick("./Rasters/Environmental_rasters/Landuse/urbanAreas_LUH2-GCB2024_1901_2024.nc")
primaryForest = brick("./Rasters/Environmental_rasters/Landuse/primaryForest_LUH2-GCB2024_1901_2024.nc")
primaryNonForest = brick("./Rasters/Environmental_rasters/Landuse/primaryNonF_LUH2-GCB2024_1901_2024.nc")
secondaryForest = brick("./Rasters/Environmental_rasters/Landuse/secondaryForest_LUH2-GCB2024_1901_2024.nc")
secondaryNonForest = brick("./Rasters/Environmental_rasters/Landuse/secondaryNonF_LUH2-GCB2024_1901_2024.nc")
population = brick("./Rasters/Environmental_rasters/population_histsoc_30arcmin_annual_1850_2021.nc", varname="total-population")

# run PCA for the environmental space sampling
eval(ES_setup)
# set parallele process
plan(multisession, workers = parallel::detectCores() - 1)

# Loop through all pseudo-absence strategies and ratios
for(s in 1:length(pa.strategy)){
  data_train = list()
  # load climate rasters
  temperature = brick(paste0("./Rasters/Environmental_rasters/ISIMIP3_past/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"))
  precipitation = brick(paste0("./Rasters/Environmental_rasters/ISIMIP3_past/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_pr_1901_2021_monmean.nc"))
  relativehumidity = brick(paste0("./Rasters/Environmental_rasters/ISIMIP3_past/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_hurs_1901_2021_monmean.nc"))
  
  # ratio loop
  for (r in 1:length(ratios)){
    data_train[[r]] = list()

    # repition loop 
    for(n in 1:nruns){
      if(s==1){eval(random_sampling)}
      if(s==2){eval(geo_density_sampling)}
      if(s==3){eval(surveillance_sampling)}
      if(s==4){ratios=ratios2; eval(ES_sampling)}
      
      ## data extraction ##
      
      data_for_brt = as.data.frame(matrix(NA,nrow=nrow(data), ncol =27))
      data_for_brt[,c(1:5)] = data[,c("nutM_code","year","year_start","year_end", "response")]
      colnames(data_for_brt) = c("nutM_code","year","year_start","year_end", "response",
                                 "temperature_winter","temperature_spring","temperature_summer","temperature_fall",
                                 "seasonal_variation","precipitation_winter","precipitation_spring","precipitation_summer","precipitation_fall",
                                 "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_fall",
                                 "croplands","pastures","rangelands","urbanAreas","primaryForest","primaryNonForest","secondaryForest",
                                 "secondaryNonForest", "population")
      
      geometry_matches = nutsM_sf$geometry[match(data_for_brt$nutM_code, nutsM_sf$LocationCo)]
      nutsM_train = st_sf(data_for_brt, geometry = geometry_matches)
      names_env = names(data_for_brt[,c(6:27)])
      
      cat("Starting data extraction", n, "of", nruns, "for PA sampling method", pa.strategy[s], "ratio", namesRatios[r])
      # parallel data extraction block 
      results = future_lapply(1:nrow(data_for_brt), function(k) {
        if (!is.na(data_for_brt$year[k])) { 
          # Temporal training - split climate variables into seasons
          message("temporal training on single year")
          year = as.numeric(data_for_brt[k,"year"]);  year2 = year
          if (year > 2021) year = 2021
          start = (year - 1901)*12+1; end = start + 11
          month_start = ((year-1) - 1901)*12+1
          winter_indices = sort(c(month_start + 11, month_start + 12, month_start + 13)) # December (previous year), January, February
          spring_indices = sort(c(month_start + 14, month_start + 15, month_start + 16)) # March, April, May
          summer_indices = sort(c(month_start + 17, month_start + 18, month_start + 19)) # June, July, August
          fall_indices   = sort(c(month_start + 20, month_start + 21, month_start + 22)) # September, October, November
          
          year_index = year2 - 1900 # landuse is updated and goes until 2024
          
          # seasonal variation
          yearly_temp = temperature[[start:end]]-273.15
          temp_raster = crop(yearly_temp, nutsM, snap = "out")
          temp_raster = mask(temp_raster, nutsM)
          
          max = cellStats(temp_raster, stat = "max", na.rm = TRUE)
          min = cellStats(temp_raster, stat = "min", na.rm = TRUE)
          
          hotest_indice = which.max(max)
          hotest_indice = start + hotest_indice - 1
          
          coldest_indice = which.min(min)
          coldest_indice = start + coldest_indice - 1
          
          
        } else {
          message("temporal training on average year range")
          year_start = as.numeric(data_for_brt[k,"year_start"]); year_end = as.numeric(data_for_brt[k,"year_end"])
          year_start_temp = year_start - 1900;  year_end_temp = year_end - 1900
          year_index = year_start_temp:year_end_temp
          
          if (year_start > 2021) year_start = 2021; if (year_end > 2021) year_end = 2021

          winter_indices=c()
          spring_indices=c()
          summer_indices=c()
          fall_indices  =c()
          
          for (yr in year_start:year_end) {
            month_start=((yr - 1) - 1901) * 12 + 1
            winter_indices=c(winter_indices, month_start + 11, month_start + 12, month_start + 13) # December (previous year), January, February
            spring_indices=c(spring_indices, month_start + 14, month_start + 15, month_start + 16) # March, April, May
            summer_indices=c(summer_indices, month_start + 17, month_start + 18, month_start + 19) # June, July, August
            fall_indices  =c(fall_indices,   month_start + 20, month_start + 21, month_start + 22) # September, October, November
          }
          
          # seasonal variation
          allyears = year_start:year_end
          hotest_indice = list()
          coldest_indice = list()
          for (i in 1:length(allyears)){
            yr = allyears[[i]]
            start = (yr - 1901)*12+1; end = start + 11
            yearly_temp = temperature[[start:end]]-273.15
            temp_raster = crop(yearly_temp, nutsM, snap = "out")
            temp_raster = mask(temp_raster, nutsM)
            
            max = cellStats(temp_raster, stat = "max", na.rm = TRUE)
            min = cellStats(temp_raster, stat = "min", na.rm = TRUE)
            
            Hotest_indice = which.max(max)
            hotest_indice[[i]] = start + Hotest_indice - 1
            Coldest_indice = which.min(min)
            coldest_indice[[i]] = start + Coldest_indice - 1
            
          }
        }
        
        # Extract environmental variables for the given year or year range
        temperature_temp_winter = mean(temperature[[winter_indices]]) - 273.15 # conversion to Celsius
        temperature_temp_spring = mean(temperature[[spring_indices]]) - 273.15
        temperature_temp_summer = mean(temperature[[summer_indices]]) - 273.15
        temperature_temp_fall   = mean(temperature[[fall_indices]]) - 273.15
        
        hotmonth_temp = mean(temperature[[hotest_indice]]) - 273.15
        coldmonth_temp = mean(temperature[[coldest_indice]]) - 273.15
        seasonal_var_temp = hotmonth_temp - coldmonth_temp
        
        precipitation_temp_winter = mean(precipitation[[winter_indices]]) * 60 * 60 * 24 # conversion to kg/m2/day
        precipitation_temp_spring = mean(precipitation[[spring_indices]]) * 60 * 60 * 24
        precipitation_temp_summer = mean(precipitation[[summer_indices]]) * 60 * 60 * 24
        precipitation_temp_fall   = mean(precipitation[[fall_indices]]) * 60 * 60 * 24
        
        relhumidity_temp_winter = mean(relativehumidity[[winter_indices]])
        relhumidity_temp_spring = mean(relativehumidity[[spring_indices]])
        relhumidity_temp_summer = mean(relativehumidity[[summer_indices]])
        relhumidity_temp_fall   = mean(relativehumidity[[fall_indices]])
        
        croplands_temp = mean(croplands[[year_index]])
        pastures_temp = mean(pastures[[year_index]])
        urbanAreas_temp = mean(urbanAreas[[year_index]])
        rangelands_temp = mean(rangelands[[year_index]])
        primaryForest_temp = mean(primaryForest[[year_index]])
        primaryNonForest_temp = mean(primaryNonForest[[year_index]])
        secondaryForest_temp = mean(secondaryForest[[year_index]])
        secondaryNonForest_temp = mean(secondaryNonForest[[year_index]])
        population_temp = mean(population[[year_index]])
        
        envVariables = list(
          temperature_temp_winter, temperature_temp_spring, temperature_temp_summer, temperature_temp_fall,
          seasonal_var_temp, precipitation_temp_winter, precipitation_temp_spring, precipitation_temp_summer, precipitation_temp_fall,
          relhumidity_temp_winter, relhumidity_temp_spring, relhumidity_temp_summer, relhumidity_temp_fall,
          croplands_temp, pastures_temp, urbanAreas_temp, rangelands_temp,
          primaryForest_temp, primaryNonForest_temp, secondaryForest_temp, secondaryNonForest_temp,
          population_temp
        )
        
        envVariables = lapply(envVariables, crop, y = nutsM, snap = "out")
        envVariables = lapply(envVariables, mask, mask = nutsM)
        
        vals = numeric(length(envVariables))
        for (i in seq_along(envVariables)) {
          if (i == length(envVariables)) {
            vals[i] = log(exact_extract(envVariables[[i]], nutsM_train[k,], fun = "sum") + 1)
          } else {
            vals[i] = exact_extract(envVariables[[i]], nutsM_train[k,], fun = "mean")
          }
        }
        
        return(vals)
      })
      data_for_brt[, names_env] = do.call(rbind, results)
      data_train[[r]][[n]] = data_for_brt
      
      cat("Data extraction run", n, "of", nruns, "for PA sampling method", 
          pa.strategy[s], "ratio", namesRatios[r], "completed.\n")
    } # end of repitiions loop for environmental data extraction
  }
  saveRDS(data_train, paste0(file = "European_model/All_the_brt_models/PA_strategies_brts/training_data", pa.strategy[s],".rds"))
}


# 2. BRT model training
#######################

ratios = c(1,3,6) # corresponds to the ratio of desired sampled PAs (1:1, 1:3, 1:6)

# Calculate the mean block size for spatial cross-validation per PA strategy and ratio
file_paths = c(
  "./European_model/All_the_brt_models/PA_strategies_brts/training_data_random.rds",
  "./European_model/All_the_brt_models/PA_strategies_brts/training_data_geo_ae.density.rds",
  "./European_model/All_the_brt_models/PA_strategies_brts/training_data_surveillance.rds",
  "./European_model/All_the_brt_models/PA_strategies_brts/training_data_espace.rds"
)
#data_train4 = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_espace.rds")

Data = lapply(file_paths, readRDS)
# keep one training dataset per PA sampling method
Data = lapply(Data, function(three_list) {
  lapply(three_list, function(sublist) {
    sublist[[1]]
  })
})

centroids = as.data.frame(matrix(nrow=length(nutsM), ncol=3)); colnames(centroids) = c("NUTM","longitude","latitude")  
centroids$NUTM = nutsM@data$LocationCo; coords = coordinates(nutsM)  
centroids$longitude = coords[,1]; centroids$latitude = coords[,2]  

plottingCorrelogram = FALSE

for (i in 1:length(Data)){
  for (j in 1:length(Data[[i]])){
    data = Data[[i]][[j]]
    matching_centroids = as.data.frame(matrix(nrow = nrow(data), ncol = 3)); colnames(matching_centroids) = c("NUTM","longitude", "latitude")
    matching_centroids$NUTM = data$nutM_code 
    matching_centroids$longitude = centroids$longitude[match(matching_centroids$NUTM, centroids$NUTM)]; matching_centroids$latitude = centroids$latitude[match(matching_centroids$NUTM, centroids$NUTM)]
    correlogram = ncf::correlog(matching_centroids[,2], matching_centroids[,3], data[,"response"], na.rm=T, increment=10, resamp=0, latlon=T)

    if (plottingCorrelogram == TRUE)
      {
        dev.new(width=4.5, height=3); par(mar=c(2.5,2.2,1.0,1.5))
        plot(correlogram$mean.of.class[-1], correlogram$correlation[-1], ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-0.4,1.0), xlim=c(0,3500))
        abline(h=0, lwd=0.5, col="red", lty=2)
        lines(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, col="gray30")
        points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.25, col="gray90", pch=16)
        points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.25, col="gray30", pch=1)
        axis(side=1, pos=-0.4, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,-0.05,0), at=seq(0,9000,1000))
        axis(side=2, pos=0, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.18,0), at=seq(-0.4,1,0.2))
        title(xlab="distance (km2)", cex.lab=0.7, mgp=c(0.3,0,0), col.lab="gray30")
        title(ylab="correlation", cex.lab=0.7, mgp=c(0.4,0,0), col.lab="gray30")
    }
  }
}

# chose block size per strategy and ratio based on plotted correlogram (when curve intercepts 0)
blocks = list(Random = list(600, 800, 800),
     Ae = list(500,700,700),
     surv = list(600, 500,800),
     es = list(1000, 1100, 1000))

for(s in 1:length(pa.strategy)){
  brt_model_scvs = list() # brt with spatial cross-validations (SCVs)
  AUCs = matrix(nrow=nruns, ncol= length(ratios)); colnames(AUCs) = paste0(pa.strategy[[s]], "_", namesRatios)
  SIppcs = matrix(nrow=nruns, ncol=length(ratios)); colnames(SIppcs) = paste0(pa.strategy[[s]], "_", namesRatios)
  thresholds = matrix(nrow=nruns, ncol=length(ratios)); colnames(thresholds) = paste0(pa.strategy[[s]], "_", namesRatios)
  tabs_list1 = list()
  
  if(s==1){data_train = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_random.rds")}
  if(s==2){data_train = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_geo_ae.density.rds")}
  if(s==3){data_train = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_surveillance.rds")}
  if(s==4){data_train = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_espace.rds")}
  
  for (r in 1:length(ratios)){
    brt_model_scvs[[r]] = list()
    tabs_list1[[r]] = list()
    current_block = blocks[[s]][[r]]
    for(b in 1:nruns){
      current_data = data_train[[r]][[b]]
      # BRT model parameters
      gbm.x = 6:27 # environmental predictors
      gbm.y = colnames(current_data)[5] # response variable
      offset = NULL
      tree.complexity = 5 # "tc" = number of nodes in the trees
      learning.rate = 0.001 # "lr" = contribution of each tree to the growing model
      bag.fraction = 0.80 # proportion of data used to train a given tree
      site.weights = rep(1, dim(current_data)[1])
      var.monotone = rep(0, length(gbm.x))
      prev.stratify = TRUE
      family = "bernoulli"
      n.trees = 100 # initial number of trees
      step.size = 10 # interval at which the predictive deviance is computed and logged
      # (at each interval, the folds are successively used as test data set
      # and the remaining folds as training data sets to compute the deviance)
      max.trees = 10000 # maximum number of trees that will be considered
      tolerance.method = "auto"
      tolerance = 0.001
      plot.main = TRUE
      plot.folds = FALSE
      verbose = TRUE
      silent = FALSE
      keep.fold.models = FALSE
      keep.fold.vector = FALSE
      keep.fold.fit = FALSE
      showingFoldsPlot = FALSE	
      n.folds = 5
      theRanges = c(current_block, current_block)*1000 # transform block size distance to meters
      
      centroids = as.data.frame(matrix(nrow=length(nutsM), ncol=3)); colnames(centroids) = c("NUTM","longitude","latitude")  
      centroids$NUTM = nutsM@data$LocationCo; coords = coordinates(nutsM)  
      centroids$longitude = coords[,1]; centroids$latitude = coords[,2]  
      
      matching_centroids = as.data.frame(matrix(nrow = nrow(current_data), ncol = 3)); colnames(matching_centroids) = c("NUTM","longitude", "latitude")
      matching_centroids$NUTM = current_data$nutM_code 
      matching_centroids$longitude = centroids$longitude[match(matching_centroids$NUTM, centroids$NUTM)]; matching_centroids$latitude = centroids$latitude[match(matching_centroids$NUTM, centroids$NUTM)]
      matching_centroids$response = current_data$response
      
      coordinates(matching_centroids) = ~longitude + latitude
      crs(matching_centroids) = crs(nutsM)
      
      spdf = SpatialPointsDataFrame(matching_centroids, current_data[,1:dim(current_data)[2]])
      myblocks = cv_spatial(spdf, column ="response", k=n.folds, size=theRanges, selection="random", plot = F)
      fold.vector = myblocks$folds_ids; n.trees = 100
      
      brt_model_scvs[[r]][[b]] =  gbm.step(current_data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
                                           var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
                                           verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit, plot.main = F)
      
      
      # 3. Model predictive performance evaluation (3 different metrics)
      ##################################################################
      
      #### Calculation of the AUC ####
      AUCs[b,r] = brt_model_scvs[[r]][[b]]$cv.statistics$discrimination.mean # Mean test AUC
      
      #### Calculation of the sorensen index to establish predictive performance of our models ####
      # Sources:
      # - computation performed according to the formulas of Leroi et al. (2018, J. Biogeography)
      # - optimisation of the threshold with a 0.01 step increment according to Li & Guo (2013, Ecography)
      
      tmp = matrix(nrow=101, ncol=2); tmp[,1] = seq(0,1,0.01)
      df = brt_model_scvs[[r]][[b]]$gbm.call$dataframe
      responses = df$response;
      data = df[,6:dim(df)[2]]
      n.trees = brt_model_scvs[[r]][[b]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
      prediction = predict.gbm(brt_model_scvs[[r]][[b]], data, n.trees, type, single.tree)
      N = dim(data)[1]; P = sum(responses==1); A = sum(responses==0)
      prev = P/(P+A) # proportion of recorded sites where the species is present
      x = (P/A)*((1-prev)/prev);
      sorensen_ppc = 0
      for (threshold in seq(0,1,0.01))
      {
        TP = length(which((responses==1)&(prediction>=threshold))) # true positives
        FN = length(which((responses==1)&(prediction<threshold))) # false negatives
        FP_pa = length(which((responses==0)&(prediction>=threshold))) # false positives
        sorensen_ppc_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
        tmp[which(tmp[,1]==threshold),2] = sorensen_ppc_tmp
        if (sorensen_ppc < sorensen_ppc_tmp)
        {
          sorensen_ppc = sorensen_ppc_tmp
          optimised_threshold = threshold
        }
      }
      tabs_list1[[r]][[b]] = tmp
      SIppcs[b,r] = sorensen_ppc
      thresholds[b,r] = optimised_threshold
      
     }
  }
  
  print(paste0("Saving results for ", pa.strategy[s], " sampling"))
  saveRDS(brt_model_scvs, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/brt_trained_model_", pa.strategy[s],".rds"))
  saveRDS(tabs_list1, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/tabs_list1_", pa.strategy[s],".rds"))
  write.csv(SIppcs, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/sippc_index_", pa.strategy[s],".csv"))
  write.csv(thresholds, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/sippc_thresholds_", pa.strategy[s],".csv"))
  write.csv(AUCs, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/AUCs_", pa.strategy[s], ".csv"))
}


# 4. Model performance figures 
##############################

# Figure X: SIppc curves
# load files 
ratios = c(1,3,6)
sorensen.data = list(random = readRDS(paste0(path, "European_model/All_the_brt_models/PA_strategies_brts/tabs_list1_random.rds")),
                     geo = readRDS(paste0(path,"European_model/All_the_brt_models/PA_strategies_brts/tabs_list1_geo_ae.density.rds")),
                     surv = readRDS(paste0(path,"European_model/All_the_brt_models/PA_strategies_brts/tabs_list1_surveillance.rds")),
                     es = readRDS(paste0(path,"European_model/All_the_brt_models/PA_strategies_brts/tabs_list1_espace.rds")))
# plot code
for (s in 1:length(sorensen.data)) {
  pdf(paste0("Figures/Supplementary/SIppc_curves", pa.strategy[s], ".pdf"), width=6, height=1.8)
  par(mfrow=c(1,3), oma=c(0,0.3,0,0), mar=c(2.5,2.5,0.5,0.5), lwd=0.4, col="gray30")
  for (i in 1:length(ratios)) {
    plot(sorensen.data[[s]][[i]][[1]], col=NA, ann=F, axes=F, xlim=c(0,1), ylim=c(0,1))
    for (j in 1:length(sorensen.data[[s]][[i]])) lines(sorensen.data[[s]][[i]][[j]], lwd=0.3, col="gray80", lty=1)
    axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.025, col.axis="gray30", mgp=c(0,0.14,0))
    axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.025, col.axis="gray30", mgp=c(0,0.35,0))
    title(ylab=expression("SI"["ppc"]), cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
    title(xlab="threshold", cex.lab=0.9, mgp=c(1.1,0,0), col.lab="gray30"); box(lwd=0.2, col="gray30")
    mtext(paste0("ratio ", namesRatios[[i]]), side=3, line=-1.3, at=0.995, cex=0.55, col="gray30", adj=1)
  }
  dev.off()
}


# Figure X : average model performance metrics table 
# Define statistical function
compute_stats = function(df, metric_name, ratios = c(1, 3, 6)) {
  summary_df = data.frame()
  
  for (i in seq_along(ratios)) {
    x = df[[i]]
    mean_val = mean(x, na.rm = TRUE)
    se = sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
    ci = 1.96 * se  # 95% CI
    
    CI_lower = mean_val - ci
    CI_upper = mean_val + ci
    
    formatted = sprintf("%.3f (%.3f – %.3f)", mean_val, CI_lower, CI_upper)
    
    summary_df = rbind(summary_df, data.frame(
      ratio = ratios[i],
      mean = mean_val,
      CI_lower = CI_lower,
      CI_upper = CI_upper,
      metric = metric_name,
      formatted = formatted
    ))
  }
  return(summary_df)
}

# Set file path
path = "European_model/All_the_brt_models/PA_strategies_brts/"

# Load all data files 
aucs = list(
  random = read.csv(paste0(path, "AUCs_random.csv"), header = TRUE),
  geo    = read.csv(paste0(path, "AUCs_geo_ae.density.csv"), header = TRUE),
  surv   = read.csv(paste0(path, "AUCs_surveillance.csv"), header = TRUE),
  es     = read.csv(paste0(path, "AUCs_espace.csv"), header = TRUE)
)

sippc = list(
  random = read.csv(paste0(path, "SIppcs_random.csv"), header = TRUE),
  geo    = read.csv(paste0(path, "SIppcs_geo_ae.density.csv"), header = TRUE),
  surv   = read.csv(paste0(path, "sippc_index_surveillance.csv"), header = TRUE),
  es     = read.csv(paste0(path, "sippc_index_espace.csv"), header = TRUE)
)

thresh = list(
  random = read.csv(paste0(path, "thresholds_random.csv"), header = TRUE),
  geo    = read.csv(paste0(path, "thresholds_geo_ae.density.csv"), header = TRUE),
  surv   = read.csv(paste0(path, "sippc_thresholds_surveillance.csv"), header = TRUE),
  es     = read.csv(paste0(path, "sippc_thresholds_espace.csv"), header = TRUE)
)

bi = list(
  random = read.csv(paste0(path, "boyce_index_random.csv"), header = TRUE),
  geo    = read.csv(paste0(path, "boyce_index_geo_ae.density.csv"), header = TRUE),
  surv   = read.csv(paste0(path, "boyce_index_surveillance.csv"), header = TRUE),
  es     = read.csv(paste0(path, "boyce_index_espace.csv"), header = TRUE)
)

# labels
pa_ratios = c("ratio 1:1", "ratio 1:3", "ratio 1:6")
strategies = c("random", "geo", "surv", "es")
strategy_names = c("Random sampling",
                   "Weigthed vector distribution sampling",
                   "Weighted surveillance capability sampling",
                   "Environmental space sampling")

# Initialize table
metrics.table = data.frame(
  `Pseudo-Absence (PA) sampling strategy` = rep(strategy_names, each = 3),
  `PA ratio` = rep(pa_ratios, times = 4),
  stringsAsFactors = FALSE
)

# Fill the table 
row_counter = 1
for (strategy in strategies) {
  auc_data = aucs[[strategy]][, -1]
  sippc_data = sippc[[strategy]][, -1]
  thresh_data = thresh[[strategy]][, -1]
  bi_data = bi[[strategy]][, -1]
  
  auc_stats = compute_stats(auc_data, "AUC")
  sippc_stats = compute_stats(sippc_data, "SIppc")
  thresh_stats = compute_stats(thresh_data, "Threshold")
  bi_stats = compute_stats(bi_data, "BI")
  
  for (i in 1:3) {
    metrics.table[row_counter, "Area Under the Curve (AUC)"] = auc_stats$formatted[i]
    metrics.table[row_counter, "Prevalence-pseudo-absence calibrated Sorensen index (SIppc)"] = sippc_stats$formatted[i]
    metrics.table[row_counter, "Threshold value maximising the SIppc"] = thresh_stats$formatted[i]
    metrics.table[row_counter, "Boyce Index (BI)"] = bi_stats$formatted[i]
    row_counter = row_counter + 1
  }
}

# write to Excel
write.xlsx(metrics.table, paste0(path,"summary_metrics_table.xlsx"), rowNames = FALSE)

################################################################################
#                       IV. SUPPLEMENTARY FIGURES
################################################################################

# 1. Density plots and pseudo-absence maps
##########################################
path = "Scripts_and_Data/all_data/all_the_brt_models/parallele/"

# Load files 
datas = list(
  random = readRDS(paste0(path, "env_training_data_random.rds")),
  geo    = readRDS(paste0(path, "env_training_data_geo_ae.density.rds")),
  surv   = readRDS(paste0(path, "env_training_data_surveillance.rds")),
  es     = readRDS(paste0(path, "env_training_data_espace.rds"))
)

rep_map = list()
for (strategy in names(datas)) {
  for (ratio in c(1, 3, 6)) {
    key = paste0(strategy, "_", ratio)
    rep_map[[key]] = sample(1:30, 1)
  }
}


# function to select the data and subset the absences
get_absences_sf = function(datas, strategy, ratio, nutsM, rep_index) {
  ratio_index = match(ratio, c(1, 3, 6))
  data = datas[[strategy]][[ratio_index]][[rep_index]]
  data = data[data$response == 0, ]
  abs_sf = nutsM[nutsM$LocationCo %in% data$nutM_code, ]
  return(st_as_sf(abs_sf))
}


# Plot the pseudo-absences
make_PA_plot = function(abs_sf) {
  ggplot(data = nutsM_sf) + 
    geom_sf(fill = "gray95", color = "gray50", size = 0.1) + 
    geom_sf(data = Ae.adj_geoms, fill = "snow3", color = "gray50", size = 0.2) +
    geom_sf(data = ae.pres_geoms, fill = "seashell4", color = "gray50", size = 0.1) +
    geom_sf(data = abs_sf, fill = "lightskyblue3", color = "gray50", size = 0.1) +
    theme_void()
}

# keep PA plots in list 
all_PA_plots = list()
for (strategy in names(datas)) {
  for (ratio in c(1, 3, 6)) {
    key = paste0(strategy, "_", ratio)
    rep_index = rep_map[[key]]
    
    abs_sf = get_absences_sf(datas, strategy, ratio, nutsM, rep_index)
    plot_name = key
    all_PA_plots[[plot_name]] = make_PA_plot(abs_sf)
  }
}

# Load data and calculate density
get_density_data = function(datas, strategy, ratio, rep_index) {
  ratio_index = match(ratio, c(1, 3, 6))
  data = datas[[strategy]][[ratio_index]][[rep_index]]
  data_density = data[, c("nutM_code", "response", "temperature_summer")]
  data_density$response = factor(data_density$response, levels = c(0, 1), labels = c("Absence", "Presence"))
  
  presence = data_density$temperature_summer[data_density$response == "Presence"]
  absence = data_density$temperature_summer[data_density$response == "Absence"]
  ov = round(overlap(list(presence, absence))$OV, 2) * 100
  
  list(data = data_density, overlap = ov)
}

# Plot density plots
make_density_plot = function(data_density, overlap_value) {
  ggplot(data_density, aes(x = temperature_summer, fill = response)) + 
    geom_density(alpha = 0.5, linewidth = 0.2, col = "gray20") +  
    scale_fill_manual(values = c("Absence" = "lightskyblue3", "Presence" = "seashell4")) +
    labs(
      x = paste0("Summer temperature – overlap area ", overlap_value, "%"),
      y = "Density",
      fill = "Category"
    ) +
    theme_light() +
    guides(fill = "none")
}


# Store plots in list
density_plots = list()

for (strategy in names(datas)) {
  for (ratio in c(1, 3, 6)) {
    key = paste0(strategy, "_", ratio)
    rep_index = rep_map[[key]]
    
    dens_result = get_density_data(datas, strategy, ratio, rep_index)
    density_plots[[key]] = make_density_plot(dens_result$data, dens_result$overlap)
  }
}

## add pseudo absence plot and density plot together
for (strategy in names(datas)) {
  pdf(file = paste0(path, "PA_densities_", strategy, ".pdf"), width = 8, height = 6)
  
  grid.arrange(
    all_PA_plots[[paste0(strategy, "_1")]],
    all_PA_plots[[paste0(strategy, "_3")]],
    all_PA_plots[[paste0(strategy, "_6")]],
    density_plots[[paste0(strategy, "_1")]],
    density_plots[[paste0(strategy, "_3")]],
    density_plots[[paste0(strategy, "_6")]],
    nrow = 2, ncol = 3,
    top = paste("Strategy:", strategy)
  )
  dev.off()
}


# 2. Current risk maps
#######################

# 1. Extract the environmental data corresponding to t0 (current decades 2000-2024)
###################################################################################

envVariableNames = c("temperature_winter","temperature_spring","temperature_summer",
                     "temperature_fall", "precipitation_winter", "precipitation_spring", "precipitation_summer",
                     "precipitation_fall", "relative_humidity_winter", "relative_humidity_spring", "relative_humidity_summer",
                     "relative_humidity_fall", "croplands", "pastures","rangelands", "urbanAreas", "primaryForest","primaryNonForest",
                     "secondaryForest","secondaryNonForest", "population")


croplands = brick("../BCLIMATE/output/landuse/LUH2-GCB2024/croplands1901_2024.nc")
pastures = brick("../BCLIMATE/output/landuse/LUH2-GCB2024/pastures1901_2024.nc")
rangelands = brick("../BCLIMATE/output/landuse/LUH2-GCB2024/rangeland1901_2024.nc")
urbanAreas = brick("../BCLIMATE/output/landuse/LUH2-GCB2024/urbanAreas1901_2024.nc")
primaryForest = brick("../BCLIMATE/output/landuse/LUH2-GCB2024/primaryForest1901_2024.nc")
primaryNonForest = brick("../BCLIMATE/output/landuse/LUH2-GCB2024/primaryNonF1901_2024.nc")
secondaryForest = brick("../BCLIMATE/output/landuse/LUH2-GCB2024/secondaryForest1901_2024.nc")
secondaryNonForest = brick("../BCLIMATE/output/landuse/LUH2-GCB2024/secondaryNonF1901_2024.nc")
population = brick("../BCLIMATE/output/population/ISIMIP3a/population_histsoc_30arcmin_annual_1850_2021.nc", varname = "total-population")
start_year = 1901; end_year = 2021; start_layer = start_year - 1850 + 1; end_layer = end_year - 1850 + 1
population = population[[start_layer:end_layer]]

temperature = brick("../BCLIMATE/output/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc")
precipitation = brick("../BCLIMATE/output/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_pr_1901_2021_monmean.nc")
relativehumidity = brick("../BCLIMATE/output/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_hurs_1901_2021_monmean.nc")


year_start = 2000 ; year_end = 2024 ; year_end2 = 2021

month_start = ((year_start-1) - 1901)*12+1; month_end = (year_end2 - 1901)*12 + 12
temperature_temp = temperature[[month_start:month_end]]
precipitation_temp = precipitation[[month_start:month_end]]
relativehumidity_temp = relativehumidity[[month_start:month_end]]

winter_indices = sort(c(12, 1, 2)) + rep((year_start:year_end2 - year_start) * 12, each = 3);  winter_indices = winter_indices[-c(1, 2, length(winter_indices))]
spring_indices = sort(c(3,4,5) + rep(((year_start+1):year_end2 - year_start) * 12, each = 3))
summer_indices = sort(c(6,7,8) + rep(((year_start+1):year_end2 - year_start) * 12, each = 3))
fall_indices = sort(c(9,10,11) + rep(((year_start+1):year_end2 - year_start) * 12, each = 3))


tempm_winter =  mean(temperature_temp[[winter_indices]]) - 273.15 # conversion to Celsius
tempm_spring = mean(temperature_temp[[spring_indices]]) - 273.15
tempm_summer = mean(temperature_temp[[summer_indices]]) - 273.15
tempm_fall = mean(temperature_temp[[fall_indices]]) - 273.15

precp_winter = mean(precipitation_temp[[winter_indices]]) * 60 * 60 * 24 # conversion to kg/m2/day
precp_spring = mean(precipitation_temp[[spring_indices]]) * 60 * 60 * 24
precp_summer = mean(precipitation_temp[[summer_indices]]) * 60 * 60 * 24
precp_fall = mean(precipitation_temp[[fall_indices]]) * 60 * 60 * 24

relh_winter = mean(relativehumidity_temp[[winter_indices]])
relh_spring = mean(relativehumidity_temp[[spring_indices]]) 
relh_summer = mean(relativehumidity_temp[[summer_indices]])
relh_fall = mean(relativehumidity_temp[[fall_indices]])


year_start_temp = year_start - 1900;  year_end_temp = year_end - 1900;year_end_temp2 = year_end2 - 1900
croplands_temp = mean(croplands[[year_start_temp:year_end_temp]])
pastures_temp = mean(pastures[[year_start_temp:year_end_temp]])
urbanAreas_temp = mean(urbanAreas[[year_start_temp:year_end_temp]])
rangelands_temp = mean(rangelands[[year_start_temp:year_end_temp]])
primaryForest_temp = mean(primaryForest[[year_start_temp:year_end_temp]])
primaryNonForest_temp = mean(primaryNonForest[[year_start_temp:year_end_temp]])
secondaryForest_temp = mean(secondaryForest[[year_start_temp:year_end_temp]])
secondaryNonForest_temp = mean(secondaryNonForest[[year_start_temp:year_end_temp]])
population_temp = mean(population[[year_start_temp:year_end_temp2]])

envVariables = list()
envVariables[[1]] = tempm_winter
envVariables[[2]] = tempm_spring
envVariables[[3]] = tempm_summer
envVariables[[4]] = tempm_fall

envVariables[[5]] = precp_winter
envVariables[[6]] = precp_spring
envVariables[[7]] = precp_summer
envVariables[[8]] = precp_fall

envVariables[[9]] =  relh_winter
envVariables[[10]] =  relh_spring
envVariables[[11]] =  relh_summer
envVariables[[12]] =  relh_fall

envVariables[[13]] =  croplands_temp 
envVariables[[14]] =  pastures_temp 
envVariables[[15]] =  urbanAreas_temp
envVariables[[16]] =  rangelands_temp
envVariables[[17]] =  primaryForest_temp 
envVariables[[18]] =  primaryNonForest_temp
envVariables[[19]] =  secondaryForest_temp 
envVariables[[20]] =  secondaryNonForest_temp
envVariables[[21]] =  population_temp

for (i in 1:length(envVariables)) envVariables[[i]] = crop(envVariables[[i]], nutsM, snap="out")
for (i in 1:length(envVariables)) envVariables[[i]] = mask(envVariables[[i]], nutsM)

nutsM_data = data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = length(envVariableNames) + 1))
colnames(nutsM_data) = c("NUTM", envVariableNames)
nutsM_data$NUTM = nutsM_IDs$nutM_code

for (i in seq_along(envVariables)) { 
  if (i == length(envVariables)) {
    extracted_values = log(exactextractr::exact_extract(envVariables[[i]], nutsM, fun = "sum") + 1)
  } else {
    extracted_values = exact_extract(envVariables[[i]], nutsM, fun = "mean")
  } 
  nutsM_data[, i + 1] = extracted_values  # + 1 because first column is NUTM
}

#NA_ids = c("EL304","EL307","EL413","EL421","EL422","EL624","ES531","ES533",
#"ES630","ES640","GG","GI","IM","JE","MEG125354","MEG125368","UKM65")

#nutsM_data = nutsM_data[!nutsM_data$NUTM %in% NA_ids, ] # remove islands

#2. Predict the current risk map using all the trained BRT models
#################################################################

# load brts
path = "Scripts_and_Data/all_data/all_the_brt_models/parallele/"
random.brt = readRDS(paste0(path,"brt_trained_model_random.rds"))
geo.brt  = readRDS(paste0(path,"brt_trained_model_geo_ae.density.rds"))
surv.brt = readRDS(paste0(path,"brt_trained_model_surveillance.rds"))                      
es.brt = readRDS(paste0(path,"brt_trained_model_espace.rds"))  

# Current predictions 
random.preds = list(); geo.preds = list(); surv.preds = list(); es.preds = list() 
for(i in 1:length(random.brt)){
  random.preds[[i]] = list()
  geo.preds[[i]] = list()
  surv.preds[[i]] =list()
  es.preds[[i]] = list()
  for(j in 1:length(random.brt[[i]])){
    
    object = random.brt[[i]][[j]]; n.trees = random.brt[[i]][[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
    object2 = geo.brt[[i]][[j]]; n.trees2 = geo.brt[[i]][[j]]$gbm.call$best.trees; type2 = "response"; single.tree2 = FALSE
    object3 = surv.brt[[i]][[j]]; n.trees3 = surv.brt[[i]][[j]]$gbm.call$best.trees; type3 = "response"; single.tree3 = FALSE
    object4 = es.brt[[i]][[j]]; n.trees4 = es.brt[[i]][[j]]$gbm.call$best.trees; type4 = "response"; single.tree4 = FALSE
    
    # Predicting ecological suitability using the trained models
    prediction = predict.gbm(object, nutsM_data[, -1], n.trees, type, single.tree)
    prediction2 = predict.gbm(object2, nutsM_data[, -1], n.trees2, type2, single.tree2)
    prediction3 = predict.gbm(object3, nutsM_data[, -1], n.trees3, type3, single.tree3)
    prediction4 = predict.gbm(object4, nutsM_data[, -1], n.trees4, type4, single.tree4)
    
    random.preds[[i]][[j]] = prediction
    geo.preds[[i]][[j]] = prediction2
    surv.preds[[i]][[j]] = prediction3
    es.preds[[i]][[j]] = prediction4
  }
}


# Group all predictions 
dfs_random = list(); dfs_geo = list(); dfs_surv = list(); dfs_es = list()
for (i in seq_along(random.preds)) {
  data.random = as.data.frame(matrix(nrow = 949, ncol = 31))
  data.geo = as.data.frame(matrix(nrow = 949, ncol = 31))
  data.surv = as.data.frame(matrix(nrow = 949, ncol = 31))
  data.es = as.data.frame(matrix(nrow = 949, ncol = 31))
  for (j in 1:30) {
    data.random[, j] = random.preds[[i]][[j]]
    data.geo[, j] = geo.preds[[i]][[j]]
    data.surv[, j] = surv.preds[[i]][[j]]
    data.es[, j] = es.preds[[i]][[j]]
  }
  
  data.random[, 31] = rowMeans(data.random[, 1:30])
  data.geo[, 31] = rowMeans(data.geo[, 1:30])
  data.surv[, 31] = rowMeans(data.surv[, 1:30])
  data.es[, 31] = rowMeans(data.es[, 1:30])
  
  dfs_random[[i]] = data.random
  dfs_geo[[i]] = data.geo
  dfs_surv[[i]] = data.surv
  dfs_es[[i]] = data.es
}



# plot Current Risk maps
nutsM_sf = st_as_sf(nutsM)
contour = st_union(nutsM_sf)
colourScale = rev(colorRampPalette(brewer.pal(11, "RdBu"))(121)[11:111])

legend1 = raster(as.matrix(c(0,1)))
models = c("Random PA sampling", "Ae. albopictus suitability \nPA sampling", "Surveillance capability\n PA sampling","Environmental space \nPA sampling")
ratios = c("1.1", "1.3", "1.6")
preds = list(dfs_random, dfs_geo, dfs_surv, dfs_es)


pdf("Scripts_and_Data/all_data/all_the_brt_models/parallele/current_predictions_aedes2.pdf", width=6, height=6)
#dev.new(width=6, height=6)
par(mfrow=c(4,3), oma=c(0,0.1,0.5,1.5), mar=c(0,2,0,0), lwd=0.1, col="gray30")
for (i in 1:length(preds)) {
  for (j in 1:length(ratios)) {
    cols = colourScale[(((preds[[i]][[j]][, 31] - 0) / (1 - 0)) * 100) + 1]
    
    plot(contour, lwd=0.4, border="gray30", col=NA)
    plot(nutsM, col = cols, border = NA, lwd = 0.1, ann = FALSE, legend = FALSE, axes = FALSE, box = FALSE, add=T)
    mtext("", side=3, line=0.3, cex=0.65, col="gray30");
    if (j == 1) {text(x = par("usr")[1] - 0.05 * diff(par("usr")[1:2]), y = mean(par("usr")[3:4]), labels = models[i], xpd = NA, srt = 90, cex = 0.85, col="gray30") }
    if (i==1){mtext(paste0("Ratio ",ratios[[j]]), side=3, line=-0.7, cex=0.65, col="gray30")}else{mtext = ""}
    if (i==4 && j ==3){plot(legend1, col= colourScale, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
                            smallplot=c(0.93,0.96,0.15,0.85), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2,
                                                                                    col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0, mgp=c(2,0.5,0), 
                                                                                    at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")), alpha=1, side=3)}
  }
}
dev.off()










