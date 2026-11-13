##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#### 00e: Process NLCD data for Tree Canopy and Impervious Surface Coverage ####

#coded by: Aja Sutton
#Climate Change Epidemiology Lab
# Scripps Institution of Oceanography, UC San Diego

#load packages; if on CRAN and not installed, install and load.
{pacman::p_load(tidyverse,sf,janitor,here,terra,tigris)#,FedData)
  
  # Automatically check, download locationallocation package from GitHub if missing, and load it
  Require::Require("giacfalk/locationallocation")
  
  options(scipen=999)
  options(tigris_use_cache = TRUE)
  
  tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)
  
  here::i_am("code/00e_get_nlcd_tree_canopy_impervious_surfaces.R")
}

#Goal: identify *where* cooling centers should be located in a CA context
#First, predict where they ought to be based on a weight we construct
#Then, analyze how these compare to where locations used as cooling centers already exist

#Meta-goal: Make this a process practitioners will be able to replicate 
#simply and easily without advanced technical skills beyond basic R skills.
#IE *highly* reproducible and fairly easy to use!

#To undertake location-allocation, we must create a weighted spatial field
#to prioritize the characteristics we think matter when deciding where cooling centres
#should be located (vs where they are, even!).


#Step 00: Gather data to produce a weighted spatial field for the locationallocation call

#Assuming these are all equally weighted to create a spatial field weighted by need
#They must all be scaled to themselves, and then summed to be meaningful.

#Data needed:
# AC coverage data (census tract-level): https://www.census.gov/data/experimental-data-products/lace.html
# Population density: US Census Bureau
# CDC Social Vulnerability Index: https://www.atsdr.cdc.gov/place-health/php/svi/svi-data-documentation-download.html
# EVI: VIIRS (not MODIS) (call API, see below)
# Tree Canopy: CALFIRE: https://www.fs.usda.gov/r05/state-private-tribal/california-urban-canopy-data
#  
# Extreme Heat data from PRISM: https://prism.oregonstate.edu/
# For now we use the 30-year annual t_max (temperature maximum) normals (average) 
# for CA at 800m resolution: https://prism.oregonstate.edu/normals/
# This won't account for anomalous extremes, but we can plan for places that are typically extreme
# Warm months vs off-season? 
# The current data choice significantly reduces analytic complexity while still being useful for prediction

#### Focus Area ####

# place_name = "San Diego city"
# place = tigris::places(state = "CA", cb = TRUE)
# place = place %>% 
#   clean_names() %>%
#   filter(str_detect(namelsad,place_name))

state = tigris::states()
state = 
  state %>%
  clean_names() %>%
  filter(name == "California") %>%
  select(name,geometry)

#### NLCD Tree Canopy and Impervious Surfaces Data ####

#Easy to download national data here: https://www.mrlc.gov/data
#Fractional Impervious Surface (CONUS) 2025 &
#NLCD Tree Canopy Cover (CONUS) 2025

#Tree Canopy
dat = rast(here('data/nlcd_tcc_conus_2025_v2025-6_wgs84/nlcd_tcc_conus_wgs84_v2025-6_20250101_20251231.tif'))

state = st_transform(state,st_crs(dat))

#remove special values: 
#254 (Non-processing area, e.g., persistent cloud/shadow)
#255 (Background / NoData)
dat = subst(dat, from = c(254,255), to = NA)
dat = crop(dat,state)
dat = mask(dat,state)

writeRaster(dat, file=here("data/clean/2025_CA_nlcd_tree_canopy.tif"), datatype = "INT2U", overwrite = TRUE)

rm(dat)
gc()

#Impervious Surfaces
dat = rast(here('data/Annual_NLCD_FctImp_2025_CU_C1V2/Annual_NLCD_FctImp_2025_CU_C1V2.tif'))

#remove special values: 250 is NA no data
dat = subst(dat, from = 250, to = NA)
dat = crop(dat,state)
dat = mask(dat,state)

writeRaster(dat, file=here("data/clean/2025_CA_nlcd_fct_impervious.tif"), overwrite = TRUE)

