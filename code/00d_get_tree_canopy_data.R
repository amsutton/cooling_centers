##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#coded by: Aja Sutton
#Climate Change Epidemiology Lab
# Scripps Institution of Oceanography, UC San Diego

#load packages; if on CRAN and not installed, install and load.
{pacman::p_load(tidyverse,sf,janitor,here,tidycensus,terra,Require,
                prism,tigris,crsuggest,data.table,exactextractr,heatindex)#,FedData)
  
  # Automatically check, download locationallocation package from GitHub if missing, and load it
  Require::Require("giacfalk/locationallocation")
  
  options(scipen=999)
  options(tigris_use_cache = TRUE)
  
  tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)
  
  here::i_am("code/00d_get_tree_canopy_data.R")
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

place_name = "San Diego city"

#### CALFIRE CA Urban Tree Canopy Data, 2018 & 2022 ####

# #most recent tree canopy data from 2022  
# trees = rast(here('data/CALFIRE_UrbanTreeCanopy_2018_2022/CALFIRE_UrbanTreeCanopy_2018_2022.gdb/'))
# trees = trees$CA_UrbanTreeCanopy_2022
# 
# #to make the cropping and masking easier and faster:
# coj = st_transform(coj,crs=st_crs(trees))
# trees = crop(trees,coj)
#trees = mask(trees,coj)
#writeRaster(trees,here('data/CALFIRE_UrbanTreeCanopy_2018_2022/urbantreecanopy_2022_clipped.tif'),overwrite=TRUE)

#reproject to match other objects
trees = rast(here('data/CALFIRE_UrbanTreeCanopy_2018_2022/urbantreecanopy_2022_clipped.tif'))

#fastest way to summarize mean in tracts for tree canopy; 
#works because canopy consists of zeroes and ones
#also handles fractional weighting!
#make sure they match crs first

tracts = tigris::tracts(state="CA")
tracts = sf::st_intersection(tracts, coj)

tracts = st_transform(tracts,st_crs(trees))
tracts = tracts %>% clean_names() %>% select(geoid,geometry) %>% distinct()

trees_mean = exact_extract(trees, tracts, c("mean","median","min","max","stdev","sum"))

trees = cbind(trees_mean,tracts)

trees = trees %>%
  select(geoid,mean:sum,geometry) %>%
  st_as_sf()

trees = st_transform(trees,crs=crs_choice)

saveRDS(trees,file= here('data/clean/sandiego_trees_by_tract.RDS'))

#Development Density Data
#for now, 2021 because later years won't call properly
# nlcd <- get_nlcd(
#     template = as(coj,"Spatial"),
#     year = 2021,
#     label = paste0(place_name)
#   )
#   
# glimpse(nlcd)

#temp = as.data.frame(nlcd)
#table(temp$Class) #to get a sense of comparison
# temp %>%
#   group_by(Class) %>%
#   summarise(n = n()) %>%
#   ungroup() %>%
#   ggplot() +
#   geom_point(aes(y=Class,x=log(n)))

