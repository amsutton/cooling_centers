##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#### 01: Load Data, Filter to Place of Interest ####

#coded by: Aja Sutton
#Climate Change Epidemiology Lab
# Scripps Institution of Oceanography, UC San Diego

#load packages; if on CRAN and not installed, install and load.
{pacman::p_load(tidyverse,sf,janitor,here,tidycensus,terra,Require,viridis,
                tigris,crsuggest,data.table,exactextractr)#,FedData)
  
  # Automatically check, download locationallocation package from GitHub if missing, and load it
  Require::Require("giacfalk/locationallocation")
  
  options(scipen=999)
  options(tigris_use_cache = TRUE)
  
  tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)
  
  here::i_am("code/01_filter_poi_data.R")


#Goal: identify *where* cooling centers should be located in a CA context
#First, predict where they ought to be based on a weight we construct
#Then, analyze how these compare to where locations used as cooling centers already exist

#Meta-goal: Make this a process practitioners will be able to replicate 
#simply and easily without advanced technical skills beyond basic R skills.
#IE *highly* reproducible and fairly easy to use!

#To undertake location-allocation, we must create a weighted spatial field
#to prioritize the characteristics we think matter when deciding where cooling centres
#should be located (vs where they are, even!).

#specify place of interest (poi)
place_name = "San Diego city"
poi = tigris::places(state = "CA", cb = TRUE)
poi = poi %>% 
  clean_names() %>%
  filter(str_detect(namelsad,place_name))

crs_choice = 
  suggest_crs(poi) %>% 
  filter(crs_units == "m") %>% 
  slice(1)
crs_choice = as.integer(crs_choice$crs_code)
}

#load and normalize data

#### impervious surfaces ####
dat = rast(here('data/clean/2025_CA_nlcd_fct_impervious.tif'))
poi = st_transform(poi,st_crs(dat))
dat = crop(dat,poi)
dat = mask(dat,poi)

#raw on 100 point scale
writeRaster(dat, here("data/clean/poi/san_diego/2025_sandiego_nlcd_fct_impervious.tif"),datatype = "INT1U", overwrite = TRUE)

#programmatic quantiles
breaks = as.numeric(global(dat, fun=quantile, probs=seq(0, 1, 0.25), na.rm=TRUE))
test = classify(dat, rcl=breaks, include.lowest=TRUE)
writeRaster(test, here("data/clean/poi/san_diego/2025_sandiego_nlcd_fct_impervious_quantile.tif"),datatype = "INT1U", overwrite = TRUE)

#log
dat = subst(dat, 0, 0.1)
dat = log(dat)
writeRaster(dat, here("data/clean/poi/san_diego/2025_sandiego_nlcd_fct_impervious_log.tif"),overwrite = TRUE)

rm(dat,test)
gc()



#### tree canopy ####
dat = rast(here('data/clean/2025_CA_nlcd_tree_canopy.tif'))
poi = st_transform(poi,st_crs(dat))
dat = crop(dat,poi)
dat = mask(dat,poi)

#raw on 100 point scale
writeRaster(dat, here("data/clean/poi/san_diego/2025_sandiego_nlcd_tree_canopy.tif"),datatype = "INT1U", overwrite = TRUE)

# svi$rplthemes = 
#   cut(svi$rpl_themes, 
#       breaks = c(0, 0.25, 0.5, 0.75, 1), 
#       include.lowest = TRUE)


#programmatic quantiles
breaks = as.numeric(global(dat, fun=quantile, probs=seq(0, 1, 0.25), na.rm=TRUE))
test = classify(dat, rcl=breaks, include.lowest=TRUE)
writeRaster(test, here("data/clean/poi/san_diego/2025_sandiego_nlcd_tree_canopy_quantile.tif"),datatype = "INT1U", overwrite = TRUE)

#log
dat = subst(dat, 0, 0.1)
dat = log(dat)
writeRaster(dat, here("data/clean/poi/san_diego/2025_sandiego_nlcd_tree_canopy_log.tif"),overwrite = TRUE)

rm(dat,test)
gc()



#### EVI ####
#Already focused on poi -- must load by city for now.
#this is for san diego:
#evi = rast("data/clean/poi/san_diego/mean_evi_october_may_2024.tif")



#### CDC Social Vulnerability Index 2022 ####
svi = st_read(here('data/clean/2022_cdc_svi_ca_tracts_rpl_themes_clean.gdb/'))

svi = st_transform(svi,crs_choice)
poi = st_transform(poi,crs_choice)

svi = st_intersection(svi,poi)
svi = st_collection_extract(svi, type = "POLYGON")

svi %>%
  #mutate(rpl_themes = quantile(rpl_themes,probs=seq(0,1,0.25),na.rm=TRUE)) %>%
  ggplot(aes(fill =  rpl_themes)) +
    geom_sf() +
    scale_fill_viridis(direction = -1)

ggsave(file=here('data/clean/poi/san_diego/cdc_svi_sandiego_rplthemes.png'),width=10,height=14, unit="in",dpi=300)

st_write(svi,
         here('data/clean/poi/san_diego/2022_cdc_svi_san_diego.gdb'),
         append = FALSE)

svi$rplthemes = 
  cut(svi$rpl_themes, 
  breaks = c(0, 0.25, 0.5, 0.75, 1), 
  include.lowest = TRUE,
  labels = c(1,2,3,4))

svi %>%
  ggplot(aes(fill =  rplthemes)) +
  geom_sf() +
  scale_fill_viridis(direction = -1,discrete = TRUE)

ggsave(file=here('data/clean/poi/san_diego/cdc_svi_sandiego_rplthemes_breaks.png'),width=10,height=14, unit="in",dpi=300)

st_write(svi,
         here('data/clean/poi/san_diego/2022_cdc_svi_san_diego.gdb'),
         layer = "2022_cdc_svi_san_diego_breaks",
         append = FALSE)

#### LACE Air Conditioning Coverage 2023 ####
#fine as is


#### PRISM heat index ####

#quantiles
