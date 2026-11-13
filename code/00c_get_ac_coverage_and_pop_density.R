##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#coded by: Aja Sutton
#Climate Change Epidemiology Lab
# Scripps Institution of Oceanography, UC San Diego

#load packages; if on CRAN and not installed, install and load.
{pacman::p_load(tidyverse,sf,janitor,here,tidycensus,terra,Require,
                tigris,crsuggest,data.table,exactextractr)#,FedData)
  
  # Automatically check, download locationallocation package from GitHub if missing, and load it
  Require::Require("giacfalk/locationallocation")
  
  options(scipen=999)
  options(tigris_use_cache = TRUE)
  
  tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)
  
  here::i_am("code/00c_get_ac_coverage_and_pop_density.R")
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
poi = tigris::places(state = "CA", cb = TRUE)
poi = poi %>% 
  clean_names() %>%
  filter(str_detect(namelsad,place_name))


#### Get AC Coverage Data from LACE 2023 ####
# national data at the census tract level, but also available at zcta if preferred later
ac = fread(here('data/local_ac_est_23_Tract.csv'))

#margin of error at 90% available, but we take point estimates for now
#filter out tracts that are 100% water
ac = ac %>%
  clean_names() %>%
  filter(str_detect(name,"California"),
         is.na(water_tract)) %>%
  rename(geoid = geo_id) %>%
  select(geoid,
         tract,
         name, #combination of census tract, county, state
         hse_occ_e, #HSE_OCC_E: Estimate for the total number of occupied households.
         ac_e, #AC_E: Estimate for the number of occupied households with any kind of air conditioning.
         ac_pe, #AC_PE: Estimate for the percentage of occupied households with any kind of air conditioning.
         no_ac_e, #without
         no_ac_pe #without
  )

#filter to area of interest
tracts = tigris::tracts(state="CA",year=2023)
tracts = tracts %>%
  clean_names() %>%
  select(geoid,geoidfq,tractce,namelsad,aland,awater,geometry)


#AC by tracts
ac = left_join(ac,tracts, by = c("geoid"="geoidfq"))
ac = st_as_sf(ac)
ac = sf::st_intersection(ac, poi)
ac = ac %>%
  select(geoid,name,aland,awater,state_name,namelsad.1,
         no_ac_pe:hse_occ_e,geometry) %>%
  rename("place"="namelsad.1")

saveRDS(ac,file=here("data/clean/ac_coverage_sandiego.rds"))

#### Get Population Density (Block Group Level) ####
popdensity = get_acs(geography = "block group",
                     table = "B01001",
                     state= "CA",
                     year = 2023,
                     survey = "acs5",
                     geometry = TRUE,
                     labels=TRUE)


#lots of pop characteristics in here by age and sex, but also total population!

#density requires land area
popdensity$larea = st_area(popdensity$geometry)/1000000#sq metres to sq kms

#takes about 10 minutes to run, careful:
popdensity = 
  popdensity %>%
  clean_names() %>%
  group_by(variable, geoid) %>%
  mutate(popdensity = estimate/larea) %>% 
  select(geoid,
         name,
         variable,
         estimate,
         popdensity,
         geometry) %>% 
  arrange(desc(popdensity)) #for fun 


saveRDS(popdensity,file=here('data/clean/ca_population_density_bysexage_5yr_acs_2023.rds'))

