##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#coded by: Aja Sutton
#Climate Change Epidemiology Lab
# Scripps Institution of Oceanography, UC San Diego

#load packages; if on CRAN and not installed, install and load.
{pacman::p_load(tidyverse,sf,janitor,here,tidycensus,terra,Require,viridis,
                prism,tigris,crsuggest,data.table,exactextractr,heatindex)#,FedData)
  
  # Automatically check, download locationallocation package from GitHub if missing, and load it
  Require::Require("giacfalk/locationallocation")
  
  options(scipen=999)
  options(tigris_use_cache = TRUE)
  
  tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)
  
  here::i_am("code/00f_get_cdc_svi2022_cchvi_cdph20172019.R")
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

#### Get and clean CDC's Social Vulernability Index data ####

svi = st_read(here('data/SVI2022_CALIFORNIA_tract.gdb'))
svi = svi %>%  
  clean_names() %>%
  mutate(rpl_themes = ifelse(rpl_themes == -999, NA, rpl_themes))

svi %>%
  #mutate(norm_rplthemes = rpl_themes/sum(rpl_themes)) %>%
  filter(county== "San Diego County") %>%
  ggplot(aes(fill =  rpl_themes)) +
  geom_sf() +
  scale_fill_viridis(direction = -1)
ggsave(file=here('data/clean/viz/cdc_svi_sandiego_county_rplthemes.png'),width=10,height=14, unit="in",dpi=300)

st_write(svi,here('data/clean/2022_cdc_svi_ca_tracts_rpl_themes_clean.gdb'))



#### CCHVI from CDPH 2017-2019 ####
dat = fread('/Users/Aja/Documents/R Directories/cooling_centers/data/ca_climate_change_health_vulnerability_indicators_cdph2026.csv')

tracts = tigris::tracts(state="CA")
tracts = sf::st_intersection(tracts, poi)

dat = dat %>% clean_names() %>% filter(county=="San Diego")

tracts = tracts %>%
  clean_names() %>%
  mutate(geoid_new = str_remove_all(geoidfq,"1400000US0"))

dat$fips = as.character(dat$fips)
dat = left_join(dat,tracts, by = c("fips"="geoid_new"))

table(dat$definition)

#indicators of potential interest
indicators = c("Percent of population aged 65 years or older",
               "Percent impervious surface cover",
               "Percent of households with no vehicle ownership",
               "Percent of households without air conditioning",
               "Percent of population age less than 5 years",
               "Percent of population aged 65 years or older",
               "Percent of population with a disability",
               "Percent of population employed and aged > 16 working outdoors",
               "Percent of population without health insurance",
               "Poverty Rate (200% FPL)",
               "Projected number of extreme heat days 2040-2060",
               "Projected number of extreme heat days 2080-2099")

dat %>%
  st_as_sf() %>%
  filter(definition %in% indicators) %>%
  st_simplify() %>% 
  ggplot(aes(fill = mean)) +
  geom_sf() +
  facet_wrap(~definition) +
  scale_fill_viridis()

ggsave(file=here('data/clean/viz/cchvi_sd_indicator_summary.png'),width=10,height=14, unit="in",dpi=300)


