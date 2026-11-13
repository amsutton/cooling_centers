##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#coded by: Aja Sutton
#Climate Change Epidemiology Lab
# Scripps Institution of Oceanography, UC San Diego

#load packages; if on CRAN and not installed, install and load.
{pacman::p_load(tidyverse,sf,janitor,here,tidycensus,terra,Require,
                prism,tigris,crsuggest,data.table,exactextractr,MODISTools)#,FedData)
  
  # Automatically check, download locationallocation package from GitHub if missing, and load it
  Require::Require("giacfalk/locationallocation")
  
  options(scipen=999)
  options(tigris_use_cache = TRUE)
  
  tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)
  
  here::i_am("code/00b_get_evi_data.R")
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




#### Get EVI data ####
#see: https://jdmwhite.github.io/Intro_to_Spatial/NDVI_time_series.html
#Using MODISTools package
#code taken from: https://github.com/mheskel/MODISTools-Lab-368/blob/master/Week%2013%20Lab%20-%20MODIS%20NDVI%20and%20EVI.R

#get site centroid, for now only for San Diego
place = tigris::places(state = "CA", cb = TRUE)
place = place %>% 
  clean_names() %>%
  filter(str_detect(namelsad,place_name))

crs_choice = suggest_crs(place) %>% filter(crs_units == "m") %>% slice(1)
crs_choice = as.integer(crs_choice$crs_code)

#get max length from centroid for EVI extraction extent in kilometres
#builds a sort of bounding box around the centroid of the place to make sure
#extent of place of interest is covered in EVI data
length = st_cast(place,"POINT")
poi = place
place = st_centroid(place)

#right length this way
length = st_distance(st_transform(length,crs = crs_choice),st_transform(place,crs_choice))
length = as.double(max(length)/1000)

#VIIRS is current product, MODIS is deprecated soon.
products <- mt_products()
products %>% filter(str_detect(description,"EVI"))

#This shows the list of available MODIS products that are extractable. 
# These include land cover type, vegatation indices EVI and NDVI, and others. 

##To check the dates available for NDVI/EVI at a certain lat/lon:
# dates<-mt_dates(product="VNP13A1", lat = st_coordinates(place)[2], lon=st_coordinates(place)[1]) #fill in lat and lon for your site
# print(dates)
# 
# lubridate::year(dates$calendar_date) %>% unique()
# 
# mt_bands(product= "VNP13A1")

#Download 16-day average EVI for most recent year:
#May to October, band = "500_m_16_days_EVI"
{
  start_time = Sys.time()
  evi = mt_subset(product = "VNP13A1",
                  lat = st_coordinates(place)[2],
                  lon = st_coordinates(place)[1], 
                  start = "2024-05-01",
                  end = "2024-10-31", 
                  km_lr = length,
                  km_ab = length,
                  band = "500_m_16_days_EVI",
                  site_name=paste0(place_name),
                  internal = TRUE,
                  progress = TRUE)
  end_time = Sys.time()
  end_time - start_time
}

#For now, we are probably most interested in extreme heat during warm months (May-October),
#and how well EVI may perform at mitigating extreme heat during these months.
#Create area-level summaries at the census tract level.


##it's very barren in SD in the summertime, no surprise there!
# mean_evi = evi %>%
#   mutate(doy = yday(as_date(calendar_date))) %>% #day of the year
#   group_by(doy) %>% 
#   summarise(evi_mean = mean(value * as.numeric(scale))) %>% 
#   distinct(doy, .keep_all = TRUE)
#   head(mean_evi)
#   
# ggplot(mean_evi, aes(x = doy, y = evi_mean)) +
#   geom_point() +
#   geom_smooth(method = 'loess') +
#   labs(x = 'Day of Year (Warm Season Only)', y = 'EVI (mean)') +
#   theme_classic() 

#rasterize and reproject
evi_split = evi %>% split(evi$calendar_date)
evi_rast_list = lapply(evi_split, function(x) {mt_to_terra(x,reproject = TRUE)})

#convert to raster stack
evi_rasts <- rast(evi_rast_list)
st_crs(evi_rasts)
poi = st_transform(poi,st_crs(evi_rasts))
evi_rasts <- mask(evi_rasts, poi)

evi = mean(evi_rasts,na.rm=TRUE)
png(here(paste0("data/clean/viz/mean_evi_october_may_2024.png")),width = 6, height = 4, units = "in", res = 300)
plot(evi,main="Mean EVI: May-October 2024")
dev.off()

writeRaster(evi, here(paste0("data/clean/poi/san_diego/mean_evi_october_may_2024.tif")), overwrite = TRUE)


# plot(evi_rasts$`2014-05-01`, main = '2014-05-01')
# plot(evi_rasts$`2024-05-24`, main = '2024-05-24')


# evi_split <- evi %>% 
#   mutate(doy = yday(as_date(calendar_date))) %>% #day of the year
#   group_by(doy) %>% 
#   mutate(evi_mean = mean(value * as.numeric(scale))) %>% 
#   distinct(doy, .keep_all = TRUE)
#?get_nlcd_annual()


