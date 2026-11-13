##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#coded by: Aja Sutton
#Climate Change Epidemiology Lab
# Scripps Institution of Oceanography, UC San Diego

#load packages; if on CRAN and not installed, install and load.
{pacman::p_load(tidyverse,janitor,here,terra,Require,
                prism,tigris,exactextractr,heatindex)
  
  
  options(scipen=999)
  options(tigris_use_cache = TRUE)

  here::i_am("code/00a_get_prism_heat_index.R")
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

#### Goal: establish heat index surface to use with weighted risk surface ####

#### Get PRISM Data: monthly max. air temperature and mean dewpoint for heat index ####

#ultimately we want to calculate the maximum heat index for our places.
#To do this, we need maximum air temperature and relative humidity.
#PRISM doesn't provide relative humidity, but this can be calculated
#from dewpoint with the air temperature data:
#First: https://qed.epa.gov/hms/meteorology/humidity/algorithms/
# https://bmcnoldy.rsmas.miami.edu/Humidity.html 
#Second: https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml



#we want monthly maximum temperature at a fine geographic granularity
prism_set_dl_dir(here('data/prism_2014_2024'))

# ##ONLY NEED TO RUN ONCE (daily values are slow to download!!!):
# ##Download national data, crop to CA only to save space (need to manually delete national files)
# get_all_dates = function(tvalue,dates,i,cashape) {
# 
#   Require("prism","terra")
# 
#   get_prism_dailys(
#     type = tvalue,
#     resolution = "800m",
#     dates = dates,
#     keepZip = FALSE)
#   
#   myraster = prism_archive_subset(tvalue,"daily", mon = 5:10, year = i, resolution = "800m")
#   myraster =  pd_to_file(myraster)
#   myraster = rast(myraster)
#   myraster = crop(myraster,cashape)
#   myraster = mask(myraster,cashape)
#   
#   writeRaster(myraster, file=here(paste0("data/ca_only_prism_2019_2024_daily_by_year/",i,"_",tvalue,"_ca_warm_season_dailys.tif")), overwrite = TRUE)
#   
# }
# 
# 
# 
# #To do: add to lab's google cloud storage and call from there instead so not saved locally.
# 
# #look at previous 5 available years' worth of' data
# years = 2019:2024
# cashape = tigris::states() #get shapefile of the state of CA for crop+mask
# cashape = cashape %>% filter(NAME=="California")
# 
# #May-October only for now
# #to avoid downloading national data, we specify CA data crop+mask in get_all_dates()
# for (i in years) {
#   
#   starttime = Sys.time()
#   print(starttime)
#   
#   dates = seq(as.Date(paste0(i,"-05-01")), as.Date(paste0(i,"-10-31")), by = "days")
#   print(dates)
#   
#   get_all_dates(tvalue="tmax",dates,i,cashape)
#   get_all_dates(tvalue="tdmean",dates,i,cashape)
#   
#   endtime = Sys.time()
#   print(paste0("Operation time for",i,": ",endtime-starttime))
#   
# }

#original code where we didn't crop in the get_all_dates function call:
# #Instead of keeping these data at the nationa level (too big), we crop down to
# #the state of California only for now:
# clean_ca_raster = function(tvalue,cashape,monthid,yearid){
#   
#   require("terra","prism")
#   
#   myraster = prism_archive_subset(tvalue,"daily", mon = 5:10, year = i, resolution = "800m")
#   myraster =  pd_to_file(myraster)
#   myraster = rast(myraster)
#   myraster = crop(myraster,cashape)
#   myraster = mask(myraster,cashape)
#   
#   writeRaster(myraster, file=here(paste0("data/ca_only_prism_2019_2024_daily_by_year/",i,"_",tvalue,"_ca_warm_season_dailys.tif")), overwrite = TRUE)
# }
# 
# #and crop
# 
# years = 2022:2023
# for (i in years) {
#   clean_ca_raster(tvalue="tmax",cashape=cashape,monthid = 5:10, yearid = i)
#   clean_ca_raster(tvalue="tdmean",cashape=cashape,monthid = 5:10, yearid = i)
# }


#Again, above should only need to be run once when the national data are downloaded,
#and then cleaned to the correct state level to keep files small (delete national files
#manually, if desired)

#clean raster to reflect ONLY place of interest (our example is San Diego)
clean_poi_raster = function(filename,tvalue,placeshape,yearid){
  
  require("terra")
  files <- list.files(here("data/ca_only_prism_2019_2024_daily_by_year"),
                      pattern = "\\.tif$", 
                      full.names = FALSE)
  files = grep(as.character(yearid),files,value = TRUE)
  files = grep(as.character(tvalue),files,value = TRUE)
  
  myraster =rast(here(paste0("data/ca_only_prism_2019_2024_daily_by_year/",files)))
  myraster = crop(myraster,placeshape)
  myraster = mask(myraster,placeshape)
  
  return(myraster)
  
}



#Filter down to POI: place of interest
place_name = "San Diego city"
poi = tigris::places(state = "CA", cb = TRUE)
poi = poi %>% 
  clean_names() %>%
  filter(str_detect(namelsad,place_name))


#So, to get a reasonable idea of extreme heat exposure in our place of interest
#we calculate the Heat Index with the above data as follows:

#For each year of data, load in some "tvalue" data of interest (tmax, tdmean),
#crop and mask to our area of interest (placeshape)
#We add monthid as a convenience to load all months' daily data; 
#can load individual months if for some reason desired.



#Cooling centers are places people go to during heatwaves. Heatwaves are 
#anomalous temperature events we can count in days to identify pixels with lots of
#very hot days. So, we need to pull daily max temperature and dewpoint data to 
#figure out which places have lots of heatwaves (e.g., max temp greater than 90%)
#to get the relative humidity, to finally calculate the HEAT INDEX for each pixel.
#Then, we can count the number of days each pixel's heat index is above some 
#threshold. For now, this threshold is the National Weather Service's heat index 
#Danger category (103F) to make it easy to use by planning folks.

#first we test with one month's worth of data

get_heat_index_raster = function(yearid,placeshape){
  
  require("terra","prism")
  
  
  #loads all years' data and masks to poi shape
  tmaxrast = clean_poi_raster(tvalue="tmax",placeshape=poi,yearid=yearid)#,monthid=5:10,yearid=yearid)
  tdmeanrast = clean_poi_raster(tvalue="tdmean",placeshape=poi,yearid=yearid)#,monthid=5:10,yearid=yearid)
  
  
  #calculate relative humidity using tmax and tdmean
  rhrast = 100*(exp((17.625*tdmeanrast)/(243.04+tdmeanrast))/exp((17.625*tmaxrast)/(243.04+tmaxrast)))
  
  #get heat index
  hi_c = function(temptmax, temprh) {
    hi = heatindex::heatindex(
      T = temptmax + 273.15, #convert to Kelvin
      rh = temprh/100 #scale of 0 to 1
    )
    #convert back to Celsius
    hi = hi - 273.35
  }
  
  #lapp hi_c across all daily layers; automatically handles parallel 
  #vector operations across spatial grids
  hiraster =  lapp(
    x = sds(tmaxrast, rhrast), #Create a SpatRasterDataset container
    fun = hi_c
  )
  
  #fix layer names to match original calendar dates
  names(hiraster) <- names(tmaxrast)
  
  #greatest number of days with extreme heat during warm months in the last 5 years.
  #figure out: for that day of the year, across years of interest,
  #how many EH days does each pixel have?
  
  #We can use https://www.weather.gov/ama/heatindex to identify a reproducible 
  #and defensible threshold with our heat index data. For now, let's say that 
  #we want to set the threshold at the National Weather Service's Danger band: 103F.
  
  #convert 103F to Celsius (data in Celsius)
  
  threshold = (103-32) * 5/9
  
  #not used because poor logic: 0.95 threshold for whole area of focus minmax(threshold)[2]
  ehraster = ifel(hiraster >= threshold,1,0)
  temp = sum(ehraster,na.rm=TRUE)
  
  png(here(paste0("data/ca_only_prism_2019_2024_daily_by_year/threshold_days/viz/",i,"heat_index_threshold_days.png")),width = 6, height = 4, units = "in", res = 300)
  plot(temp,main="N days >=103F \n NWS threshold")
  dev.off()
  
  writeRaster(temp, here(paste0("data/ca_only_prism_2019_2024_daily_by_year/threshold_days/",i,"_heat_index_n_days_103fthreshold.tif")), overwrite = TRUE)
  
  
}



#build total heat 
years = 2019:2024
for (i in years){
  get_heat_index_raster(yearid=i,placeshape=poi)
}

files = list.files(here("data/ca_only_prism_2019_2024_daily_by_year/threshold_days"),
                    pattern = "\\.tif$", 
                    full.names = TRUE)
rasters = rast(files)

rasters2 = sum(rasters, na.rm = TRUE)
png(here(paste0("data/clean/viz/heat_index_threshold_days.png")),width = 6, height = 4, units = "in", res = 300)
plot(rasters2,main="N days >=103F \n NWS threshold")
dev.off() 
writeRaster(rasters2, here(paste0("data/clean/heat_index_n_days_103fthreshold.tif")), overwrite = TRUE)

#average across years 
rasters = sum(rasters, na.rm = TRUE)/n_distinct(years)

png(here(paste0("data/clean/viz/heat_index_threshold_days_average.png")),width = 6, height = 4, units = "in", res = 300)
plot(rasters,main="N days >=103F \n NWS threshold / N years")
dev.off()

writeRaster(rasters, here("data/clean/heat_index_n_days_103fthreshold_average.tif"), overwrite = TRUE)


