##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#### 03: Build Location-Allocation ####

#coded by: Aja Sutton
#Climate Change Epidemiology Lab
# Scripps Institution of Oceanography, UC San Diego

#load packages; if on CRAN and not installed, install and load.
{pacman::p_load(tidyverse,sf,janitor,here,tidycensus,terra,Require,viridis,
                tigris,crsuggest,data.table,exactextractr,locationallocation)
  
  # Automatically check, download locationallocation package from GitHub if missing, and load it
  Require::Require("giacfalk/locationallocation")
  
  options(scipen=999)
  options(tigris_use_cache = TRUE)
  
  tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)
  
  here::i_am("code/03_build_location_allocation.R")
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

{
#define place of interest
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


#spatial risk weight surface
wt =rast(here("data/clean/poi/san_diego/san_diego_weighted_risk_surface.tif"))
wt = project(wt,"EPSG:4269")
plot(wt)

#heat surface
heat = rast(here("data/clean/heat_index_n_days_103fthreshold_average.tif"))
heat = project(heat,"EPSG:4269")
plot(heat)

#population density (general)
dens = rast(here('data/clean/poi/san_diego/san_diego_general_population_density.tif'))
dens = project(dens,"EPSG:4269")
plot(dens)

# #cumulative heat and risk
# heatalign = resample(heat,wt, method = "bilinear")
# stack = c(wt,heatalign)
# stack = sum(stack,na.rm=TRUE)
# plot(stack)
# hist(stack)

#popdensity and risk
densalign = resample(dens,wt, method = "bilinear")
stack = c(wt,densalign)
stack = sum(stack,na.rm=TRUE)
plot(dens)
plot(wt)
plot(stack-dens)
hist(stack)

rm(heatalign,densalign)
gc()

#cooldens#cool zones (2023 SD ones for now)
#more: https://www.news.caloes.ca.gov/cooling-centers-and-resources-2/
cz = st_read(here('data/HLTH_COOL_ZONES_PUBLIC.geojson'))
# st_crs(cz)
# st_crs(poi)
cz = st_transform(cz,st_crs(poi))
cz = st_intersection(cz,poi)
# st_crs(cz)
# st_crs(poi)

#the traveltime() call can take a long time! be prepared!
#also, it can time out -- if you get a "bad gateway" 
#message it usually means little and you should try again.

#some basic steps it's undertaking:
#calling in the friction and osm data (big download, not that slow)
#cleaning everything up, and calling gdistance::accCost() which 
#calculate the accumulated cost of (below, in this case) walking 
#using road network data from OSM.
#In our case below, it's a linear model (lm) -- this is more than adequate for
#the purposes of selecting cooling centers, and is faster and easier than 
#the other option, random forest (rf).
#When cache = TRUE, it caches the OSM and Malaria Atlas friction file for later.


#### Set Method ####

method = "fastest" # "walk" or "fastest"

#THIS IS SLOW. It's ok! Take a break! Drink some water.
#Using the San Diego city example, it takes about 10 minutes.

if (method == "walk"){
  starttime = Sys.time()
  print(starttime)

  traveltime =
    cz |>
    traveltime(
      bb_area = poi,
      mode = "walk",
      res_output = 100,
      cache = TRUE
    )

  endtime = Sys.time()
  endtime - starttime
  save(traveltime,file=here('data/clean/poi/san_diego/sd_traveltime_localloc_walking.RData'))
}

if (method == "fastest"){
  starttime = Sys.time()
  print(starttime)

  traveltime =
    cz |>
    traveltime(
      bb_area = poi,
      mode = "fastest",
      res_output = 100,
      cache = TRUE
    )

  endtime = Sys.time()
  endtime - starttime

  save(traveltime,file=here('data/clean/poi/san_diego/sd_traveltime_localloc_fastestmethod.RData'))
}



#walking or fastest method available?
#load(here('data/clean/poi/san_diego/sd_traveltime_localloc_walking.RData'))
#load(here('data/clean/poi/san_diego/sd_traveltime_localloc_fastestmethod.RData'))



#visualize and get basic stats
traveltime |>
  traveltime_plot(
    bb_area = poi,
    facilities = cz,
    contour_traveltime = 10) 
ggsave(file=here(paste0('data/clean/poi/san_diego/traveltime_viz',method,'.png')),height=8,width=8, unit="in",dpi=300)


#slow! not working well...
traveltime |>
    traveltime_stats(
      demand = heat,
      breaks = c(5, 10, 15, 30), #minutes
      objectiveminutes = 10)
ggsave(file=here(paste0('data/clean/poi/san_diego/traveltime_stats_',method,'.png')),height=8,width=8, unit="in",dpi=300)

crs = terra::crs(dens)
dens = raster::raster(dens)
raster::crs(dens) = crs

wt= raster::raster(wt)
raster::crs(wt) = crs

raster::crs(traveltime$travel_time) = crs

stack = raster::raster(stack)
raster::crs(stack) = crs

raster::compareRaster(traveltime$travel_time,dens,wt,stack, 
                      extent = TRUE, rowcol = FALSE, crs = TRUE, 
                      res = FALSE, orig = FALSE, stopiffalse = FALSE)

poi = st_crop(poi,raster::extent(traveltime$travel_time))
cz = st_crop(cz,raster::extent(traveltime$travel_time))

dens = resample(dens, traveltime$travel_time, method = "bilinear")
wt = resample(wt, traveltime$travel_time, method = "bilinear")
stack = resample(stack, traveltime$travel_time, method = "bilinear")

#they all match now
raster::compareRaster(traveltime$travel_time,dens,wt,stack, 
                      extent = TRUE, rowcol = FALSE, crs = TRUE, 
              res = FALSE, orig = FALSE, stopiffalse = FALSE)

cz = st_transform(cz,st_crs(poi))

plot(heat)

allocation_data =
  stack |>
    allocation(
      bb_area = poi,
      facilities = cz,
      traveltime = traveltime,
      weights = heat,
      objectiveminutes = 10,
      objectiveshare = 0.99,
      heur = "max",
      approach = "norm",
      exp_demand = 1,
      exp_weights = 1
    )

allocation_data |> 
  allocation_plot(poi) 
ggsave(file=here(paste0('data/clean/poi/san_diego/san_diego_allocation_continuous_map_',method,'.png')),height=8,width=8, unit="in",dpi=300)

temp + geom_sf(data=cz,colour = "blue")

plot(dens)
ggsave(file=here('data/clean/poi/san_diego/san_diego_allocation_continuous_map.png'),height=8,width=8, unit="in",dpi=300)



#discrete
locationallocation::
allocation_data <-
  dens |>
  allocation_discrete(
    bb_area = poi,
    candidate = poi |> st_sample(20),
    facilities = cz,
    n_fac = 5,
    n_samples = 100,
    traveltime = traveltime_data,
    weights = NULL,
    objectiveminutes = 15,
    objectiveshare = NULL,
    approach = "norm",
    exp_demand = 1,
    exp_weights = 1,
    par = FALSE
  )


allocation_data <-
  dens |>
  allocation(
    bb_area = poi,
    facilities = cz,
    traveltime = traveltime,
    weights = NULL,
    objectiveminutes = 15,
    objectiveshare = 0.99,
    heur = "max",
    approach = "norm",
    exp_demand = 1,
    exp_weights = 1
  )
ggsave(file=here('data/clean/poi/san_diego/san_diego_allocation_continuous.png'),height=8,width=8, unit="in",dpi=300)

