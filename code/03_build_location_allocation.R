##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#### 03: Build Location-Allocation ####

#coded by: Aja Sutton
#Climate Change Epidemiology Lab
# Scripps Institution of Oceanography, UC San Diego

#load packages; if on CRAN and not installed, install and load.
{pacman::p_load(tidyverse,sf,janitor,here,tidycensus,terra,Require,viridis,lwgeom,
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
heat = rast(here("data/clean/poi/san_diego/san_diego_heat_index_n_days_103fthreshold_binned.tif"))
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
droplowpopdensity = densalign
droplowpopdensity[droplowpopdensity < 2] = 0
hist(droplowpopdensity)
stack = (densalign * wt) #/ minmax(wt)[2]
plot(stack)
#stack[stack == -Inf] = 0

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
ggsave(file=here(paste0('data/clean/poi/san_diego/traveltime_viz_',method,'.png')),height=8,width=8, unit="in",dpi=300)

#greedy heuristic wants to put points in places with few people,
#so we must force drop the geographies with places without people!
heat = resample(heat,stack,method="bilinear")
droplowpopdensity = densalign
droplowpopdensity[droplowpopdensity < 2] = 0
heatpopdensity =   droplowpopdensity*heat
plot(  heatpopdensity )

# stackheat = mean(c(stack,heat))
# hist(stackheat)



traveltime |>
    traveltime_stats(
      demand = heatpopdensity,
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
plot()
allocation_data =
  stack |>
    allocation(
      bb_area = poi,
      facilities = cz,
      traveltime = traveltime,
      weights = (heatpopdensity),
      objectiveminutes = 10,
      objectiveshare = 0.99,
      heur = "max",
      approach = "norm",
      exp_demand = 1,
      exp_weights = 1
    )

allocation_data |> 
  allocation_plot(poi, facilities = cz) 
ggsave(file=here(paste0('data/clean/poi/san_diego/san_diego_allocation_continuous_map_',method,'.png')),height=8,width=8, unit="in",dpi=300)



#discrete
{
  starttime = Sys.time()
  print(starttime)

  
  poi %>%
  ggplot() +
  geom_sf() +
  geom_sf(data=poi %>% st_sample(1000),col="red")


  discrete_allocation_data <-
    stack |> 
    allocation_discrete(
      bb_area = poi,
      candidate = poi |> st_sample(2000),
      facilities = NULL,
      n_fac = 12,
      traveltime = traveltime,
      weights = heatpopdensity,
      objectiveminutes = 10,
      objectiveshare = NULL,
      approach = "norm",
      #exp_demand = 1,
      #exp_weights = 1,
      par = TRUE
    )
  
  endtime = Sys.time()
  endtime - starttime 
  
  discrete_allocation_data |>
    allocation_plot(poi, facilities = cz)
  
  ggsave(file=here(paste0('data/clean/poi/san_diego/san_diego_allocation_discrete_with_populationheatweight_',method,'.png')),height=8,width=8, unit="in",dpi=300)
  

}



{
#find facilities around existing ones
discrete_allocation_data <-
  stack |> 
  allocation_discrete(
    bb_area = poi,
    candidate = poi |> st_sample(2000),
    facilities = cz,
    n_fac = 2,
    traveltime = traveltime,
    weights = heatpopdensity,
    objectiveminutes = 10,
    objectiveshare = NULL,
    approach = "norm",
    #exp_demand = 1,
    #exp_weights = 1,
    par = TRUE
  )

endtime = Sys.time()
endtime - starttime 

discrete_allocation_data |>
  allocation_plot(poi, facilities = cz)

ggsave(file=here(paste0('data/clean/poi/san_diego/san_diego_allocation_discrete_with_populationheatweight_new_vs_existing_',method,'.png')),height=8,width=8, unit="in",dpi=300)


}




#### Sensitivity of heatpopdensity as (heat * )

