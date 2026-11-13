##### Cooling Centers in California ####
##### A Location-Allocation Problem ####

#### 02: Build Spatial Risk Surface (i.e. Weights for Location-Allocation) ####

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
  
  here::i_am("code/02_build_spatial_weights.R")
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


#load data


#### Raster Data ####

#Tree Canopy:
trees = rast(here("data/clean/poi/san_diego/2025_sandiego_nlcd_tree_canopy.tif"))
trees[trees>100] = NA
breaks <- global(trees, fun = quantile, probs = seq(0, 1,0.25), na.rm = TRUE)
breaks <- as.numeric(breaks)

#rev: simplify and downweight tree canopy -- we want less to be weighted higher for NEED
newvals = rev(seq(1,n_distinct(breaks),by=1))

retrees = cbind(lag(unique(breaks),n=1),unique(breaks),newvals)
retrees[1,1] = 0
retrees = as.matrix(retrees,byrow=TRUE)
trees = classify(trees,retrees,include.lowest = TRUE)

plot(trees)

#log-Impervious Surface:
imper = rast(here("data/clean/poi/san_diego/2025_sandiego_nlcd_fct_impervious.tif"))
breaks <- global(imper, fun = quantile, probs = seq(0, 1,0.25), na.rm = TRUE)
breaks <- as.numeric(breaks)
breaks = unique(breaks)
#simplify impervious surface -- we want less to be weighted higher for NEED
newvals = seq(1,n_distinct(breaks),by=1)
reimper = cbind(lag(unique(breaks),n=1),unique(breaks),newvals)
reimper = ifelse(is.na(reimper),0,reimper)
reimper = as.matrix(reimper,byrow=TRUE)
imper = classify(imper,reimper,include.lowest = TRUE)
plot(imper)


#EVI (enhanced vegetation index)
evi = rast("data/clean/poi/san_diego/mean_evi_october_may_2024.tif")
evi[evi < -1 | evi > 1] <- NA
evi = log(evi)
breaks <- global(evi, fun = quantile, probs = seq(0, 1,0.25), na.rm = TRUE)
breaks <- as.numeric(breaks)
breaks = unique(breaks)
breaks
hist(evi)
#simplify and downweight evi -- we want less to be weighted higher for NEED
newvals = rev(seq(1,n_distinct(breaks),by=1))
reevi = cbind(unique(breaks),newvals)
lag = lag(unique(breaks),n=1)
reevi = cbind(lag,reevi)
reevi = ifelse(is.na(reevi),reevi[1,2],reevi)
reevi = as.matrix(reevi,byrow=TRUE)

#could blanket convert all negative values to 0 and binarize?
evi = classify(evi,reevi,include.lowest = TRUE)
plot(evi)


#create weight surface with these values
st_crs(evi)
st_crs(trees)
st_crs(imper)

evi = project(evi, imper, method = "bilinear")


#### Census Tract- and Block Group-level Data (polygons) ####

#CDC SVI: breaks
svi = st_read(here('data/clean/poi/san_diego/2022_cdc_svi_san_diego.gdb'),
              layer = "_2022_cdc_svi_san_diego_breaks")
svi = st_transform(svi,crs_choice)

svi = svi %>%
  clean_names() %>%
  select(namelsad,rplthemes,fips,location,SHAPE) %>%
  rename(geometry = SHAPE)

#AC Coverage:
ac = readRDS(here("data/clean/ac_coverage_sandiego.rds"))
ac = st_transform(ac,crs_choice)

#prioritize NO AC access (percent)
ac = ac %>%
  select(place,name,geoid,no_ac_pe,hse_occ_e,geometry) %>%
  filter(no_ac_pe > 0) %>%
  mutate(no_ac_pe = no_ac_pe/100)

ac$no_ac_br = 
  cut(ac$no_ac_pe, 
      breaks = seq(0,max(ac$no_ac_pe),by=0.1), 
      include.lowest = TRUE,
      labels = c(1,2,3,4,5,6))

ac %>%
  ggplot() +
  geom_sf(aes(fill=no_ac_br)) +
  scale_fill_viridis(discrete=TRUE,drop = FALSE)

ac %>%
  ggplot() +
  geom_sf(aes(fill=no_ac_pe)) +
  scale_fill_viridis()

#vulnerable-to-heat population density: ie >65 years and <18 years

dens =readRDS(here('data/clean/ca_population_density_bysexage_5yr_acs_2023.rds'))

#filter to place of interest
dens = st_intersection(dens,poi)

geoidtotal = dens %>%
  select(geoid,variable,estimate) %>%
  filter(variable =="total") %>%
  rename(totals =estimate) %>%
  distinct() %>%
  select(-variable) %>%
  st_drop_geometry()

dens =
  dens %>%
  select(geoid,name,variable,estimate,popdensity,aland,namelsad,geometry) %>%
  group_by(geoid) %>%
  filter(str_detect(variable,"65|67|70|75|80|85|under 5|10|15")) %>%
  mutate(vuln_pop = sum(estimate),
         vulnerable_popdens = vuln_pop/(aland/1000000),
         ) %>%
  select(geoid,name,namelsad,vuln_pop,vulnerable_popdens,geometry) %>%
  distinct() 

dens = left_join(dens,geoidtotal, "geoid")

dens = 
  dens %>%
  group_by(geoid) %>%
  mutate(vulnerable_popprop = vuln_pop/totals) %>%
  ungroup() %>%
  select(geoid,name,namelsad,
         vuln_pop,vulnerable_popdens,vulnerable_popprop,
         geometry)

dens %>%
  ggplot() +
  geom_sf(aes(fill=vulnerable_popdens))
dens %>%
  ggplot() +
  geom_sf(aes(fill=vulnerable_popprop))

#quantiles to put them on a similar scales to others
quantiles = quantile(dens$vulnerable_popdens, c(0,0.25, 0.5, 0.75), q = c(0,0.25, 0.5, 0.75))

dens = 
  dens %>%
  mutate(vulnpop_q =case_when(vulnerable_popdens < quantiles[2] ~ 1,
                      vulnerable_popdens < quantiles[3] ~ 2,
                      vulnerable_popdens < quantiles[4] ~ 3,
                     vulnerable_popdens >= quantiles[4] ~ 4,
                      is.na(vulnerable_popdens) ~ 1,
                     .default=999)) 

dens %>% filter(vulnpop_q == 999)

#join polygons
# plot(ac$geometry,col=NA, border="black",type="l")
# plot(svi$geometry,col=NA, border="blue",type="l",add=TRUE)

svi = 
  svi %>% 
  st_drop_geometry()

ac = ac %>%
  mutate(geoid = str_remove_all(geoid,"1400000US"),
         name = str_replace_all(name,",",";")) %>%
  select(-no_ac_pe)

wt_polygon = left_join(svi,ac, by = c("fips"="geoid","namelsad"="place","location"="name"))
# glimpse(wt_polygon)
# glimpse(dens)
# wt_polygon = left_join(wt_polygon,dens,by=c("fips"="geoid","namelsad"="namelsad","location"="name"))


wt_polygon = 
  wt_polygon %>%
  st_as_sf() %>%
  filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON"))

wt_polygon = 
  wt_polygon %>%
  select(namelsad,rplthemes,no_ac_br,geometry) %>% #vulnerable_popprop,
  rename(place = namelsad,
         svi = rplthemes,
         no_ac = no_ac_br) %>%
  #mutate(svi = factor(svi,levels=seq(1,4,by=1))) %>% 
  mutate(no_ac = as.double(no_ac),
         svi = as.double(svi)) %>%
  st_as_sf() %>%
  filter(!st_is_empty(.)) 



#impute missing SVI and AC values from nearest neighbours

#only a few missing
sum(is.na(wt_polygon$no_ac))
sum(is.na(wt_polygon$svi))

index = st_touches(wt_polygon,wt_polygon)
wt_polygon =
  wt_polygon %>% 
  mutate(id = 1:nrow(.),
         no_ac_impute = ifelse(is.na(no_ac),
                               apply(index, 1, function(i){mean(.$no_ac[i])}),
                               no_ac)) %>%
  group_by(id) %>%
  mutate(weight = rowMeans(across(c(no_ac_impute,svi)),na.rm=TRUE)) %>%
  ungroup() 

wt_polygon %>%
  ggplot() +
  geom_sf(aes(fill=weight)) +
  scale_fill_viridis()

glimpse(wt_polygon)
#rasterize wt_polygon and join to stack
wt = rasterize(wt_polygon, trees, field = "weight")
dens = dens %>%  st_collection_extract("POLYGON")
dens = st_transform(dens,st_crs(wt))
wt_dens = rasterize(dens, trees, field = "vulnpop_q")

plot(trees)
plot(evi)



#mean of:
#tree canopy, evi, impervious surfaces, 
#cdc svi, no AC, vulnerable population density (<18 years and 65 and older)
stack = c(trees,evi,imper,wt,wt_dens)
stack = terra::mean(stack,na.rm=TRUE)
plot(stack)
hist(stack)  #normal distribution
freq(stack) #5 is rare, 4 is also uncommon -- good allocation weighting

#what if we tae 
# stack2 = c(trees,evi,imper,wt)#,wt_dens)
# stack2 = terra::mean(stack2,na.rm=TRUE)
# plot(stack2)
# 
# diffstack = stack2-stack
# plot(diffstack)
# 
# hist(stack2)
# hist(stack)

writeRaster(stack, here("data/clean/poi/san_diego/san_diego_weighted_risk_surface.tif"),overwrite = TRUE)


### Save general population density as poi-specific

dens =readRDS(here('data/clean/ca_population_density_bysexage_5yr_acs_2023.rds'))

dens = st_transform(dens,st_crs(poi))

dens = 
  dens %>%
  filter(variable=="total") %>%
  select(geoid,popdensity,geometry)

#filter to place of interest
dens = st_intersection(dens,poi)
dens = dens[st_geometry_type(dens) %in% c("POLYGON", "MULTIPOLYGON"),]

dens$popdensity = as.double(dens$popdensity)
quantiles = quantile(dens$popdensity, seq(0,1,by=0.1), q = c(1:10))

dens = 
  dens %>%
  mutate(popdens_q =case_when(popdensity == 0 ~ 0,
                              popdensity < quantiles[2] ~ 1,
                              popdensity < quantiles[3] ~ 2,
                              popdensity < quantiles[4] ~ 3,
                              popdensity < quantiles[5] ~ 4,
                              popdensity < quantiles[6] ~ 5,
                              popdensity < quantiles[7] ~ 6,
                              popdensity < quantiles[8] ~ 7,
                              popdensity < quantiles[9] ~ 8,
                              popdensity < quantiles[10] ~ 9,
                              popdensity >= quantiles[10] ~ 10,
                              is.na(popdensity) ~ 1,
                              .default=999)) 

dens %>%
  ggplot() +
  geom_sf(aes(fill=popdens_q,color = after_scale(fill)),linewidth=0.1) +
  scale_fill_viridis()
ggsave(file=here('data/clean/poi/san_diego/populationdensity_quantiles.png'))

dens = st_transform(dens,st_crs(trees))
densrast = rasterize(dens, trees, field = "popdens_q")

writeRaster(densrast, here("data/clean/poi/san_diego/san_diego_general_population_density.tif"),overwrite = TRUE)



