library(SpaDES)
library(raster)
library(sf)
library(magrittr)

setPaths(modulePath = file.path("../"),
         inputPath = file.path("inputs"),
         outputPath = file.path("ouputs"),
         cachePath = file.path("cache"))

paths <- getPaths() # shows where the 4 relevant paths are
times <- list(start = 0, end = 2)

parameters <- list(
  #.progress = list(type = "text", interval = 1), # for a progress bar
  ## If there are further modules, each can have its own set of parameters:
  #module1 = list(param1 = value1, param2 = value2),
  simpleHarvest = list(.plotInitialTime = 1,
                       .plotInterval = 1)
)
studyArea <- shapefile("C:/Ian/Campbell/RIA/Land-R/inputs/ftStJohn_studyArea.shp")
harvestAreas <- studyArea %>%
  rgeos::gUnaryUnion(spgeom = ., id = 'TSNMBRDSCR')
harvestAreas$AnnualBiomassHarvestTarget <- 2115000 * 0.45  #cubic metres x average Mg/m3 of wood https://cedarstripkayak.wordpress.com/lumber-selection/162-2/ 
#Probably need a more refined approach for estimating volume - This is lumber when LandR has full trees
#For future reference
#Mackenzie is 4500000, Dawson is 1860000, Ft Nelson is 2582350, Prince George is 8350000
#Note that Ft Nelson needs to be labeled as Ft Nelson in the actual Final RIA TSA file (it was dropped when merging)

dem <- raster("C:/Ian/Data/Elevation/GMTED2010 West Canada 500/GMTED2010N50W150_150/50n150w_20101117_gmted_med150.tif")
parks <- shapefile("C:/Ian/Data/Protected Areas/BC/TA_PARK_ECORES_PA_SVW/TA_PEP_SVW_polygon.shp") %>%
  spTransform(., crs(studyArea)) %>%
  crop(., studyArea) 

#Try combining this with buffered rivers
rivers <- Cache(shapefile, "C:/Ian/Data/Hydrography/USA_CAN_rivers.shp/") %>%
  spTransform(., CRSobj = crs(parks)) 
rivers <- Cache(rgeos::gIntersection, rivers, studyArea) %>%
  sf::st_as_sf(.) %>%
  sf::st_buffer(., 250)

parks <- sf::st_as_sf(parks)
protAreas <- Cache(sf::st_union, parks, rivers) %>% 
  sf::as_Spatial(.)

modules <- list("simpleHarvest")
objects <- list("DEMraster" = dem,
                'areasToExclude' = protAreas)
inputs <- list()
outputs <- list()

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects)

mySimOut <- spades(mySim)
