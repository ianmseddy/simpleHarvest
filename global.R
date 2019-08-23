library(SpaDES)
library(raster)

setPaths(modulePath = file.path("../"),
         inputPath = file.path("inputs"),
         outputPath = file.path("ouputs"),
         cachePath = file.path("cache"))

getPaths() # shows where the 4 relevant paths are

times <- list(start = 0, end = 2)

parameters <- list(
  #.progress = list(type = "text", interval = 1), # for a progress bar
  ## If there are further modules, each can have its own set of parameters:
  #module1 = list(param1 = value1, param2 = value2),
  simpleHarvest = list(.plotInitialTime = 1,
                       .plotInterval = 1)
)
dem <- raster("C:/Ian/Data/Elevation/GMTED2010 West Canada 500/GMTED2010N50W150_150/50n150w_20101117_gmted_med150.tif")
parks <- shapefile("C:/Ian/Data/Protected Areas/BC/TA_PARK_ECORES_PA_SVW/TA_PEP_SVW_polygon.shp")
#Try combining this with buffered rivers
modules <- list("simpleHarvest")
objects <- list("DEMraster" = dem,
                'areasToExclude' = parks)
inputs <- list()
outputs <- list()

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects)

mySimOut <- spades(mySim)
