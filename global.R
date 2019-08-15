library(SpaDES)

setPaths(modulePath = file.path("/home/ieddy/Git"),
         inputPath = file.path("../../temp"),
         outputPath = file.path("../../temp"),
         cachePath = file.path("../../temp"))

getPaths() # shows where the 4 relevant paths are

times <- list(start = 0, end = 10)

parameters <- list(
  #.progress = list(type = "text", interval = 1), # for a progress bar
  ## If there are further modules, each can have its own set of parameters:
  #module1 = list(param1 = value1, param2 = value2),
  #module2 = list(param1 = value1, param2 = value2)
)
modules <- list("simpleHarvest")
objects <- list()
inputs <- list()
outputs <- list()

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects)

mySimOut <- spades(mySim)