
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "simpleHarvest",
  description = paste("This is a very simplistic harvest module designed to interface with the LandR suite of modules",
                      "It will create a raster of harvested patches, but will not simulate actual harvest.",
                      "Should be paired with LandR_reforestation",
  keywords = c("harvest", "LandR", "rstCurrentHarvest"),
  authors = c(person(c("Ian"), "Eddy", email = "ian.eddy@canada.ca", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.5.9008", simpleHarvest = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "simpleHarvest.Rmd"),
  reqdPkgs = list("PredictiveEcology/LandR@development", "sp", "raster", 'sf', 'magrittr', 'fasterize'),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated?
                    This is generally intended for data-type modules, where stochasticity and time are not relevant"),
    defineParameter("harvestTarget", numeric, 0.01, 0, 1,
                    desc= "proportion of harvestable area to harvest each timestep"),
    defineParameter("minAgesToHarvest", "numeric", 50, 1, NA, desc =  "minimum ages of trees to harvest"),
    defineParameter("maxPatchSizetoHarvest", numeric, 10, 1, NA,
                    desc = "maximum size for harvestable patches, in pixels")
  ),
  inputObjects = bind_rows(
    expectsInput('cohortData', objectClass = 'data.table', desc = "table with pixelGroup, age, species, and biomass of cohorts"),
    expectsInput(objectName = 'demRaster', objectClass = 'RasterLayer', desc = 'a DEM used to exclude areas from harvest. If not provided,
                 there will be no elevation restriction applied to harvest', sourceURL = NA),
    expectsInput("pixelGroupMap", objectClass = "RasterLayer", desc = "Raster giving locations of pixelGroups"),
    expectsInput('rasterToMatch', objectClass = 'RasterLayer',
                 desc = 'a template raster for all raster operations. Cannot have lat/long projection'),
    expectsInput("studyArea", objectClass = "SpatialPolygonsDataFrame", desc = "study area used to crop spatial inputs, if applicable",
                 sourceURL = "https://drive.google.com/open?id=1TlBfGfes_6UQW4M3jib8zgY5sGd_yjmY"),
    expectsInput("thlb", objectClass = "RasterLayer",
                 desc = "binary raster with 1 denoting harvestable pixels. Must match spatial attributes of rasterToMatch")
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = 'rstCurrentHarvest', objectClass = 'RasterLayer', desc = 'Binary raster representing annual harvested areas'),
    createsOutput(objectName = 'harvestExclusionRaster', objectClass = "RasterLayer", desc = "Binary raster representing areas that will never be harvested"),
    createsOutput(objectName = "harvestIndex", objectClass = "data.table", desc = "data.table with cell indices for harvest areas")
  )
))

## event types
#   - type `init` is required for initialization

doEvent.simpleHarvest = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim), "simpleHarvest", "harvest")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "simpleHarvest", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "simpleHarvest", "save")

    },
    plot = {
      plot(sim$rstCurrentHarvest)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "simpleHarvest", "plot")

    },
    save = {

      sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "simpleHarvest", "save")

    },
    harvest = {

      if (LandR::scheduleDisturbance(sim$rstCurrentHarvest, time(sim))) {
        sim$rstCurrentHarvest <- harvestSpreadInputs(pixelGroupMap = sim$pixelGroupMap,
                                                     cohortData = sim$cohortData,
                                                     thlb = sim$thlb,
                                                     maxCutSize = P(sim)$maxPatchSizetoHarvest,
                                                     target = P(sim)$harvsetTarget,
                                                     minAgesToHarvest =  = P(sim)$minAgesToHarvest,
                                                     harvestIndex = sim$harvestIndex)
      }

      sim <- scheduleEvent(sim, time(sim) + 1,  "simpleHarvest", "harvest")

    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {

  return(invisible(sim))
}


### template for save events
Save <- function(sim) {

  sim <- saveFiles(sim)

  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {

  plot(sim$rstCurrentHarvest)

  return(invisible(sim))
}

harvestSpreadInputs <- function(pixelGroupMap,
                                cohortData,
                                exclusionAreas,
                                harvestAreas,
                                maxCutSize,
                                minAgesToHarvest,
                                target,
                                harvestIndex) {

  browser()
  #This function calculates which areas are available to cut, how many cuts of mean size would be needed to achieve the target
  # based on the mean biomass in a pixel.
  #Make an ageMap
  cohortData <- copy(cohortData)
  standAges <- cohortData[, .(BweightedAge = sum(B * age)/sum(B)), .(pixelGroup)]

  pixID <- data.table('pixelGroup' = getValues(pixelGroupMap),
                      "pixelIndex" = 1:ncell(pixelGroupMap),
                      "thlb" = getValues(thlb))
  pixID <- na.omit(pixID) #drop non-forested pixels, even if they're in thlb
  pixID <- pixID[thlb == 1,] #drop pixels that aren't in thlb

  landStats <- standAges[pixID, on = c("pixelGroup")]
  landStats <- landStats[BweightedAge >= minAgesToHarvest,]
  #I assume this method is faster than matching the pixelGroup raster values with a vector
  harvestableAreas <- raster(thlb)
  thlb[landStats$pixelID] <- 1 #this will be spread probability


  #calculate harvest target
  harvestTarget <- round(nrow(landStats) * target)
  minCuts <- round(harvestTarget/maxCutSize)

  #if every cut reaches maximum size, you need minCuts to begin with
  initialCuts <- sample(landStats$pixelIndex, size = minCuts, replace = FALSE)

  firstIteration <- spread2(landscape = harvestableAreas,
                            start = minCuts,
                            spreadProb = 0.3,
                            maxSize = maxCutSize)

  return(harvest)
}

.inputObjects <- function(sim) {

  if (!suppliedElsewhere(object = c("rasterToMatch", "studyArea", "cohortData", "thlb", "pixelGroupMap"), sim)) {
    stop("please supply all objects - this module has no defaults at this time")
  }
}
