
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "simpleHarvest",
  description = "This is a simple harvest module designed to interface with the LandR suite of modules.", 
  keywords = c("harvest", "LandR", "rstCurrentHarvest"),
  authors = c(person(c("Ian"), "Eddy", email = "ian.eddy@canada.ca", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.5.9008", simpleHarvest = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "simpleHarvest.Rmd"),
  reqdPkgs = list("PredictiveEcology/LandR@development", "sp"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? 
                    This is generally intended for data-type modules, where stochasticity and time are not relevant"),
    defineParameter("ElevationToExclude", "numeric", NA, NA, NA, "Elevation threshold above which areas are excluded from harvest")
    
  ),
  inputObjects = bind_rows(
    expectsInput('areasToExclude', objectClass = 'SpatialPolygonsDataFrame', desc = "A shapefile with all areas to exclude from harvest, e.g. National parks,
                 riparian areas. Defaults to a series of relevant areas in British Columbia"),
    expectsInput('cohortData', objectClass = 'data.table', desc = "table with pixelGroup, age, species, and biomass of cohorts"),
    expectsInput(objectName = 'demRaster', objectClass = 'RasterLayer', desc = 'a DEM used to exclude areas from harvest. If not provided, 
                 there will be no elevation restriction applied to harvest', sourceURL = NA),
    expectsInput("harvestAreas", objectClass = "SpatialPolygonsDataFrame", desc = 'A shapefile with layers representing unique harvest areas that 
    contains (at least) two *IMPORTANT* attributes:
                 - 1) AnnualBiomassHarvestTarget: the annual biomass (Mg) to remove, 
                 - 2) MaximumAllowableCutSize: the maximum allowable size for a single contiguous harvest (ha).'),
    expectsInput("pixelGroupMap", objectClass = "RasterLayer", desc = "Raster giving locations of pixelGroups"),
    expectsInput('rasterToMatch', objectClass = 'RasterLayer', desc = 'a template raster for all raster operations. Cannot have lat/long projection', 
                 sourceURL = NA),
    expectsInput("studyArea", objectClass = "SpatialPolygonsDataFrame", desc = "study area used to crop spatial inputs, if applicable")
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = 'rstCurrentHarvest', objectClass = 'RasterLayer', desc = 'Binary raster representing annual harvested areas'),
    createsOutput(objectName = 'harvestExclusionRaster', objectClass = "RasterLayer", desc = "Binary raster representing areas that will never be harvested")
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
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "simpleHarvest", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "simpleHarvest", "save")
    },
    plot = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      #plotFun(sim) # uncomment this, replace with object to plot
      # schedule future event(s)

      # e.g.,
      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "simpleHarvest", "plot")

      # ! ----- STOP EDITING ----- ! #
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "simpleHarvest", "save")

      # ! ----- STOP EDITING ----- ! #
    },
    event1 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "simpleHarvest", "templateEvent")

      # ! ----- STOP EDITING ----- ! #
    },
    event2 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "simpleHarvest", "templateEvent")

      # ! ----- STOP EDITING ----- ! #
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
  
 sim$harvestExclusionRaster <- generateExclusionRaster()#args to come
  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  #Plot(sim$object)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
Event1 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event2
Event2 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  browser()
  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  
  if (!suppliedElsewhere("studyArea", sim)) {
    message("study area not supplied. Defaulting to Ft. St. John TSA")
    sim$studyArea <- prepInputs(url = extractURL("studyArea", sim),
                                fun = "shapefile", 
                                destinationPath = dataPath(sim),
                                useCache = TRUE, 
                                targetFile = 'ftStJohn_studyArea.shp')
  }
  
  if (!suppliedElsewhere("rasterToMatch", sim)) {
    message("rasterToMatch not supplied. Defaulting to LCC2005. CRS will overide that of studyArea, if one was provided")
    sim$rasterToMatch <- prepInputsLCC(studyArea = sim$studyArea, filename2 = NULL, destinationPath = tempdir(), useCache = TRUE)
    sim$studyArea <- spTransform(sim$studyArea, CRSobj = crs(sim$rasterToMatch))
  }

  if (!suppliedElsewhere("cohortData", sim) & !suppliedElsewhere('pixelGroupMap', sim)) {
    message("either cohortData, pixelGroup, or both are not supplied. Using simulated data. Consider running the module LBMR")
    sim$cohortData <- data.table("pixelGroup" = c(1, 2), "speciesCode" = c("Pice_mar", "Pice_gla"), B = c(1000, 500), age = c(50, 10))
    temp <- sim$rasterToMatch
    temp[] <- sample(c(1,2), size = ncell(temp), replace = TRUE)
  }
  if (!suppliedElsewhere("harvestAreas", sim)) {
    message("harvestAreas not supplied. You should really supply harvestAreas so simulated units correspond to raster. Modifying studyArea instead.")
    sim$harvestAreas <- sim$studyArea
    sim$harvestAreas$AnnualBiomassHarvestTarget <- 1e5
    sim$MaximumAllowableCutSize <- 6.25 * 10 #62.5 hectares. This is for easy calculation with cell size
  }
  
  
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above