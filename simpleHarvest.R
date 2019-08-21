
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "simpleHarvest",
  description = "This is a simple harvest module designed to interface with the LandR suite of modules. It will create a harvest map, but 
  will not simulate actual harvest. Should be paired with LandR_reforestation", 
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
    defineParameter("ElevationToExclude", "numeric", 1500, NA, NA, "Elevation threshold above which areas are excluded from harvest"),
    defineParameter("minAndMaxAgesToHarvest", "numeric", c(40, 100), NA, NA, desc =  "minimum and maximum ages of trees to harvest")
    #defieParameter("speciesToHarvest", ... need to consier allowing only certain species
    
  ),
  inputObjects = bind_rows(
    expectsInput('areasToExclude', objectClass = 'SpatialPolygonsDataFrame', desc = "A shapefile with all areas to exclude from harvest, e.g. National parks,
                 riparian areas"),
    expectsInput('cohortData', objectClass = 'data.table', desc = "table with pixelGroup, age, species, and biomass of cohorts"),
    expectsInput(objectName = 'demRaster', objectClass = 'RasterLayer', desc = 'a DEM used to exclude areas from harvest. If not provided, 
                 there will be no elevation restriction applied to harvest', sourceURL = NA),
    expectsInput("harvestAreas", objectClass = "SpatialPolygonsDataFrame", desc = 'A shapefile with layers representing unique harvest areas that 
    contains (at least) 3 *IMPORTANT* attributes:
                 - 1) AnnualBiomassHarvestTarget: the annual biomass (Mg) to remove, 
                 - 2) MaximumAllowableCutSize: the maximum allowable size for a single contiguous harvest (ha).
                 - 3) meanCutSize: the mean size of cut blocks (ha)'),
    expectsInput("pixelGroupMap", objectClass = "RasterLayer", desc = "Raster giving locations of pixelGroups"),
    expectsInput('rasterToMatch', objectClass = 'RasterLayer', desc = 'a template raster for all raster operations. Cannot have lat/long projection', 
                 sourceURL = NA),
    expectsInput("studyArea", objectClass = "SpatialPolygonsDataFrame", desc = "study area used to crop spatial inputs, if applicable",
                 sourceURL = "https://drive.google.com/open?id=1TlBfGfes_6UQW4M3jib8zgY5sGd_yjmY")
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
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "simpleHarvest", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "simpleHarvest", "save")
      sim <- scheduleEvent(sim, time(sim), "simpleHarvest", "harvest")
    },
    plot = {
      
      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "simpleHarvest", "plot")

    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "simpleHarvest", "save")

      # ! ----- STOP EDITING ----- ! #
    },
    harvest = {
      # ! ----- EDIT BELOW ----- ! #
      if (LandR::scheduleDisturbance(sim$rstCurrentHarvest, time(sim))) {
        sim$rstCurrentHarvest <- harvestTrees(pixelGroupMap = sim$pixelGroupMap,
                                              cohortData = sim$cohortData,
                                              exclusionAreas = sim$harvestExclusionRaster,
                                              harvestAreas = sim$harvestAreas,
                                              ageWindow = P(sim)$minAndMaxAgesToHarvest,
                                              harvestIndex = sim$harvestIndex
                                              )
      }

      sim <- scheduleEvent(sim, time(sim) + 1,  "simpleHarvest", "harvest")

      # ! ----- STOP EDITING ----- ! #
    },
    event2 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event


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
 
  if (!is.null(sim$areasToExclude) | !is.null(sim$DEMraster)) {
    if (!is.null(sim$areasToExclude)){
      
      protectedAreas <- spTransform(sim$areasToExclude, CRSobj = crs(sim$rasterToMatch)) %>%
        st_as_sf(.) %>%
        fasterize(sf = ., raster = sim$rasterToMatch, field = NULL, background = 0)
    }
  
    if (!is.null(sim$DEMraster)) {
      dem <- postProcess(sim$DEMraster, rasterToMatch = sim$rasterToMatch,
                         studyArea = sim$studyArea, useCache = TRUE, filename2 = NULL)
      dem[dem[] < P(sim)$ElevationToExclude] <- 0
      dem[dem[] >= P(sim)$ElevationToExclude] <- 1
      if (!is.null(protectedAreas)) {
        protVal <- getValues(protectedAreas)
        demVal <- getValues(dem)
        protectedAreas[] <-protVal + demVal
        protectedAreas[protectedAreas > 0] <- 1
        protectedAreas[is.na(dem)] <- NA
      }
      
    }
    sim$harvestExclusionRaster <- protectedAreas
  }
  
  #I dont' know if extract is faster than fasterize + extracting cell index for each unique value
  #Generate harvestLandscapeIndex needed for spread events
  sim$harvestAreas$newID <- as.numeric(row.names(sim$harvestAreas))
  harvestAreas <- sf::st_as_sf(sim$harvestAreas)
  harvestLand <- fasterize(sf = harvestAreas, raster = sim$rasterToMatch, field = 'newID')
  harvestIndex <- data.table("pixelIndex" = 1:ncell(harvestLand), "newID" = getValues(harvestLand))
  sim$harvestIndex <- harvestIndex[!is.na(newID)]
  
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

  return(invisible(sim))
}

harvestTrees <- function(pixelGroupMap, cohortData, exclusionAreas, harvestAreas, ageWindow, harvestIndex) {

  #Make an ageMap
  maxAges <- cohortData[, .(totalB = sum(B), age = max(age)), .(pixelGroup)]
  pixID <- data.table('pixelGroup' = getValues(pixelGroupMap), "pixelIndex" = 1:ncell(pixelGroupMap))
  #Join with pixel ID
  landStats <- maxAges[pixID, on = c("pixelGroup")]
  
  #I assume this method is faster than matching the pixelGroup raster values with a vector
  ageMap <- setValues(pixelGroupMap, landStats$age)
  
  mat <- c(-Inf, min(ageWindow), 1, min(ageWindow), max(ageWindow), 0, max(ageWindow), Inf, 1) %>%
    matrix(., ncol = 3, byrow = TRUE)
  ageMap <- reclassify(ageMap, rcl = mat)
  
  #harvestable ages are now 1
  
  #Now add ageMap and non-harvestable areas. 
  exlVals <- getValues(exclusionAreas)
  ageVals <- getValues(ageMap)
  nonharvest <- setValues(ageMap, c(exlVals + ageVals))
  
  #only zeros can be harvested
  spreadP <- nonharvest

  #build biomass table
  harvestTable <- data.table(id = harvestIndex$newId, index = harvestIndex$pixelIndex, biomass = landStats$totalB,
                  harvestStatus = getValues(nonharvest))
  
  
  numberOfCuts <- lapply(sim$harvestAreas, FUN = function(id, hT = harvestTable) {
    browser()
  })
  
  #Need some kind of mechanism to decide where to cut. Randomly sample without replacement... but how many spots
  cutAmounts <- lapply(sim$non)
  
  SpaDES.tools::spread2(landscape = nonharvest)
  spreadP[nonharvest[] > 0] <- 0
  spreadP[nonharvest = 0] <- 0.25 #This is the pertinent question. 
  
  #Need to use maxSize to make certain max of biomass and harvest size
  
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("studyArea", sim)) {
    message("study area not supplied. Defaulting to Ft. St. John TSA")
    sim$studyArea <- prepInputs(url = extractURL(objectName = "studyArea", sim = sim),
                                fun = "shapefile", 
                                destinationPath = tempdir(),
                                useCache = TRUE)
  }
  
  if (!suppliedElsewhere("rasterToMatch", sim)) {
    message("rasterToMatch not supplied. Defaulting to LCC2005. CRS will overide that of studyArea, if one was provided")
    sim$rasterToMatch <- prepInputsLCC(studyArea = sim$studyArea, filename2 = NULL, destinationPath = tempdir(), useCache = TRUE)
    sim$studyArea <- spTransform(sim$studyArea, CRSobj = crs(sim$rasterToMatch))
    sim$rasterToMatch <- mask(sim$rasterToMatch, sim$studyArea)
  }

  if (!suppliedElsewhere("cohortData", sim) & !suppliedElsewhere('pixelGroupMap', sim)) {
    message("either cohortData, pixelGroup, or both are not supplied. Using simulated data. Consider running the module LBMR")
    sim$cohortData <- data.table("pixelGroup" = c(1, 2), "speciesCode" = c("Pice_mar", "Pice_gla"), B = c(1000, 500), age = c(50, 10))
    temp <- sim$rasterToMatch
    temp[temp < 15] <- sample(c(1,2), size = length(temp[temp < 15]), replace = TRUE)
    temp[temp > 2] <- NA
    sim$pixelGroupMap <- temp
  }
  if (!suppliedElsewhere("harvestAreas", sim)) {
    message("harvestAreas not supplied. You should really supply harvestAreas so simulated units correspond to raster. Modifying studyArea instead.")
    harvestAreas <- sim$studyArea
    harvestAreas$AnnualBiomassHarvestTarget <- 1e5
    harvestAreas$MaximumAllowableCutSize <- 6.25 * 10 #62.5 hectares. This is for easy calculation with cell size
    harvestAreas$meanCutSize <- 6.25 * 5
    sim$harvestAreas <- harvestAreas
  }
  
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
