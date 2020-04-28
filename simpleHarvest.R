
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
                                           exclusionAreas = sim$harvestExclusionRaster,
                                           harvestAreas = sim$harvestAreas,
                                           ageWindow = P(sim)$minAndMaxAgesToHarvest,
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

  sim <- saveFiles(sim)

  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  
  plot(sim$rstCurrentHarvest)

  return(invisible(sim))
}

harvestSpreadInputs <- function(pixelGroupMap, cohortData, exclusionAreas, harvestAreas, ageWindow, harvestIndex) {

  #This function calculates which areas are available to cut, how many cuts of mean size would be needed to achieve the target
  # based on the mean biomass in a pixel.
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

  #Now add ageMap and non-harvestable areas
  ageVals <- getValues(ageMap)
  if (!is.null(exclusionAreas)) {
    exlVals <- getValues(exclusionAreas)
  } else {
    exlVals <- rep.int(0, times = length(ageVals))
  }

  harvest <- setValues(ageMap, c(exlVals + ageVals))
  #only zeros can be harvested

  #Remove NAs, build table of biomass, harvest status, and id of each pixel
  harvestStatus = getValues(harvest) %>%
    .[!is.na(.)]
  landStats <- landStats[!is.na(age)]
  harvestIndex <- harvestIndex[pixelIndex %in% landStats$pixelIndex]
  #build biomass table
  harvestTable <- data.table(newID = harvestIndex$newID, index = harvestIndex$pixelIndex, biomass = landStats$totalB,
                  harvestStatus = harvestStatus)
  
  #Join with the data in harvestAreas
  harvestGIS <- data.table(harvestAreas@data)
  harvestGIS <- harvestGIS[, .(newID, AnnualBiomassHarvestTarget, MaximumAllowableCutSize, meanCutSize)]
  harvestTable <- harvestTable[harvestGIS, on = c("newID")]
  
  #Find ha in a pixel and biomass conversion from g/m2 to tonnes/ha
  UCF <- 0.1 #Unit Conversion Factor: x * 1 g/m2 -> x * 10000 g/ha ->  x * 0.01 Mg/ha
  pixelSize <- prod(res(pixelGroupMap))/10000
  #Calculate the no. of harvest events to simulate in each harvestID (AnnualBiomassharvestTarget/meanB/pixel*meanCutsize)
  cutParams <- harvestTable[, .(pixelHarvestTarget = mean(AnnualBiomassHarvestTarget)/(mean(biomass) * UCF * pixelSize),
                                meanPixelsPerCut = mean(meanCutSize)/pixelSize,
                                maxCutPixels = round(mean(MaximumAllowableCutSize)/pixelSize, digits = 0),
                                meanPixelsPerCut = round(mean(meanCutSize/pixelSize), digits = 0)),
                                .(newID)]
  CutsToSimulate <- cutParams[, .(totalCuts = round(pixelHarvestTarget/meanPixelsPerCut), maxCutPixels, meanPixelsPerCut), .(newID)]
  
  
  #From harvestTable, sample pixels that are available to cut - could use focal to pick blobs but is it more trouble than worth?
  harvestLocs <- lapply(CutsToSimulate$newID, FUN = function(i, HarvestTable = harvestTable, CS = CutsToSimulate){
    HarvestTable <- HarvestTable[newID == i & harvestStatus == 0]
    CS <- CS[newID == i,]
    pix <- sample(x = HarvestTable$index, size = CutsToSimulate$totalCuts, replace = FALSE)
    return(pix)
    })

  names(harvestLocs) <- CutsToSimulate$newID
  
  #calculate harvest in each area separately - this is because maxCut may differ..
  harvestInEachPolygon <- lapply(names(harvestLocs), FUN = function(i, 
                                                                    loc = harvestLocs, 
                                                                    landscape = harvest, 
                                                                    CutTable = CutsToSimulate, 
                                                                    HarvestTable = harvestTable){
    initialPixels <- harvestLocs[i][[i]]
    CutTable <- CutTable[newID == i]
    HarvestTable = HarvestTable[newID == i]
    HarvestTable$harvestable <- 1
    HarvestTable[harvestStatus > 0]$harvestable <- 0
    
    #Calculate landscape properties - assume all forests have 8 valid neighbours for now 
    #TODO: fix this assumption
    propHarvest <- sum(HarvestTable$harvestable)/nrow(HarvestTable$harvestable)
    meanMaxRatio <- CutTable$meanPixelsPerCut/CutTable$maxCutPixels
    probLandscape <- setValues(landscape, NA)
    probLandscape[HarvestTable$index] <- HarvestTable$harvestable
    
    
    #If maxCuts > 9, this will require multiple spread iterations
    if (CutTable$maxCutPixels > 9) {
      #this is the cumulative distribution of getting the mean number of cells
      sampleProbs <- dbinom(CutTable$meanPixelsPerCut, size = CutTable$maxCutPixels, c(1:100)/100)
      bestProb <- c(1:100/100)[sampleProbs == max(sampleProbs)]
      randomMaxes <- rbinom(length(initialPixels), size = CutTable$maxCutPixels, prob = bestProb)
      
      #mean of randomMaxes should be close to CutTable$meanPixelsPerCut but will underestimate because not all neighbours harvestable
      #this is hard to control for because it depends on variegation of the harvestable cells
      harvest <- spread2(landscape = probLandscape, 
                    start = initialPixels, 
                    spreadProb = probLandscape, #this is 1, so harvest will always spread to random maxes
                    maxSize = randomMaxes,
                    asRaster = FALSE)
    } else {
      sampleProbs <- dbinom(CutTable$meanPixelsPerCut, size = CutTable$maxCutPixels, c(1:100)/100)
      bestProb <- c(1:100/100)[sampleProbs == max(sampleProbs)]
      probLandscape[probLandscape == 1] <- bestProb
      
      harvest <- spread2(landscape = probLandscape, 
                         start = initialPixels, 
                         spreadProb = probLandscape, #this is 1, so harvest will always spread to random maxes
                         maxSize = CutTable$maxCutPixels,
                         asRaster = FALSE)
    }
    return(harvest)
  })
  
  harvestValues <- do.call(rbind, harvestInEachPolygon)
  harvest[!is.na(harvest[])] <- 0
  harvest[harvestValues$pixels] <- 1
  
  return(harvest)
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
