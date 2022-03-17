
#### Download Tutorial Dataset and Setup Settings System ####
from <- "https://github.com/bpollner/data/raw/main/flowdexTutorial/flowdexTutorial.zip"
destination <- "~/desktop"
targetZip <- paste0(destination, "/flowdexTutorial.zip")
download.file(from, targetZip, mode="wb")
unzip(targetZip, exdir = destination)


### Initialise Settings File System ###
library(flowdex)
setupSettings(destination)
# possibly restart R



### Quckstart ###
# If you want to immediately see what flowdex can do:
destination <- "~/desktop" # again, as R was possibly restarted
setwd(paste0(destination, "/flowdexTutorial"))
library(flowdex)
fdmat <- flowdexit() # this might take a few seconds
fdmat@pData # to inspect volume and sample ID data
fdmat@cyTags # to inspect class- and numerical variables assigned to each sample
fdmat[[1]] # to inspect fluorescence distribution



#### Create home-folder and folder structure ####
destination <- "~/desktop" # again, as R was possibly restarted
expHome <- paste0(destination, "/tapWater@home")
dir.create(expHome)
setwd(expHome)
library(flowdex) # as you might have had to restart R above
genfs() # is generating the required folder structure within "tapWater@home"

## now would be data acquisition

# move or copy fcs files
from <- list.files(paste0(destination, "/flowdexTutorial/fcsFiles"), full.names = TRUE)
to <- paste0(destination, "/tapWater@home/fcsFiles")
file.copy(from, to, overwrite = TRUE)

# move or copy dictionary
from <- list.files(paste0(destination, "/flowdexTutorial/dictionary"), full.names = TRUE)
to <- paste0(destination, "/tapWater@home/dictionary")
file.copy(from, to, overwrite = TRUE)


#### Workflow A: Define Gating Strategy and Polygon Gates ####

## make Gating Set and plot ##
gsAll <- makeGatingSet()
gsRed <- makeGatingSet(patt="GPos_T6_th2")
gsRed
plotgates(gsRed, toPdf = FALSE, x="FITC.A", y="PerCP.A")

## draw a gate ##
drawGate(gsRed, flf=5, gn="root", pggId="BactStainV1", channels=".")
# point and click to draw a polygon gate around the population of interest

# draw gate again
drawGate(gsRed, flf=5, gn="root", pggId="BactStainV1", show="BactStainV1")

# copy BactStainV1 from tutorial
from <- paste0(destination, "/flowdexTutorial/gating/BactStainV1")
to <- paste0(destination, "/tapWater@home/gating")
file.copy(from, to, overwrite = TRUE)

# copy gating strategy template
from <- paste0(destination, "/tapWater@home/templates/gateStrat.xlsx")
to <- paste0(destination, "/tapWater@home/gating")
file.copy(from, to)

# copy the gating strategy
from <- paste0(destination, "/flowdexTutorial/gating/gateStrat.xlsx")
to <- paste0(destination, "/tapWater@home/gating")
file.copy(from, to, overwrite=TRUE)

# ad gate to gsRed
gsRed_ga <- addGates(gsRed)
flowWorkspace::plot(gsRed_ga) # to view the gate hierarchy
plotgates(gsRed_ga, toPdf = FALSE) # to view the gated data


# second iteration: create a nested gate
drawGate(gsRed_ga, flf=5, gn="DNA+", pggId="pg2", channels = c("FITC.A", "SSC.A"))
# point and click to draw a polygon gate around the population of interest

# add the gate in the gating strategy file:
# GateName: FooGate
# Parent: DNA+
# GateOnX: FITC.A
# GateOnY: SSC.A
# GateDefinition: pg2
# extractOn: SSC.A
# minRange: 0
# maxRange: 4000
# keepData: TRUE

gsRed_ga <- addGates(gsRed_ga) # add the gate to the gating set
flowWorkspace::plot(gsRed_ga) # visualise gating hierarchy
plotgates(gsRed_ga, toPdf = F)

## copy gating strategy 2
from <- list.files(paste0(destination, "/flowdexTutorial/gating"), full.names = TRUE)
to <- paste0(destination, "/tapWater@home/gating")
file.copy(from, to, overwrite=TRUE) # this might overwrite the polygon gate definition you created above

# and make again reduced gating set
gsRed2 <- makeAddGatingSet(patt="GPos_T6_th2", gateStrat = "gateStrat_2")
flowWorkspace::plot(gsRed2) # to view the gating hierarchy
plotgates(gsRed2, toPdf = F) # to view the gated data





#### Workflow B: Extract Fluorescence Distributions and Visualize ####
# copy all relevant files from the flowdex tutorial
destination <- "~/desktop"
fromFcs <- list.files(paste0(destination, "/flowdexTutorial/fcsFiles"), full.names=TRUE)
fromDict <- list.files(paste0(destination, "/flowdexTutorial/dictionary"), full.names=TRUE)
fromGating <- list.files(paste0(destination, "/flowdexTutorial/gating"), full.names=TRUE)
toFcs <- paste0(destination, "/tapWater@home/fcsFiles")
toDict <- paste0(destination, "/tapWater@home/dictionary")
toGating <- paste0(destination, "/tapWater@home/gating")
#
file.copy(fromFcs, toFcs, overwrite=TRUE)
file.copy(fromDict, toDict, overwrite=TRUE)
file.copy(fromGating, toGating, overwrite=TRUE)


# if all requirements are met, we can call
fdmat1 <- flowdexit()
# and to inspect the result:
fdmat1 # displaying cyTags and sample IDs (in the slot pDAta)
fdmat1[[1]] # displaying data from the first (and only) gate

# access gating set
gs1 <- gsenv$gatingSet
# inspect the gating set:
gs1
flowWorkspace::plot(gs1)


fdmat2 <- flowdexit(gateStrat = "gateStrat_2")
# and to inspect the result:
fdmat2 # displaying cyTags and sample IDs (in the slot pDAta)
fdmat2[[1]] # displaying data from the first gate
fdmat2[[2]] # displaying data from the second gate

# access gating set
gs2 <- gsenv$gatingSet
# inspect the gating set:
gs2
flowWorkspace::plot(gs2)



# visualise gated data
plotgates(gs1, toPdf = FALSE) # a bit much on one graphic

# use split in plotgates
colnames(fdmat1@cyTags)
plotgates(gs1, spl="Y_Time.d", toPdf = TRUE)

# read in fcs files from only a single day and split by C_treatment
# to obtain a more use-friendly graphic
gs_d4 <- makeAddGatingSet(patt = "T4", gateStrat = "gateStrat_2")
plotgates(gs_d4, spl="C_treatment", fns="_day4")
gs_d5 <- makeAddGatingSet(patt = "T5", gateStrat = "gateStrat_2")
plotgates(gs_d5, spl="C_treatment", fns="_day5")
gs_d6 <- makeAddGatingSet(patt = "T6", gateStrat = "gateStrat_2")
plotgates(gs_d6, spl="C_treatment", fns="_day6")

# and plot all gates:
plotgates(gs_d6, spl="C_treatment", fns="_day6_allGates", plotAll = TRUE)



# visualize fluorescence distributions
plotFlscDist(fdmat1)


# create a subset of data
fdmat_s <- flowdexit(patt = "T4_th1")
fdmat_s # has now only 12 samples
gs_s <- gsenv$gatingSet # not required, just to keep it
#
plotFlscDist(fdmat_s, toPdf = FALSE) # much nicer

colnames(fdmat_s@cyTags)
plotFlscDist(fdmat_s, spl = "C_treatment",  ti="Day 4, first third",  toPdf = FALSE)
plotFlscDist(fdmat_s, spl = "C_treatment", fns="_d4_th1", ti="Day 4, first third") # export to pdf




#### Accessory Functions ####


## Check and Repair fcs files ##
# copy erroneous fcs files to experiment home
destination <- "~/desktop"
from <- list.files(paste0(destination, "/flowdexTutorial/fcsF_E_rep"), full.names=TRUE)
to <- paste0(destination, "/tapWater@home/fcsF_E_rep")
dir.create(to)
file.copy(from, to, overwrite=TRUE)

# check fcs files
checkRepairFcsFiles(fn="fcsF_E_rep")

# and repair
checkRepairFcsFiles(fn="fcsF_E_rep", fcsRepair = TRUE, confirm = FALSE)
# check again, all should be good now:
checkRepairFcsFiles(fn="fcsF_E_rep")



## Repair Volumes ##
# copy erroneous fcs files to experiment home
destination <- "~/desktop"
from <- list.files(paste0(destination, "/flowdexTutorial/fcsF_E_vol_sid"), full.names=TRUE)
to <- paste0(destination, "/tapWater@home/fcsF_E_vol_sid")
dir.create(to)
file.copy(from, to, overwrite=TRUE)

# repair volume data
repairVolumes(fn = "fcsF_E_vol_sid", vol=1234567)
# press enter to confirm
# and check again:
repairVolumes(fn = "fcsF_E_vol_sid", vol=1234567) # all should be good
file.copy(from, to, overwrite=TRUE) # restore the "bad" files
repairVolumes(fn = "fcsF_E_vol_sid", vol=1234567, confirm=FALSE)
# to run it without having to confirm

# force all to have the same volume data
repairVolumes(fn = "fcsF_E_vol_sid", vol=1010101, confirm=FALSE, includeAll = TRUE)



## repair SID ##
flowset <- repairSID(fn = "fcsF_E_vol_sid")
flowset@phenoData@data # very bad sample ID in the fourth sample
# view the correct sample IDs of the other samples
# copy one of those correct sample IDs
# paste and modify it - it should be beaker #3:
nsid <- "tr: GPos; Td: 5; wt: nativ; ap: no; th: th1; ha: ha1; bk: b3"
# also copy and paste the sample name
sana <- "N_na_GPos_T5_th1_b3.fcs" # the  name of the sample having the faulty sample ID
# now put all together and write fcs file with correct sample ID back to disk
repairSID(fs=flowset, fn="fcsF_E_vol_sid", name=sana, newSID = nsid)
# press enter to confirm
#
# and check again:
flowset <- repairSID(fn = "fcsF_E_vol_sid")
flowset@phenoData@data # all is good



## Apply Bandpass Filter ##
fdmat_s <- flowdexit(patt = "T4_th1")
plotFlscDist(fdmat_s, toPdf = FALSE)

fdmat_s_bp <- applyBandpass(fdmat_s, bandpass = c(1600, 2400))
fdmat_s[[1]] # compare
fdmat_s_bp[[1]] #
ncol(fdmat_s[[1]])
ncol(fdmat_s_bp[[1]])

# visualize difference
plotFlscDist(fdmat_s_bp, toPdf = FALSE)

# export rawdata with applied bandpass
exportFdmatData(fdmat_s_bp, expo.name = "flscData_d4_th1")



