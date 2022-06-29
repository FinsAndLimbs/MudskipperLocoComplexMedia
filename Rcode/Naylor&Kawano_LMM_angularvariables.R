################################################################################
##  Naylor& Kawano - Mudskipper locomotion on complex media: angular variables ##
##                                                              June 13, 2022  ##
################################################################################
#Adapted from 'mudskipperKinematics' repo (https://github.com/FinsAndLimbs/mudskipperKinematics)
  #and code from Z. Quigley
      #note that some terms reflect salamander anatomy


rm(list = ls()) #clear environment 


# 1) Load packages and data -----------------------------------------------

library(devtools)    # for install_github()
library(signal)      # for interp1()
library(tidyverse)   # for summarize()
library(reshape2)    # for melt()
library(grid)        # for textGrob()
library(gridExtra)   # for arrangeGrob() and grid.arrange()
library(stats)
library(lme4)
library(performance)
library(emmeans)
library(MuMIn)

## use devtools to load the kraken repo from GitHub
#install_github("MorphoFun/kraken", dependencies = TRUE)

##load kraken repo as a library to use the repo functions
library(kraken)



# 2) Prepare raw 2D coordinates for MATLAB 3D transformation ----------------------------------------------------

#***(skip this step if already using 3D coordinate files)***

#a) Gather files, limit and rename columns

# Choose the directory containing the xypts files exported from DLTdv
xy_path <- choose.dir(caption = "raw_DLTdv8_output")
#list files
xy_list <- list.files(xy_path, pattern="xypts.csv", full=TRUE)

# Name matrices in array with trial name
# e.g., FS_10l_pb06f01_d
Trial <- list(substring(basename(xy_list), 1, 16))

# identify treatment group
Group <- list(substring(basename(xy_list), 1, 6))

# Get data for 5 anatomical landmarks 
## The order of landmarks in the xypts.csv is: 
# 1. pectoral fin tip
# 2. intra-fin joint ("elbow")
# 3. shoulder
# 4. midline point 1 (anterior dorsal fin)
# 5. midline point 2 (above shoulder)
# -> limit the dataset to the first 10 columns (5 landmarks) 
  #& remove NaN rows (i.e.,frames that are not digitized within videos)
  #rename columns for compatibility with MATLAB code (Butcher & Blob 2008)
Pb_RawFiles <- array(lapply(xy_list, read.csv, header=T, na.strings = c("NaN")), dimnames=Trial)
for (j in 1:length(Pb_RawFiles)) {
  colnames(Pb_RawFiles[[j]]) <- c("ToeTip_X", "ToeTip_Y", "KneeElbow_X", "KneeElbow_Y", "HipShoulder_X", "HipShoulder_Y",
                                   "dorsal1_X", "dorsal1_Y", "dorsal2_X", "dorsal2_Y")
}

for (i in 1:length(Pb_RawFiles)) {
  Pb_RawFiles[[i]] <- Pb_RawFiles[[i]][,1:10]
  Pb_RawFiles[[i]] <- na.omit(Pb_RawFiles[[i]])
}


#b) Add dummy columns for compatibility with MATLAB code

## Order of landmarks for putting into MATLAB: 
# 1. tip of toe (fin)
# 2. Elbow (intra-fin joint)
# 3. Shoulder
# 4. (Wrist)                     ---------------> not tracked for this project
# 5. (Metacarpophalangeal joint) ---------------> not tracked for this project
# 6. Midline 1 
# 7. Midline 2 

# add dummy columns as place holders for the wrist and MP points
# the dummy columns must have numeric data (the files will NOT run if filled with NA's)
Pb_RawFiles_MATLAB <- lapply(Pb_RawFiles, FUN = function(x) data.frame(x[,1:6], 0, 0, 0, 0, x[,7:10]))
for (j in 1:length(Pb_RawFiles_MATLAB)) {
  colnames(Pb_RawFiles_MATLAB[[j]]) <- c("ToeTip_X", "ToeTip_Y", "KneeElbow_X", "KneeElbow_Y", "HipShoulder_X", "HipShoulder_Y",
                                         "AnkleWrist_X", "AnkleWrist_Y", "MetaPhalJoint_X", "MetaPhalJoint_Y",
                                         "dorsal1_X", "dorsal1_Y", "dorsal2_X", "dorsal2_Y")
}

# c) Save each file as a separate .txt file to be analyzed with MATLAB code

RawFiles_MATLAB_Path <- paste(dirname(xy_path), "/Step2_RawFiles", sep = "")
lapply(1:length(Pb_RawFiles_MATLAB), FUN = function(x) write.table(data.frame(Pb_RawFiles_MATLAB[[x]]), file = paste(RawFiles_MATLAB_Path, "/", names(Pb_RawFiles_MATLAB[x]), "_Raw.txt", sep = ""), sep = "\t", row.names = FALSE))


## -> MATLAB code will return .csv's with xyz (3D) coordinates




# 3) Prepare 3D coordinate data for analyses ----------------------------------------------------

#### a) Read in XYZ (3D) coordinates from MATLAB ####


#Choose the directory containing the _XYZ.csv files
#
pb_kin_path <- choose.dir(caption = "MATLAB_output_XYZ")
#list files
pb_kin_list <- list.files(pb_kin_path, pattern="XYZ.csv", full=TRUE)

##Landmark order for XYZ files outputted from MATLAB* 
  #*this is different than the inputted order
#Columns:
# 1:3 =   anterior midline point
# 4:6 =   posterior midline point
# 7:9 =   shoulder
# 10:12 = intra-fin
# 13:15 = DUMMY
# 16:18 = DUMMY
# 19:21 = pectoral fin tip 


#Name matrices in array with trial name (e.g., FS_001_pb06f01)
# remove the directory info from the directory path*
  #this duplicates some of the steps taken before, 
  #but it helps ensure that we correctly match kinematic files with trial names
pb_kin_filenames <- basename(pb_kin_list)
pb_kin_trial <- substring(pb_kin_filenames, 1, 14)


# Get data 
kinFiles_pb_raw <- array(lapply(pb_kin_list, read.csv, header=F))
#assign trial names
names(kinFiles_pb_raw) <- pb_kin_trial



#### b) Interpolate data to 101 points ####


# Standardize the number of frames per trial via interpolation
# interpolate to 101 points; each frame represents 1% (from 0 to 100%) of stride
  #use kraken::interpolateR() function

# apply interpolation to all trials            
Pb_kine_interp <- lapply(kinFiles_pb_raw, FUN = function(x) interpolateR(x, 101))




# 4) Calculate angular variables --------------------------------------------------------

#Angular variables:
# 1) intra-fin angle
# 2) protraction/retraction of radials (about the shoulder)
# 3) abduction/adduction of radials (about the shoulder)
# 4) yaw (of the body)
# 5) pitch (of the body)

#Variable abbreviations:
#elbow: intra-fin
#proret: protraction / retraction
#abad: abduction / adduction

#modified code from Z. Quigley
anglecalc <- function(trial){
  #Single function to calculate angles
  #Trial data must have landmark XYZ coordinates for columns and frame numbers for rows
  #(col1 = 1x, col2 = 1y, col3= 1z, col4 = 2x ... col21 = 7z)
  #(row1 = frame 1, row2 = frame 2 ... row101 = frame 101)
 
 #Use jointAngle() from kraken
  elbow <- jointAngle(trial[,  19:21], trial[, 10:12], trial[, 7:9])
  
 #Use yawAngle() from kraken
  yaw_reg <- yawAngle(trial[, 4:6], trial[, 1:3])
  #get absolute value for magnitude of yaw
  yaw <- abs(yaw_reg)
  
 #Protraction()** from kraken
  #*use modified function defined later in the script
  proret <- protraction_updated(trial[, 7:9], trial[, 10:12], yaw_reg)
  
 #Abduction code based on Butcher & Blob (2008) MATLAB code - see updated below
  #multiply by -1 as a correction (for coordinate system)
  abad <- abduction(trial[, 7:9], trial[, 10:12], yaw_reg)*-1
  
 #Pitch** code based on Butcher & Blob (2008) MATLAB code
  #**use modified function defined later in the script
  pitch <- pitch3D(trial)

 #create data frame with calculations  
  return(data.frame(elbow, yaw, proret, abad, pitch))
}

 #**updated protraction function created by Z. Quigley
 #(based on Butcher & Blob 2008 MATLAB code)
protraction_updated <- function(P1,P2, Yaw, ...) {
  HipPoint <- P1
  KneePoint <- P2
  
  FemurVector <- HipPoint - KneePoint
  HipPointTrans <- HipPoint - HipPoint
  KneePointTrans1A <- KneePoint - HipPoint 
  #@KneePointTrans1B <- KneePoint - HipPoint
  FemurVectorTrans1A <- HipPointTrans - KneePointTrans1A
  #@FemurVectorTrans1B <- HipPointTrans - KneePointTrans1B
  
  #@TibiaVector <- KneePoint - AnklePoint
  #@KneePointTrans2A <- KneePoint - KneePoint
  #@AnklePointTransA <- AnklePoint - KneePoint
  #@KneePointTrans2B <- KneePoint - KneePoint
  #@AnklePointTransB <- AnklePoint-KneePoint
  #@KneePointTrans2C <- KneePoint-KneePoint
  #@AnklePointTransC <- AnklePoint-KneePoint
  #@TibiaVectorTransA <- KneePointTrans2A - AnklePointTransA
  #@TibiaVectorTransB <- KneePointTrans2B-AnklePointTransB
  #@TibiaVectorTransC <- KneePointTrans2C-AnklePointTransC
  
  ## dot product between two row vectors
  wdot <- function(a, b) {
    y <- a*b
    y <- t(sum(t(y)))
    return(y)
  }
  
  ## Cross product between two row vectors 
  #library(pracma) # for the cross() function that returns a vector
  wcross <- function(a, b) {
    c <- t(pracma::cross(t(a),t(b)))
    return(c)
  }
  
  ## vlength
  vlength <- function(x) {
    v <- (wdot(x, x))^0.5
    return(v)
  }
  
  ## Setting up Transverse Plane
  TVPlane1 <- c(0, 0, 0)
  TVPlane2 <- c(0, .1, 0)
  TVPlane3 <- c(0, 0, .1)
  TVVector1 <- TVPlane2 - TVPlane1
  TVVector2 <- TVPlane3 - TVPlane1
  # had to take the transpose to get it into row-form instead of column form
  TVNorm <- t(wcross(TVVector1,TVVector2))
  
  FemTVAngInit <- matrix()
  for (i in 1:nrow(KneePointTrans1A)) {
    FemurVectorTrans1AA <- FemurVectorTrans1A[i,]
    dotFemurVectorTV <- wdot(FemurVectorTrans1AA,TVNorm)
    MagFemurVectorTV <- vlength(FemurVectorTrans1AA)
    MagTVNorm <- vlength(TVNorm)
    MagFemurVectorTVNorm <- MagFemurVectorTV*MagTVNorm
    CosdotFemurVectorTV <- dotFemurVectorTV/MagFemurVectorTVNorm
    FemTVAngA <- (acos(CosdotFemurVectorTV))*(180/pi) #converting radians to degree
    # %makes a one column matrix of angles of the protraction/retraction angle
    FemTVAngInit[i] <- FemTVAngA
  }
  
  FemTVAng <- FemTVAngInit
  #setting the zero of the angles to be perpendicular to x axis
  FemTVAng <- (90-FemTVAng)
  FemTVAng <- (90-FemTVAng)+Yaw-90
  
  return(FemTVAng)
}

 #*updated abduction function created by Z. Quigley
 #(based on Butcher & Blob 2008 MATLAB code)
abduction <- function(prox, dist, m1, m2){
  
  #input variables are n by 3 dataframes containing the X, Y, Z coordinates for 
  #the proximal joint (shoulder/ hip), distal joint (elbow/knee), and two points 
  ##along the midline, where n equals the number of digitized trial frames
  
  
  
  HipPoint <- prox
  KneePoint <- dist
  
  FemurVector <- HipPoint - KneePoint
  HipPointTrans <- HipPoint - HipPoint
  KneePointTrans1A <- KneePoint - HipPoint 
  FemurVectorTrans1A <- HipPointTrans - KneePointTrans1A
  
  ## dot product between two row vectors
  wdot <- function(a, b) {
    y <- a*b
    y <- t(sum(t(y)))
    return(y)
  }
  
  ## Cross product between two row vectors 
  #library(pracma) # for the cross() function that returns a vector
  wcross <- function(a, b) {
    c <- t(pracma::cross(t(a),t(b)))
    return(c)
  }
  
  ## vlength
  vlength <- function(x) {
    v <- (wdot(x, x))^0.5
    return(v)
  }
  
  #Generate frontal plane at the height of the proximal joint by moving m1 and m2 
  #to height of prox, then generating plane beween points
  ## Setting up frontal Plane
  planeP1 <- c(0, 0, 0)
  #planeP1[, 1] <-  prox[, 1]#set m1 to height of prox
  planeP2 <- c(.1, 0, 0)
  #planeP2[, 1] <-  prox[, 1]#set m2 to height of prox
  planeP3 <- c(0, .1, 0)
  planeVec1 <- planeP2 - planeP1
  planeVec2 <- planeP3 - planeP1
  PlaneNorm <- t(wcross(planeVec1,planeVec2))#normal to plane
  # had to take the transpose to get it into row-form instead of column form
  
  
  FemFRAngInit <- matrix()
  for (i in 1:nrow(KneePointTrans1A)) {
    FemurVectorTrans1AA <- FemurVectorTrans1A[i,]
    dotFemurVectorFR <- wdot(FemurVectorTrans1AA,PlaneNorm)
    MagFemurVectorFR <- vlength(FemurVectorTrans1AA)
    MagPlaneNorm <- vlength(PlaneNorm)
    MagFemurVectorPlaneNorm <- MagFemurVectorFR*MagPlaneNorm
    CosdotFemurVectorFR <- dotFemurVectorFR/MagFemurVectorPlaneNorm
    FemFRAngA <- (acos(CosdotFemurVectorFR))*(180/pi) #converting radians to degree
    # %makes a one column matrix of angles of the abduction/adduction angle
    FemFRAngInit[i] <- FemFRAngA
  }
  
  FemFRAng <- abs(FemFRAngInit)
  #setting the zero of the angles to be perpendicular to x axis
  FemFRAng <- -(90-FemFRAng) 
  #I double checked my calculations with mudskipper data that was also processed
  #in MATLAB, and I was getting the proper magnitudes, but they were all positive
  #instead of negative? I added a "-" to the last line to correct
  
  #could correct for roll, but not sure how to calculate and Rick didn't do it 
  #in his code. Roll looks pretty negligible anyways
  #Roll <- 0
  # FemFRAng <- (90-FemFRAng) -90
  
  return(FemFRAng)
}

 #updated pitch function created by Z. Quigley
pitch3D <- function(trial){
  p1 <- trial[, 4:6] #posterior midpoint 
  p2 <- trial[, 1:3] #anterior midpoint 
  p3 <- cbind(p2[,1], p2[, 2], p1[, 3]) #generate third point at height of p1, directly above/below p2 (depending on +/- pitch)
  #
  #        p2
  #       /
  #      /
  #     /
  #    /
  #   /
  #p1/_____p3
  # (points form right triangle, and we'll calculate angle defined by (p3,p1,p2))
  #calculate differences (vectors) for xy and z coordinates
  dx21 <- p2[,1]-p1[,1] #dX p1 and p2
  dy21 <- p2[,2]-p1[,2] #dY p1 and p2
  dx31 <- p3[,1]-p1[,1] #dX p1 and p3
  dy31 <- p3[,2]-p1[,2] #dY p1 and p3
  dz   <- p2[,3]-p1[,3] #dZ p1 and p2
  pitch <- vector()
  for(i in c(1:nrow(p1))){
    pitch <- c(pitch, atan2(dz[i], sqrt(dx21[i]*dx31[i] + dy21[i]*dy31[i]))*180/pi)
  }
  return(pitch)
}



###Apply all calculations          
Pb_kine_angles <-lapply(Pb_kine_interp, anglecalc)



# 5) Extract stride velocity & compare across conditions --------------------------------------------------------


#### a) Use MorphoFun function* to obtain stride speed ####
 #*(modified for 3D coordinates) 
 
CrudeSpd <- function(x, rate, calib=1, distance.scale="units", time.scale="s") {
  xyz <- matrix(NA, 1, 5) #create empty matrix
  xyz[,1] <- nrow(x)*(1/rate) # first column = time elapsed based on frame rate & number of frames in trial
  xyz[,2] <- sqrt( ((x[nrow(x),1]-x[1,1])^2)+((x[nrow(x),2]-x[1,2])^2) + ((x[nrow(x),3]-x[1,3])^2) ) #column 2 = 3D displacement
  xyz[,3] <- xyz[,2]/calib #column 3 = calibrated 3D displacement (*does not apply for our data; already calibrated in MATLAB)
  xyz[,4] <- xyz[,3]/xyz[,1] #column 4 = distance/time => velocity
  xyz[,5] <- abs(xyz[,4]) #column 5 = absolute value of velocity => speed
  xyz <- data.frame(xyz) #turn matrix into data frame with column labels
  names(xyz) <- c(paste("Time",time.scale, sep="."), "Distance.Uncalib", paste("Distance",distance.scale,sep="."), 
                 paste("Velocity",distance.scale,time.scale, sep="."), paste("Speed",distance.scale,time.scale, sep="."))
  return(xyz)
}
#frame rate
rate <- 100

#Add trial names to speed values and create a data frame
Pb_kine_spd <- lapply(kinFiles_pb_raw, FUN = function(x) CrudeSpd(x, rate = rate))
Pb_kine_spd_df <- data.frame(data.matrix(do.call(rbind, Pb_kine_spd)))
Pb_kine_spd_df$trial <- names(Pb_kine_spd)
row.names(Pb_kine_spd_df) <- NULL


#### b) Test if stride speed is different across conditions ####

#Subset data based on substrate and incline
Pb_kine_spd_df$substrate <- substring(Pb_kine_spd_df$trial, 1, 2)
Pb_kine_spd_df$incline <- substring(Pb_kine_spd_df$trial, 4, 5)
Pb_kine_spd_df$ind <- substring(Pb_kine_spd_df$trial, 8, 11)

#Run linear-mixed effect model (LMM) with main fixed effect (enviro condition), and individual as random effect
  #and examine estimated marginal means (EMMs) - see Section #7 for details on approach
pb_kine_lmer_spd <- lmer(Speed.units.s ~ substrate*incline + (1|ind), data = Pb_kine_spd_df)
#calculate EMMs
emmeans(pb_kine_lmer_spd, ~substrate*incline) #speed is different between conditions



# 6) Prepare dataset for analysis --------------------------------------------------------


#### a) Correct pitch angles for conditions with 10 and 20 degree inclines ####

#Subtract incline values from pitch angles (column 5)
for(i in 1:length(Pb_kine_angles)){
  
  Pb_kine_angles[[i]][5] <- Pb_kine_angles[[i]][5] - as.numeric(Pb_kine_spd_df$incline[[i]]) 
  
}


#### b) Extract maximum, minimum, and mean values for each angle from each focal stride ####

Pb_kine_angles_max <- NULL
Pb_kine_angles_min <- NULL
Pb_kine_angles_mean <- NULL

for (i in 1:length(Pb_kine_angles)) {
  Pb_kine_angles_max[[i]] <- lapply(Pb_kine_angles[[i]], FUN = function(x) max(x))
}
Pb_kine_angles_max_df <- do.call(rbind.data.frame, Pb_kine_angles_max)
Pb_kine_angles_max_df$trial <- names(Pb_kine_angles)

for (i in 1:length(Pb_kine_angles)) {
  Pb_kine_angles_min[[i]] <- lapply(Pb_kine_angles[[i]], FUN = function(x) min(x))
}
Pb_kine_angles_min_df <- do.call(rbind.data.frame, Pb_kine_angles_min)
Pb_kine_angles_min_df$trial <- names(Pb_kine_angles)

for (i in 1:length(Pb_kine_angles)) {
  Pb_kine_angles_mean[[i]] <- lapply(Pb_kine_angles[[i]], FUN = function(x) mean(x))
}
Pb_kine_angles_mean_df <- do.call(rbind.data.frame, Pb_kine_angles_mean)
Pb_kine_angles_mean_df$trial <- names(Pb_kine_angles)


#### c) Combine extracted angles, speed, and other covariates ####

###Combine angles with speed data

Pb_kine_angles_max_spd <- merge(Pb_kine_spd_df, Pb_kine_angles_max_df, by="trial" )
Pb_kine_angles_min_spd <- merge(Pb_kine_spd_df, Pb_kine_angles_min_df, by="trial" )
Pb_kine_angles_mean_spd <- merge(Pb_kine_spd_df, Pb_kine_angles_mean_df, by="trial" )


###Combine these data frames with categorical covariate data (e.g., tail/body use, premovement)

## read in .csv for contact variables
contactvarsdf <- read.csv(choose.files(caption = "contactvariable data file"))
#rename first column of contactvariables df to match angles+speed df
names(contactvarsdf)[1] <- "trial"

#trim down contactvarsdf to just the 101 trials with angular data
Pb_timing_angles_files <- contactvarsdf[which(contactvarsdf$trial %in% Pb_kine_angles_max_spd$trial), ]
Pb_timing_angles_files

#merge contact df with max, mean, and min angular variable dfs
Pb_kine_angles_max_spd_categ <- merge(Pb_kine_angles_max_spd, Pb_timing_angles_files, by=c("trial"))
Pb_kine_angles_max_spd_categ

Pb_kine_angles_mean_spd_categ <- merge(Pb_kine_angles_mean_spd, Pb_timing_angles_files, by=c("trial"))
Pb_kine_angles_mean_spd_categ

Pb_kine_angles_min_spd_categ <- merge(Pb_kine_angles_min_spd, Pb_timing_angles_files, by=c("trial"))
Pb_kine_angles_min_spd_categ



# 7) Find best model --------------------------------------------------------


#model averaging approach:
#compare linear mixed-effects models (LMMs) with different covariate combinations
#base structure: outcome variable ~ substrate*incline + (1|Indiv)
#REML = FALSE (maximum likelihood)
#covariates: Speed.units.s (stride speed; continuous), tail.body.use (binary), 
# Premovement (binary)


#a) Run LMMs for each outcome variable (max, mean, and minimum) 
   

##### intra-fin angle ####

## Max


#no additional fixed terms/covariates
intrafin_max0 <- lmer(elbow ~ substrate*incline + (1|ind), 
                   data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#one additional fixed term/covariate
intrafin_max1a <- lmer(elbow ~ substrate*incline + Speed.units.s + (1|ind), 
                   data = Pb_kine_angles_max_spd_categ, REML=FALSE)

intrafin_max1b <- lmer(elbow ~ substrate*incline + tail.body.use + (1|ind), 
                    data = Pb_kine_angles_max_spd_categ, REML=FALSE)

intrafin_max1c <- lmer(elbow ~ substrate*incline + Premovement + (1|ind), 
                    data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
intrafin_max2a <- lmer(elbow ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                    data = Pb_kine_angles_max_spd_categ, REML=FALSE)

intrafin_max2b <- lmer(elbow ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                    data = Pb_kine_angles_max_spd_categ, REML=FALSE)

intrafin_max2c <- lmer(elbow ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                    data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
intrafin_max3 <- lmer(elbow ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                   data = Pb_kine_angles_max_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(intrafin_max0, intrafin_max1a, intrafin_max1b, intrafin_max1c, intrafin_max2a, intrafin_max2b, intrafin_max2c, intrafin_max3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

#2c, 1b, 3, 1c, 2a have delta AICc< 2; select 1b

anova(intrafin_max1b, intrafin_max2c)
anova(intrafin_max1b, intrafin_max3)
anova(intrafin_max1b, intrafin_max2a)

anova(intrafin_max1c, intrafin_max2a)
anova(intrafin_max1c, intrafin_max2c)
anova(intrafin_max1c, intrafin_max3)

anova(intrafin_max2a, intrafin_max3)
anova(intrafin_max2c, intrafin_max3)



## Min


#no additional fixed terms/covariates
intrafin_min0 <- lmer(elbow ~ substrate*incline + (1|ind), 
                   data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#one additional fixed term/covariate
intrafin_min1a <- lmer(elbow ~ substrate*incline + Speed.units.s + (1|ind), 
                   data = Pb_kine_angles_min_spd_categ, REML=FALSE)

intrafin_min1b <- lmer(elbow ~ substrate*incline + tail.body.use + (1|ind), 
                    data = Pb_kine_angles_min_spd_categ, REML=FALSE)

intrafin_min1c <- lmer(elbow ~ substrate*incline + Premovement + (1|ind), 
                    data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
intrafin_min2a <- lmer(elbow ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                    data = Pb_kine_angles_min_spd_categ, REML=FALSE)

intrafin_min2b <- lmer(elbow ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                    data = Pb_kine_angles_min_spd_categ, REML=FALSE)

intrafin_min2c <- lmer(elbow ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                    data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
intrafin_min3 <- lmer(elbow ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                   data = Pb_kine_angles_min_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(intrafin_min0, intrafin_min1a, intrafin_min1b, intrafin_min1c, intrafin_min2a, intrafin_min2b, intrafin_min2c, intrafin_min3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 0, 1c, 1b, 2c have delta AICc< 2; select 0

anova(intrafin_min0,intrafin_min1b) 
anova(intrafin_min0,intrafin_min1c)
anova(intrafin_min0, intrafin_min2c)

anova(intrafin_min1b, intrafin_min2c) 

anova(intrafin_min1c, intrafin_min2c)



## Mean


#no additional fixed terms/covariates
intrafin_mean0 <- lmer(elbow ~ substrate*incline + (1|ind), 
                    data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#1 additional fixed term/covariate
intrafin_mean1a <- lmer(elbow ~ substrate*incline + Speed.units.s + (1|ind), 
                    data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

intrafin_mean1b <- lmer(elbow ~ substrate*incline + tail.body.use + (1|ind), 
                     data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

intrafin_mean1c <- lmer(elbow ~ substrate*incline + Premovement + (1|ind), 
                     data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#2 additional fixed terms/covariates
intrafin_mean2a <- lmer(elbow ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                     data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

intrafin_mean2b <- lmer(elbow ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                     data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

intrafin_mean2c <- lmer(elbow ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                     data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#3 additional fixed terms/covariates
intrafin_mean3 <- lmer(elbow ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                    data = Pb_kine_angles_mean_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(intrafin_mean0, intrafin_mean1a, intrafin_mean1b, intrafin_mean1c, intrafin_mean2a, intrafin_mean2b, intrafin_mean2c, intrafin_mean3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 1b has delta AICc< 2; select 1b



##### protraction-retraction ####

## Max


#no additional fixed terms/covariates
proret_max0 <- lmer(proret ~ substrate*incline + (1|ind), 
                      data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#one additional fixed term/covariate
proret_max1a <- lmer(proret ~ substrate*incline + Speed.units.s + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

proret_max1b <- lmer(proret ~ substrate*incline + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

proret_max1c <- lmer(proret ~ substrate*incline + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
proret_max2a <- lmer(proret ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

proret_max2b <- lmer(proret ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

proret_max2c <- lmer(proret ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
proret_max3 <- lmer(proret ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                      data = Pb_kine_angles_max_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(proret_max0, proret_max1a, proret_max1b, proret_max1c, proret_max2a, proret_max2b, proret_max2c, proret_max3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

#0, 1a have delta AICc< 2; select 0

anova(proret_max0, proret_max1a)



## Min


#no additional fixed terms/covariates
proret_min0 <- lmer(proret ~ substrate*incline + (1|ind), 
                      data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#one additional fixed term/covariate
proret_min1a <- lmer(proret ~ substrate*incline + Speed.units.s + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

proret_min1b <- lmer(proret ~ substrate*incline + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

proret_min1c <- lmer(proret ~ substrate*incline + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
proret_min2a <- lmer(proret ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

proret_min2b <- lmer(proret ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

proret_min2c <- lmer(proret ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
proret_min3 <- lmer(proret ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                      data = Pb_kine_angles_min_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(proret_min0, proret_min1a, proret_min1b, proret_min1c, proret_min2a, proret_min2b, proret_min2c, proret_min3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 0, 1b have delta AICc< 2; select 0

anova(proret_min0,proret_min1b) 



## Mean


#no additional fixed terms/covariates
proret_mean0 <- lmer(proret ~ substrate*incline + (1|ind), 
                       data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#1 additional fixed term/covariate
proret_mean1a <- lmer(proret ~ substrate*incline + Speed.units.s + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

proret_mean1b <- lmer(proret ~ substrate*incline + tail.body.use + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

proret_mean1c <- lmer(proret ~ substrate*incline + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#2 additional fixed terms/covariates
proret_mean2a <- lmer(proret ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

proret_mean2b <- lmer(proret ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

proret_mean2c <- lmer(proret ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#3 additional fixed terms/covariates
proret_mean3 <- lmer(proret ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_mean_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(proret_mean0, proret_mean1a, proret_mean1b, proret_mean1c, proret_mean2a, proret_mean2b, proret_mean2c, proret_mean3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 0, 1a have delta AICc< 2; select 0

anova(proret_mean0,proret_mean1a) 



##### abduction-adduction ####

## Max


#no additional fixed terms/covariates
abad_max0 <- lmer(abad ~ substrate*incline + (1|ind), 
                      data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#one additional fixed term/covariate
abad_max1a <- lmer(abad ~ substrate*incline + Speed.units.s + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

abad_max1b <- lmer(abad ~ substrate*incline + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

abad_max1c <- lmer(abad ~ substrate*incline + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
abad_max2a <- lmer(abad ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

abad_max2b <- lmer(abad ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

abad_max2c <- lmer(abad ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
abad_max3 <- lmer(abad ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                      data = Pb_kine_angles_max_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(abad_max0, abad_max1a, abad_max1b, abad_max1c, abad_max2a, abad_max2b, abad_max2c, abad_max3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

#2c, 3, 1a, 2c, 2a, 1c, 1b, 0 have delta AICc< 2; select 1a

anova(abad_max0,abad_max1) 
anova(abad_max0,abad_max1b) 
anova(abad_max0,abad_max1c)
anova(abad_max0, abad_max2a)
anova(abad_max0, abad_max2b)
anova(abad_max0, abad_max2c)
anova(abad_max0, abad_max3)

anova(abad_max1a, abad_max2a)
anova(abad_max1a, abad_max2b)
anova(abad_max1a, abad_max2c)
anova(abad_max1a, abad_max3)

anova(abad_max1b, abad_max2a)
anova(abad_max1b, abad_max2b)
anova(abad_max1b, abad_max2c) 
anova(abad_max1b, abad_max3)

anova(abad_max1c, abad_max2a)
anova(abad_max1c, abad_max2b)
anova(abad_max1c, abad_max2c)
anova(abad_max1c, abad_max3)

anova(abad_max2a, abad_max3)
anova(abad_max2b, abad_max3)
anova(abad_max2c, abad_max3)




## Min


#no additional fixed terms/covariates
abad_min0 <- lmer(abad ~ substrate*incline + (1|ind), 
                      data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#one additional fixed term/covariate
abad_min1a <- lmer(abad ~ substrate*incline + Speed.units.s + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

abad_min1b <- lmer(abad ~ substrate*incline + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

abad_min1c <- lmer(abad ~ substrate*incline + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
abad_min2a <- lmer(abad ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

abad_min2b <- lmer(abad ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

abad_min2c <- lmer(abad ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
abad_min3 <- lmer(abad ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                      data = Pb_kine_angles_min_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(abad_min0, abad_min1a, abad_min1b, abad_min1c, abad_min2a, abad_min2b, abad_min2c, abad_min3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 1a, 0 have delta AICc< 2; select 1a

anova(abad_min0,abad_min1a) 



## Mean


#no additional fixed terms/covariates
abad_mean0 <- lmer(abad ~ substrate*incline + (1|ind), 
                       data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#1 additional fixed term/covariate
abad_mean1a <- lmer(abad ~ substrate*incline + Speed.units.s + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

abad_mean1b <- lmer(abad ~ substrate*incline + tail.body.use + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

abad_mean1c <- lmer(abad ~ substrate*incline + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#2 additional fixed terms/covariates
abad_mean2a <- lmer(abad ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

abad_mean2b <- lmer(abad ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

abad_mean2c <- lmer(abad ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#3 additional fixed terms/covariates
abad_mean3 <- lmer(abad ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_mean_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(abad_mean0, abad_mean1a, abad_mean1b, abad_mean1c, abad_mean2a, abad_mean2b, abad_mean2c, abad_mean3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 2c, 1b have delta AICc< 2; select 1b

anova(abad_mean1b,abad_mean2c) 




##### yaw #####

## Max


#no additional fixed terms/covariates
yaw_max0 <- lmer(yaw ~ substrate*incline + (1|ind), 
                      data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#one additional fixed term/covariate
yaw_max1a <- lmer(yaw ~ substrate*incline + Speed.units.s + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

yaw_max1b <- lmer(yaw ~ substrate*incline + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

yaw_max1c <- lmer(yaw ~ substrate*incline + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
yaw_max2a <- lmer(yaw ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

yaw_max2b <- lmer(yaw ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

yaw_max2c <- lmer(yaw ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
yaw_max3 <- lmer(yaw ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                      data = Pb_kine_angles_max_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(yaw_max0, yaw_max1a, yaw_max1b, yaw_max1c, yaw_max2a, yaw_max2b, yaw_max2c, yaw_max3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

#0, 1b have delta AICc< 2; select 0

anova(yaw_max0, yaw_max1b)



## Min


#no additional fixed terms/covariates
yaw_min0 <- lmer(yaw ~ substrate*incline + (1|ind), 
                      data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#one additional fixed term/covariate
yaw_min1a <- lmer(yaw ~ substrate*incline + Speed.units.s + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

yaw_min1b <- lmer(yaw ~ substrate*incline + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

yaw_min1c <- lmer(yaw ~ substrate*incline + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
yaw_min2a <- lmer(yaw ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

yaw_min2b <- lmer(yaw ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

yaw_min2c <- lmer(yaw ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
yaw_min3 <- lmer(yaw ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                      data = Pb_kine_angles_min_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(yaw_min0, yaw_min1a, yaw_min1b, yaw_min1c, yaw_min2a, yaw_min2b, yaw_min2c, yaw_min3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 0, 1a, 1b have delta AICc< 2; select 0

anova(yaw_min0,yaw_min1a) 
anova(yaw_min0,yaw_min1b)



## Mean


#no additional fixed terms/covariates
yaw_mean0 <- lmer(yaw ~ substrate*incline + (1|ind), 
                       data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#1 additional fixed term/covariate
yaw_mean1a <- lmer(yaw ~ substrate*incline + Speed.units.s + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

yaw_mean1b <- lmer(yaw ~ substrate*incline + tail.body.use + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

yaw_mean1c <- lmer(yaw ~ substrate*incline + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#2 additional fixed terms/covariates
yaw_mean2a <- lmer(yaw ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

yaw_mean2b <- lmer(yaw ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

yaw_mean2c <- lmer(yaw ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#3 additional fixed terms/covariates
yaw_mean3 <- lmer(yaw ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_mean_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(yaw_mean0, yaw_mean1a, yaw_mean1b, yaw_mean1c, yaw_mean2a, yaw_mean2b, yaw_mean2c, yaw_mean3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 0, 1a have delta AICc< 2; select 0

anova(yaw_mean0,yaw_mean1a)



#### pitch ####

## Max


#no additional fixed terms/covariates
pitch_max0 <- lmer(pitch ~ substrate*incline + (1|ind), 
                      data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#one additional fixed term/covariate
pitch_max1a <- lmer(pitch ~ substrate*incline + Speed.units.s + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

pitch_max1b <- lmer(pitch ~ substrate*incline + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

pitch_max1c <- lmer(pitch ~ substrate*incline + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
pitch_max2a <- lmer(pitch ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

pitch_max2b <- lmer(pitch ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

pitch_max2c <- lmer(pitch ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_max_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
pitch_max3 <- lmer(pitch ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                      data = Pb_kine_angles_max_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(pitch_max0, pitch_max1a, pitch_max1b, pitch_max1c, pitch_max2a, pitch_max2b, pitch_max2c, pitch_max3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

#0 has delta AICc< 2; select 0



## Min


#no additional fixed terms/covariates
pitch_min0 <- lmer(pitch ~ substrate*incline + (1|ind), 
                      data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#one additional fixed term/covariate
pitch_min1a <- lmer(pitch ~ substrate*incline + Speed.units.s + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

pitch_min1b <- lmer(pitch ~ substrate*incline + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

pitch_min1c <- lmer(pitch ~ substrate*incline + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#two additional fixed terms/covariates
pitch_min2a <- lmer(pitch ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

pitch_min2b <- lmer(pitch ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

pitch_min2c <- lmer(pitch ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_min_spd_categ, REML=FALSE)

#three additional fixed terms/covariates
pitch_min3 <- lmer(pitch ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                      data = Pb_kine_angles_min_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(pitch_min0, pitch_min1a, pitch_min1b, pitch_min1c, pitch_min2a, pitch_min2b, pitch_min2c, pitch_min3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 2c, 1b have delta AICc< 2; select 2c

anova(pitch_min1b, pitch_min2c) 



## Mean


#no additional fixed terms/covariates
pitch_mean0 <- lmer(pitch ~ substrate*incline + (1|ind), 
                       data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#1 additional fixed term/covariate
pitch_mean1a <- lmer(pitch ~ substrate*incline + Speed.units.s + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

pitch_mean1b <- lmer(pitch ~ substrate*incline + tail.body.use + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

pitch_mean1c <- lmer(pitch ~ substrate*incline + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#2 additional fixed terms/covariates
pitch_mean2a <- lmer(pitch ~ substrate*incline + Speed.units.s + tail.body.use + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

pitch_mean2b <- lmer(pitch ~ substrate*incline + Speed.units.s + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

pitch_mean2c <- lmer(pitch ~ substrate*incline + tail.body.use + Premovement + (1|ind), 
                        data = Pb_kine_angles_mean_spd_categ, REML=FALSE)

#3 additional fixed terms/covariates
pitch_mean3 <- lmer(pitch ~ substrate*incline + Speed.units.s + tail.body.use + Premovement + (1|ind), 
                       data = Pb_kine_angles_mean_spd_categ, REML=FALSE)


#b) rank models by AICc score
model.sel(pitch_mean0, pitch_mean1a, pitch_mean1b, pitch_mean1c, pitch_mean2a, pitch_mean2b, pitch_mean2c, pitch_mean3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

# 1b, 2a, 2c have delta AICc< 2; select 1b

anova(pitch_mean1b, pitch_mean2a) 
anova(pitch_mean1b, pitch_mean2c) 




# 8) Calculate point and interval estimates -----------------------------------------------------------------------------------------------------


#### intra-fin angle ####

## Max

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

intrafin_max1b <- lmer(elbow ~  substrate * incline + tail.body.use + (1|ind),
                       Pb_kine_angles_max_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(intrafin_max1b, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(intrafin_max1b))
qqline(residuals(intrafin_max1b), col="blue")
hist(residuals(intrafin_max1b), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(intrafin_max1b, ~ substrate*incline + tail.body.use)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(intrafin_max1b)


## Min

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

intrafin_min0 <- lmer(elbow ~  substrate * incline + (1|ind),
                      data = Pb_kine_angles_min_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(intrafin_min0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(intrafin_min0))
qqline(residuals(intrafin_min0), col="blue")
hist(residuals(intrafin_min0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(intrafin_min0, ~ substrate*incline)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(intrafin_min0)


## Mean

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

intrafin_mean1b <- lmer(elbow ~  substrate * incline + tail.body.use + (1|ind),
                        data = Pb_kine_angles_mean_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(intrafin_mean1b, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(intrafin_mean1b))
qqline(residuals(intrafin_mean1b), col="blue")
hist(residuals(intrafin_mean1b), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(intrafin_mean1b, ~ substrate*incline + tail.body.use)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(intrafin_mean1b)



#### protraction-retraction ####

## Max

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

proret_max0 <- lmer(elbow ~  substrate * incline + (1|ind),
                    Pb_kine_angles_max_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(proret_max0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(proret_max0))
qqline(residuals(proret_max0), col="blue")
hist(residuals(proret_max0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(proret_max0, ~ substrate*incline)


#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(proret_max0)


## Min

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

proret_min0 <- lmer(elbow ~  substrate * incline + (1|ind),
                    data = Pb_kine_angles_min_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(proret_min0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(proret_min0))
qqline(residuals(proret_min0), col="blue")
hist(residuals(proret_min0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(proret_min0, ~ substrate*incline)


#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(proret_min0)


## Mean

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

proret_mean1b <- lmer(elbow ~  substrate * incline + tail.body.use + (1|ind),
                      data = Pb_kine_angles_mean_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(proret_mean1b, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(proret_mean1b))
qqline(residuals(proret_mean1b), col="blue")
hist(residuals(proret_mean1b), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(proret_mean1b, ~ substrate*incline + tail.body.use)


#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(proret_mean1b)



#### abduction-adduction ####

## Max

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

abad_max1a <- lmer(elbow ~  substrate * incline + Speed.units.s.y + (1|ind),
                   Pb_kine_angles_max_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(abad_max1a, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(abad_max1a))
qqline(residuals(abad_max1a), col="blue")
hist(residuals(abad_max1a), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(abad_max1a, ~ substrate*incline + Speed.units.s.y)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(abad_max1a)


## Min

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

abad_min1a <- lmer(elbow ~  substrate * incline + Speed.units.s.y + (1|ind),
                   data = Pb_kine_angles_min_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(abad_min1a, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(abad_min1a))
qqline(residuals(abad_min1a), col="blue")
hist(residuals(abad_min1a), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(abad_min1a, ~ substrate*incline + Speed.units.s.y)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(abad_min1a)


## Mean

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

abad_mean1b <- lmer(elbow ~  substrate * incline + tail.body.use + (1|ind),
                    data = Pb_kine_angles_mean_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(abad_mean1b, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(abad_mean1b))
qqline(residuals(abad_mean1b), col="blue")
hist(residuals(abad_mean1b), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(abad_mean1b, ~ substrate*incline + tail.body.use)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(abad_mean1b)



#### yaw ####

## Max

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

yaw_max0 <- lmer(elbow ~  substrate * incline + (1|ind),
                 Pb_kine_angles_max_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(yaw_max0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(yaw_max0))
qqline(residuals(yaw_max0), col="blue")
hist(residuals(yaw_max0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(yaw_max0, ~ substrate*incline)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(yaw_max0)


## Min

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

yaw_min0 <- lmer(elbow ~  substrate * incline + (1|ind),
                 data = Pb_kine_angles_min_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(yaw_min0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(yaw_min0))
qqline(residuals(yaw_min0), col="blue")
hist(residuals(yaw_min0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(yaw_min0, ~ substrate*incline)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(yaw_min0)


## Mean

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

yaw_mean0 <- lmer(elbow ~  substrate * incline + (1|ind),
                  data = Pb_kine_angles_mean_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(yaw_mean0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(yaw_mean0))
qqline(residuals(yaw_mean0), col="blue")
hist(residuals(yaw_mean0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(yaw_mean0, ~ substrate*incline)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(yaw_mean0)



#### pitch ####

## Max

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

pitch_max0 <- lmer(elbow ~  substrate * incline + (1|ind),
                   Pb_kine_angles_max_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(pitch_max0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(pitch_max0))
qqline(residuals(pitch_max0), col="blue")
hist(residuals(pitch_max0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(pitch_max0, ~ substrate*incline)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(pitch_max0)


## Min

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

pitch_min2c <- lmer(elbow ~  substrate * incline + tail.body.use + Premovement + (1|ind),
                    data = Pb_kine_angles_min_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(pitch_min2c, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(pitch_min2c))
qqline(residuals(pitch_min2c), col="blue")
hist(residuals(pitch_min2c), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(pitch_min2c, ~ substrate*incline + tail.body.use + Premovement)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(pitch_min2c)


## Mean

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

pitch_mean1b <- lmer(elbow ~  substrate * incline + tail.body.use + (1|ind),
                     data = Pb_kine_angles_mean_spd_categ, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(pitch_mean1b, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(pitch_mean1b))
qqline(residuals(pitch_mean1b), col="blue")
hist(residuals(pitch_mean1b), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(pitch_mean1b, ~ substrate*incline + tail.body.use)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(pitch_mean1b)

