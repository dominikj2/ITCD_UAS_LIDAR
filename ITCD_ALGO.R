rm(list=ls())
options(digits = 12)
options(warn=1) 

library(plyr)
library(lidR)
library (pastecs)
library(raster)
library(sp)
library(dplyr)
library(fields)
library(rgl)
library(SDraw)
library(stats)
library(VoxR)
library(stringr)
library(gstat)
library(data.table)
library(fossil)
library(units)
library(EBImage)
library(sf)
library(stars)

# DEFINE DIRECTORIES
Output_Dir <- # DIRECTORY FOR OUTPUTTING THE SEGMENTATION POINT CLOUD
Data_Dir <- # DIRECTORY FOR DATA THAT PROVIDES THE CENTROIDS OF EACH PLOT WITHIN UAS FLIGHTS
            # PLOT CENTROIDS FILE (i.e. PLOT_CENTRE_COORDS.csv) HAS THREE COLUMNS FOR UTM COORDINATES:  c("Flight", "POINT_X", "POINT_Y")
LAS_Dir <- # DIRECTORY FOR LAS FILES REPRESENTING EACH UAS FLIGHT
FUN_Dir <- # DIRECTORY POINTING TO THE FUNCTIONS_ITCD.R FILE WITH ALL THE CUSTOMISED FUNCTIONS USED IN THE ALGORITHM

# SOURCE FUNCTIONS
source(paste(FUN_Dir, "FUNCTIONS_ITCD.R", sep="")) 

# READ IN FLIGHT CETROID DATA
Flight_Centroids <- read.csv(paste(Data_Dir,"PLOT_CENTRE_COORDS.csv", sep=""), row.names = NULL)

# PROJECTION SYSTEM
Proj_Sys <- st_crs(4326) # MELBOURNE WATER CATCHMENTS, VICTORIA, AUSTRALIA
Proj_Sys <- Proj_Sys$input

Flights  <- c(1)
Folder_Version <- "_V1" 

for(f in 1:length(Flights)){

  # GET FLIGHT ID
  Flight_Number <- Flights[f]

  # CREATE OUTPUT FOLDERS
  Version_Output <-paste(Output_Dir, "OUTPUT",Folder_Version, sep="")
  dir.create(file.path(Output_Dir, paste("OUTPUT",Folder_Version, sep="")), showWarnings = FALSE)

  dir.create(file.path(Version_Output, paste("Flight_", Flight_Number, sep="")), showWarnings = FALSE)
  Flight_Output <- paste(Version_Output, "/", "Flight_", Flight_Number, sep="")

  Output_LAS <-paste(Flight_Output, "/LAS" , sep="")
  dir.create(file.path(Flight_Output, "LAS"), showWarnings = FALSE)
  Output_CSV <-paste(Flight_Output, "/CSV" , sep="")
  dir.create(file.path(Flight_Output, "CSV"), showWarnings = FALSE)

  ##########
  # READ LAS
  ##########
  LAS_Flight_Data <- lidR::readLAS(paste(LAS_Dir, "/Flight", Flight_Number, ".laz", sep=""), select = "xyzp")
  crs(LAS_Flight_Data) <- Proj_Sys
  
  ##################################
  # SUBSET LAS FLIGHT TO SELECT PLOT
  ##################################
  SUBSET_Half_Dim <- 40 # NOTE THAT SIZE OF PLOT INCLUDES A 10 M BUFFER THAT IS CLIPPED OUT aT END TO REMOVE EDGE EFFECT

  Flight_Centroid <- Flight_Centroids[which(as.numeric(Flight_Centroids$Flight) == Flight_Number),]
  Subset_Extent <- matrix(c(c((Flight_Centroid$POINT_X - SUBSET_Half_Dim),(Flight_Centroid$POINT_X + SUBSET_Half_Dim)), c((Flight_Centroid$POINT_Y - SUBSET_Half_Dim),(Flight_Centroid$POINT_Y + SUBSET_Half_Dim))), 2,2)
  LAS_subsetF <- clip_rectangle(LAS_Flight_Data, Subset_Extent[1,1], Subset_Extent[1,2], Subset_Extent[2,1], Subset_Extent[2,2])
  crs(LAS_subsetF) <- Proj_Sys 

  # CREATE ATTRIBUTES
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.numeric(0), name="TreeID", desc ="TreeID")
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.numeric(0), name="CentID", desc ="CentID")
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.integer(seq(1, nrow(LAS_subsetF@data), 1)), name="PID", desc ="PID")
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.integer(0), name="VoxID", desc ="VoxID")
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.integer(0), name="VoxCl", desc ="VoxCl")
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.numeric(0), name="Under", desc ="Under")
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.integer(0), name="Strata_ID", desc ="Strata_ID")

  # RELEASE MEMORY
  rm(LAS_Flight_Data)
  gc()

  ###########################################################################################################################################
  ###########################################################################################################################################
  # PARAMETERS
  ###########################################################################################################################################
  ###########################################################################################################################################

  # UNDERSTOREY PARAMETRES
  Para_Grid_Profile_Resolution <- 3
  Para_Under_GridSearchWindow <- 11

  Para_K_BandWidth <- 0.4
  Para_Gap_Threshold_Percent_U <- 0.2
  Para_Gap_Threshold_Percent_U2 <- 0.1
  Para_Gap_Threshold_Percent_U3 <- 0.04

  Para_Vx_R_Res <- 0.2  
  Para_Vx_R_Res2 <- 0.5

  Para_Sl_Z <- Para_Vx_R_Res*2
  Para_EDT_H <- Para_Vx_R_Res*8

  Para_EDT_V <- Para_Vx_R_Res*4
  Para_EDT_V2 <- Para_Vx_R_Res*6

  Para_PCA_Height_Inc <- 6   
  Para_PCA_UniqueZ <- 4  
  Para_PCA_ZRange <- Para_Sl_Z*3 

  Para_Min_TID_Rng <- 3
  Para_Min_DistU1_Rng <- 2
  Para_Half_StZ <- 2

  Para_Min_PolyArea <- 0.04

  Para_Rm_SmRng1 <- 3  
  Para_Rm_MaxZ <- 8 
  Para_Rm_VCnt <- 30   
  Para_Rm_SmRng2 <- 5     

  Para_Rm_DistU1_VoxCnt <- 100

  Para_Threshold_Cls2Top_minZ <- 1/5
  Para_Threshold_Assess_Merger <- 1/3
  Para_Threshold_Assess_Merge2 <- 1/4

  # VIDEO & TREE CHARACTERISATION PARAMETERS.
  Para_Vid_Surrt <- 0.3
  Para_St_Per_Vid <- 5

  Para_Canopy_Base_BW <- 0.5
  Para_Canopy_Base_Thresh <- 0.2

  Para_1st_Consec_Range_Reduce <- 2
  Para_CutOff_Grid_Multiplier <- 2
  Para_CutOff_Grid_Multiplier2 <- 3
  Para_DiamMax <- 2
  Para_Count_Poly <- 1

  Para_Dist_Btw_Ranges <- 2
  Para_Min_TrunkHeight <- 3

  # CREATE EMPTY VECTORS
  Code_Number <- c()
  Count_LAS_TIDs <- c()
  Count_Vox_TIDs <- c()

  TIME_Section <- c()
  TIME_Section_ID <- c()
  TIME__Count_Vox_TID <- c()
  TIME__Count_TID <- c()

  ########################################################################################################################################### 1
  ########################################################################################################################################### 1
  # u LOOP: Removing Understorey
  ########################################################################################################################################### 1
  ########################################################################################################################################### 1

  # CREATE GRID for Computing profiles
  Grid = grid_metrics(LAS_subsetF, length(Z), res= Para_Grid_Profile_Resolution)
  Grid_df <- data.frame(xyFromCell(Grid,seq(1, ncell(Grid),1)))
  Grid_df <- data.frame(Grid_df, Hits = Grid@data@values)
  colnames(Grid_df)[1:2] <- c("X", "Y")

  # EMPTY CSV FOR VOXEL GAP RANGE
  Vox_Output <- data.frame(X=numeric(),
                           Y=numeric(),
                           Strata_ID=numeric(),
                           Start_Largest_Gap=numeric(),
                           End_Largest_Gap=numeric())
  count = 0

  ######################################################################################################################### 1b
  # GRID VEG PROFILE USING DENSITY PLOTS (GAP RANGES FOR EACH GRID)
  #################################################################

  # LOOP THROUGH GRIDS IN y DIRECTION
  Unique_X <- unique(Grid_df$X)
  for(d in 1:length(Unique_X)){
    X_Table <- Grid_df[which(Grid_df$X >= (Unique_X[d] -Para_Grid_Profile_Resolution/2) & Grid_df$X < (Unique_X[d] +Para_Grid_Profile_Resolution/2)),]
    Unique_Y <- unique(X_Table$Y)
    # LOOP THROUGH GRIDS IN Y DIRECTION
    for(e in 1:length(Unique_Y)){
      count = count +1
      Unique_Grid <- LAS_subsetF@data[which(LAS_subsetF$Y >= (Unique_Y[e] - Para_Grid_Profile_Resolution/2) & LAS_subsetF$Y < (Unique_Y[e] + Para_Grid_Profile_Resolution/2) &
                                              LAS_subsetF$X >= (Unique_X[d] - Para_Grid_Profile_Resolution/2) & LAS_subsetF$X < (Unique_X[d] + Para_Grid_Profile_Resolution/2)  ),]
      
      if(nrow(Unique_Grid)> 0){
        # FIRST PROFILE PROCEDURE
        OneGrid_Gap <- GAP_DENSITY_FUN(Unique_Grid$Z,
                                                 Para_BW = Para_K_BandWidth,
                                                 Para_Threshold_Percent = Para_Gap_Threshold_Percent_U,
                                                 Plot = "No",
                                                 TreeID = 1)
  
        Start_Largest_Gap<- OneGrid_Gap$Start_Largest_Gap
        End_Largest_Gap<- OneGrid_Gap$End_Largest_Gap
  
        if(OneGrid_Gap$End_Largest_Gap == -Inf){End_Largest_Gap <- max(Unique_Grid$Z)}
  
        # SECOND PROFILE PROCEDURE
        Reduced_Z <- Unique_Grid$Z[which(Unique_Grid$Z < End_Largest_Gap)]
        OneGrid_Gap <- GAP_DENSITY_FUN(Reduced_Z, Para_BW = Para_K_BandWidth,
                                                 Para_Threshold_Percent = Para_Gap_Threshold_Percent_U2,
                                                 Plot = "No",
                                                 TreeID = 1)
  
        Start_Largest_Gap<- OneGrid_Gap$Start_Largest_Gap
        End_Largest_Gap<- OneGrid_Gap$End_Largest_Gap
  
        # STORE OUTPUT (Start and end of largest gap in veg profile) FOR CSV
        Vox_Data <- c(Unique_X[d], Unique_Y[e], count, Start_Largest_Gap , End_Largest_Gap)
  
        # IF START OR END OF GAP IS INF USING GAP_DENSITY_FUN, THERE IS NO GAP SO GIVE SAME VALUE AS START OR END
        if(length(which(Vox_Data == -Inf)) >0){
          Vox_Data[which(Vox_Data == -Inf)] <- max(Vox_Data[4:5])
        }
        if(length(which(Vox_Data == Inf)) >0){
          Vox_Data[which(Vox_Data == Inf)] <- min(Vox_Data[4:5])
        }
        Vox_Output<- rbind(Vox_Output, Vox_Data)
      }
    }
  } # D LOOP

  # OUTPUT CSV OF GAP RANGES FOR EACH GRID
  colnames(Vox_Output) <- c("X", "Y", "ID", "Top_Under", "Cutoff_Under")

  ######################################################################################################################### 1c
  # GENERATE IDW RASTERS FOR TOP OF UNDERSTOREY
  #############################################
  # CREATE EMPTY RASTER OBJECT AND CONVERT INTO GRID FOR INTERPOLATION
  coordinates(Vox_Output) = ~X+Y
 
  point_extent <- extent(Vox_Output) # 
  point_extent[1] <- point_extent@xmin-5
  point_extent[2] <- point_extent@xmax+5
  point_extent[3] <- point_extent@ymin-5
  point_extent[4] <- point_extent@ymax+5
  
  rast <- raster(ext = point_extent, resolution= Para_Vx_R_Res) # HIGH RESOLUTION
  
  x.range <- as.numeric(c(point_extent@xmin, point_extent@xmax))  
  y.range <- as.numeric(c(point_extent@ymin, point_extent@ymax))
  grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = Para_Vx_R_Res), 
                     y = seq(from = y.range[1], to = y.range[2], by = Para_Vx_R_Res))
  coordinates(grd) <- ~x + y
  gridded(grd) <- TRUE
  
  # INTERPOLATE UNDERSTOREY USING LARGER GRID CELLS UNDERSTOREY TOPS
  idw_Top_Under<-gstat::idw(formula=Top_Under ~ 1, locations=Vox_Output, newdata=grd)
  colnames(idw_Top_Under@data)[1] <- "Under"
  idw_Top_Under_R <- raster(idw_Top_Under)
  
  LAS_subsetF <- merge_spatial(LAS_subsetF, idw_Top_Under_R, attribute = "Under")
 
  ######################################################################################################################## 1d
  # CLEANING UP UNDERSTOREY (VORONOI METHOD) AND CLASSIFYING LIDAR (UNDER(20), STEM (40), CANOPY (30))
  ####################################################################################################

  LAS_subsetF@data$Strata_ID[which(LAS_subsetF$Z <= LAS_subsetF@data$Under)] <- as.integer(20)

  # SEPERATE UNDER AND NOT_UNDER POINTS
  LAS_IDW_Under <- filter_poi(LAS_subsetF, LAS_subsetF@data$Strata_ID==  as.integer(20))
  LAS_IDW_NotUnder <- filter_poi(LAS_subsetF, LAS_subsetF@data$Strata_ID!=  as.integer(20))
  
  # VORONOI METHOD FOR CLEANING UNDERSTOREY
  R_IDW_Under_Max <-   grid_metrics(LAS_IDW_Under, max(Z) , res = 1) 
  Voronoi_Maximas <- VORONOI_FOCAL_FUN(R_IDW_Under_Max,
                                     Win_row=Para_Under_GridSearchWindow,
                                     Win_col=Para_Under_GridSearchWindow, min_or_max = "max")
  Voronoi_Median_Median <- median(Voronoi_Maximas$max_XYZ_MedianVoriPoly$Median_Height_MaximLoc_VoriPoly)

  # LOOP THROUGH VORONOI'S AND GET START/END OF GAP TO RECLASS SOME UNDER INTO "STEM" LAYER
  Local_Maximas_Voronoi_Poly <- Voronoi_Maximas$Maxima_Vori_Poly
  if(length(Local_Maximas_Voronoi_Poly) >= 1){

    for(v in 1:length(Local_Maximas_Voronoi_Poly)){

      # GETTING COORDS OF VORONOI POLYGON AND SUBSETTING THE LIDAR USING VORONOI POLYGON AREA
      One_Voronoi <- Local_Maximas_Voronoi_Poly[v,]
      XY_One_Voronoi <- One_Voronoi@polygons[[1]]@Polygons[[1]]@coords
      LAS_Voronoi <-clip_roi(LAS_IDW_Under, Polygon(data.frame(XY_One_Voronoi))) # NOTE THAT I NEEDED TO PRODUCE NEW CLEAN POLYGON

      # APPLY GAP DENSITY FUNCTION TO EACH VORONOI
      Gap_OneVoronoi <- GAP_DENSITY_FUN(LAS_Voronoi$Z,
                                                  Para_BW = Para_K_BandWidth,
                                                  Para_Threshold_Percent = Para_Gap_Threshold_Percent_U3,
                                                  Plot = "No",
                                                  TreeID = 1)

      # HEIGTH CUT OFF POINT FOR EACH VORONOI POLYGON
      Cutoff_Height <- Gap_OneVoronoi$Start_Largest_Gap
      if(is.infinite(Cutoff_Height)){Cutoff_Height <- Gap_OneVoronoi$End_Largest_Gap}
      if(is.infinite(Cutoff_Height)){Cutoff_Height <- Voronoi_Maximas$max_XYZ_MedianVoriPoly$Start_Largest_Gap[[v]]}

      # ID CHANGES OF Strata_ID FOR POINTS CORRECTED WITH VORONOI METHOD
      LAS_Voronoi_Remove <- filter_poi(LAS_Voronoi, Z > Cutoff_Height)
      if(nrow(LAS_Voronoi_Remove@data) > 0){ 
        LAS_IDW_Under@data$Strata_ID[which(LAS_IDW_Under@data$PID %in% LAS_Voronoi_Remove@data$PID)] <- as.integer(40)
        LAS_IDW_Stem_Add <- filter_poi(LAS_IDW_Under, Strata_ID==  as.integer(40))
        LAS_IDW_Under <- filter_poi(LAS_IDW_Under, Strata_ID==  as.integer(20))
        if(exists("LAS_IDW_Stem_Add_All")){
          LAS_IDW_Stem_Add_All@data <- rbind(LAS_IDW_Stem_Add_All@data, LAS_IDW_Stem_Add@data)
        }else{
          LAS_IDW_Stem_Add_All <- LAS_IDW_Stem_Add
        }
      } # REMOVE POINTS IF ABOVE START LARGEST GAP
    } # LOOP THROUGH VORONOI PLOTS

    # SMALL WORK AROUND ... ADDING THE very tops of Under (i.e. 0.5% of points) back into Overstorey layer
    Quantile_995 <- quantile(LAS_IDW_Under$Z, 0.995)
    LAS_IDW_Under_GtrQ995 <- filter_poi(LAS_IDW_Under, Z > Quantile_995)
    LAS_IDW_NotUnder@data <- rbind(LAS_IDW_NotUnder@data, LAS_IDW_Under_GtrQ995@data)
    LAS_IDW_Under <- filter_poi(LAS_IDW_Under, Z <= Quantile_995)

    # UPDATING STRATA_ID
    LAS_IDW_Under <- filter_poi(LAS_IDW_Under, LAS_IDW_Under@data$Strata_ID==  as.integer(20))
    if(nrow(LAS_IDW_Stem_Add_All@data) > 0){ 
      LAS_IDW_NotUnder@data <- rbind(LAS_IDW_NotUnder@data, LAS_IDW_Stem_Add_All@data)
    }
  } # ONLY BEGIN LOOP IF THERE ARE POLYGONS WITH MAXIMAS OR MINIMAS

  # GENERATE RASTER OF THE UNDERSTOREY
  R_Under <- grid_canopy(LAS_IDW_Under, Para_Vx_R_Res, algorithm = p2r()) 
  R_Under[is.na(R_Under[])] <- 0

  # FINAL FLIGHT SUBSET WITH UNDERSTOREY REMOVED
  LAS_subsetF <- LAS_IDW_NotUnder
  
  ########################################################################################################################################### 2
  ########################################################################################################################################### 2
  # hh LOOP: GENERATES St CENTROID CLUSTERS USING PCA, WS DEN, BRANCHES, AND DISTANCE THRESHOLDS.
  ########################################################################################################################################### 2
  ########################################################################################################################################### 2
  
  ptm <- proc.time()

  # ADD REQUIRED ATTRIBUTES
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.numeric(0), name="CentID", desc ="CentID")
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.numeric(1), name="BranchID", desc ="BranchID") 
  LAS_subsetF <- add_lasattribute(LAS_subsetF, x=as.integer(1), name="BranchSource", desc ="BranchSource") 
  
  # SEQUENCE OF SLICE HEIGHTS
  Z_Inc_Slice <- seq(floor(min(LAS_subsetF$Z)/Para_Sl_Z)*Para_Sl_Z, ceiling(max(LAS_subsetF$Z)/Para_Sl_Z)*Para_Sl_Z, Para_Sl_Z)
  
  # EMPTY CENTROID DF
  Cent_All <- data.frame(TreeID = numeric(),
                         X = numeric(),
                         Y = numeric(),
                         Z = numeric(),
                         CentID = numeric(),
                         BranchID = numeric(),
                         BranchSource = numeric(),
                         VoxCl = numeric())

  # LOOP THROUGH SLICES
  maxCentID <- 0
  for(hh in 1:length(Z_Inc_Slice)){
    
    LAS_oneSl <- filter_poi(LAS_subsetF, Z > (Z_Inc_Slice[hh]) & Z <= (Z_Inc_Slice[hh] + Para_Sl_Z))
    
    # IF THERE IS LIDAR IN THE SLICE
    if(nrow(LAS_oneSl@data) > 0){  
      
      ################################################################################################################################################## 2
      ############################################################################# 2a
      # SLICE DENSITY WATERSHED RASTER AND CENTROID COMPUTATION (WS_CLUSTERING_FUN)
      ############################################################################# 2a
      
      # RASTERISING THE SLICE POINT DENSITY
      R_oneSl_Den <-   grid_metrics(LAS_oneSl, length(Z), res = Para_Vx_R_Res)
      
      if(R_oneSl_Den@ncols > 1 & R_oneSl_Den@nrows > 1){

        # WS CLUSTERING OF SLICE
        Output_WS_Cl <- WS_CLUSTERING_FUN(LAS_oneSl,  Para_Vx_R_Res = Para_Vx_R_Res, LAS_or_Vox = "LAS")
        R_WS <- Output_WS_Cl$R_Stack
        # PROCEED ONLY IF R_WS IS NOT NULL
        if(!is.null(R_WS)){
          LAS_oneSl <- Output_WS_Cl$LAS_WS_Cl
          LAS_oneSl@data$CentID <- LAS_oneSl@data$CentID + maxCentID
          
          WS_XY_maxDen <- Output_WS_Cl$WS_XY_maxDen[,-6]
          WS_XY_maxDen$CentID <- WS_XY_maxDen$CentID + maxCentID 
          maxCentID <- max(WS_XY_maxDen$CentID)
          
          WS_XY_maxDen$TreeID <- 0
          WS_XY_maxDen <- WS_XY_maxDen[, match(c("TreeID", "CentID", "X", "Y", "Z"), colnames(WS_XY_maxDen))]
         
          # #WS_XY_maxDen$BranchID <- 0
          if(nrow(Cent_All) == 0){
            
            WS_XY_maxDen$BranchID <- 1
            WS_XY_maxDen$BranchSource <- 0
            WS_XY_maxDen$TreeID <- seq(1,nrow(WS_XY_maxDen), 1)
            Cent_All <- WS_XY_maxDen
            
            LAS_allSt <- LAS_oneSl
            LAS_oneSl@data$TreeID <- LAS_oneSl@data$CentID
            
          }else{
            
          
            ######################################################################################################################################################## 2
            ############################################################################## 2b
            # DETERMINE WHICH STEM CENTROID CLUSTERS ARE SUITABLE FOR PCA USING CONTRAINTS 
            ############################################################################## 2b

            # CENTROIDS BELOW SLICE hh (ONE SLICE BELOW)
            Cent_OneSl_Belowhh <- Cent_All[which(Cent_All$Z >= Z_Inc_Slice[hh]- Para_Sl_Z & Cent_All$Z < (Z_Inc_Slice[hh])),]
            if(nrow(Cent_OneSl_Belowhh) >= 1){
              
              # CENTROIDS BELOW SLICE hh (Para_PCA_Height_Inc SLICES BELOW) FOR PCA CALCULATION OF STEM/BRANCH POSITIONS
              Cent_nSl <- Cent_All[which(Cent_All$Z >= (Z_Inc_Slice[hh]- Para_Sl_Z*Para_PCA_Height_Inc ) & Cent_All$Z < (Z_Inc_Slice[hh])),]
             
              Summary_Cent_nSl <- as.data.frame(Cent_nSl %>%
                                                  dplyr::select(Z, TreeID) %>%
                                                  dplyr::group_by(TreeID) %>%
                                                  dplyr::summarise(Z_Inc_Slice = length(unique(Z)),
                                                                   Dist_Z = max(Z)-min(Z), .groups = 'drop'))
              
              TID_St_nCl_nRng <- Summary_Cent_nSl$TreeID[which(Summary_Cent_nSl$Z_Inc_Slice >= Para_PCA_UniqueZ  &
                                                                 Summary_Cent_nSl$Dist_Z >=Para_PCA_ZRange )] 

              if(length(TID_St_nCl_nRng) > 0){
                
                # LOOP THROUGH TID
                for(PC in 1:length(TID_St_nCl_nRng)){
                  
                  # # GET ONE TID AND LOOP THROUGH BRANCHES
                  Index_oneSt_PCA <- which(Cent_OneSl_Belowhh$TreeID == TID_St_nCl_nRng[PC])  
                  
                  if(length(Index_oneSt_PCA) > 0){ 
                    
                    # GET ALL CENTROIDS IN TOP FIVE SLICES
                    OneSt_topSl_cent <- as.data.frame(Cent_All[which(Cent_All$TreeID == TID_St_nCl_nRng[PC]),]) # CENTROID_ID
                    OneSt_topSl_cent <- OneSt_topSl_cent[ with(OneSt_topSl_cent, rev(order(Z))),] # ORDER FOR WHILE LOOP

                    ###########################################################
                    # WHILE LOOP EACH BRANCH AND MAKE A NEW PCA FOR EACH BRANCH
                    ###########################################################
                    Unique_BranchID <- rev(sort(unique(OneSt_topSl_cent$BranchID)))
                    
                    #LOOP THROUGH EACH BRANCH OF STEM 
                    for(UB in 1:length(Unique_BranchID)){
                      
                      # FOR EACH TOP BRANCH WORK DOWN STEM TO GET THE TOP CENTROIDS
                      OneBr_OneSt_topSl_cent<- OneSt_topSl_cent[which(OneSt_topSl_cent$BranchID == Unique_BranchID[UB]),] 
                      MinBranchID <- min(OneSt_topSl_cent$BranchID)
                      Count_While <- 0
                      while(length(unique(OneBr_OneSt_topSl_cent$Z)) < Para_PCA_UniqueZ |  
                            (range(OneBr_OneSt_topSl_cent$Z)[2] - range(OneBr_OneSt_topSl_cent$Z)[1]) < Para_PCA_ZRange ){
                        MinBranchSource <- min(OneBr_OneSt_topSl_cent$BranchSource)
                        
                        OneBr_OneSt_topSl_cent <- rbind(OneBr_OneSt_topSl_cent,
                                                        OneSt_topSl_cent[which(OneSt_topSl_cent$BranchID ==  MinBranchSource),])
                        if(min(OneBr_OneSt_topSl_cent$BranchID) == MinBranchID){break()}
                        Count_While <- Count_While + 1
                        print(paste("Count_While................", Count_While, sep=""))
                      }

                      ######################################################################################################################################################## 2
                      #################### 2d
                      # PCA CALC OF BRANCH 
                      #################### 2d
                      
                      # REDUCE BRANCH AND ASSOCIATED PARENT BRANCHES TO TOP FOUR CENTROIDS 
                      OneBr_OneSt_topSl_cent <- OneBr_OneSt_topSl_cent[rev(order(OneBr_OneSt_topSl_cent$Z)),]
                      
                      #RESTRICT CENTROIDS TO A 1.5 m RANGE or SELECTING TOP 4
                      Index_Top_Centroids<- which(OneBr_OneSt_topSl_cent$Z >=(max(OneBr_OneSt_topSl_cent$Z)- Para_PCA_ZRange ))
                      if(length(Index_Top_Centroids) > Para_PCA_UniqueZ ){
                        Index_Top_Centroids <- 1:min(Para_PCA_UniqueZ ,nrow(OneBr_OneSt_topSl_cent)) 
                      }
                      OneBr_OneSt_topSl_cent <- OneBr_OneSt_topSl_cent[Index_Top_Centroids,]
                      
                      
                      if(nrow(OneBr_OneSt_topSl_cent) >= Para_PCA_UniqueZ){
                        ###############################
                        # CALCULATE PCA FOR EACH BRANCH (AND ASSOCIATED PARENT BRANCHES WHEN BRANCH IS TOO SHORT)
                        ###############################
                        Index_Col <- match(c("TreeID", "X", "Y", "Z"),colnames(OneBr_OneSt_topSl_cent))
                        PCA_OneSt <- STEM_PCA_FUN(OneBr_OneSt_topSl_cent[,Index_Col], Para_Slice_Height =  Para_Sl_Z, Up_Down = "UP")

                        ######################################################################################################################################################## 2
                        ################################################## 2e
                        # CHANGING CENTROIDS POSITION WITH PCA CALCULATION
                        ################################################## 2e
                        
                        # UPDATE St POSITION AT hh in Cent_OneSl_Belowhh
                        St_PCA <- PCA_OneSt$STEM_PCA[2,]

                        # DISTANCE BETWEEN PCA LOCATION AND CENTROID LOCATION
                        Index_Top_Cent <- which.min(abs(St_PCA$Z-OneBr_OneSt_topSl_cent$Z))
                        Nearest_to_St_PCA <- OneBr_OneSt_topSl_cent[Index_Top_Cent,]
                        Dist_PCA_to_TopCent<- as.vector(rdist(St_PCA[,match(c("X", "Y"), colnames(St_PCA))],
                                                              Nearest_to_St_PCA[,match(c("X", "Y"), colnames(Nearest_to_St_PCA))]))

                        # UPDATE IF PCA LOCATION IS WITHIN Para_EDT_V M OF TOP CENTOID
                        if(Dist_PCA_to_TopCent < Para_EDT_V ){ 
                          Index_Replace_PCA <- which(Cent_OneSl_Belowhh$CentID == Nearest_to_St_PCA$CentID)
                          Cent_OneSl_Belowhh[Index_Replace_PCA,which(colnames(OneBr_OneSt_topSl_cent) %in% c("X", "Y"))] <- St_PCA[which(colnames(St_PCA) %in% c("X", "Y"))]
                        } # IF PCA CENTROID IS WITHIN 2 M OF BELOW CENTROID
                      } # IF BRANCH HAS ENOUGH CENTROIDS FOR PCA CALCULATION
                    } # LOOP UB (LOOP THROUGH EACH TOP BRANCH)
                  } # IF THE St HAS A TOP CENTROID JUST BELOW H
                } # LOOP THROUGH St (PC)
              }# IF THERE ARE AT LEAST 4 CENTROIDS AT THREE DIFFERENT HEIGHTS
              
              ######################################################################
              # DISTANCE BETWEEN BELOW SLICE AND NEW WS SLICE CLUSTERS
              ######################################################################
              Index_Col_belowhh <- match(c("X", "Y"),colnames(Cent_OneSl_Belowhh)) #c("X", "Y", "Z")
              Index_Col_WS<- match(c("X", "Y"),colnames(Cent_OneSl_Belowhh)) #c("X", "Y", "Z")
              Dist_WS_oneBelowhh <- rdist(Cent_OneSl_Belowhh[,Index_Col_belowhh], WS_XY_maxDen[,Index_Col_WS]) 
              if(all(dim(Dist_WS_oneBelowhh) == c(1,1))){ # (1)
                Dist_WS_oneBelowhh <- matrix(c(as.vector(Dist_WS_oneBelowhh), Para_EDT_V, Para_EDT_V, Para_EDT_V), nrow=2, ncol=2) 
              }
              if(nrow(Dist_WS_oneBelowhh) == 1 & ncol(Dist_WS_oneBelowhh) > 1){ # (2)
                Dist_WS_oneBelowhh <- rbind(Dist_WS_oneBelowhh, rep(Para_EDT_V, ncol(Dist_WS_oneBelowhh))) 
              }
              
              ################################################################################################################################
              # PAIR DISTANCE AND REMOVING PAIRS WHERE TreeID IS NOT THE CLOSEST TO WSMaxDEN (Note: WSMaxDEN can only be close to one TID)
              ################################################################################################################################
              Index_Dist_WS_oneBelowhh <- data.frame(which(apply(Dist_WS_oneBelowhh, 2, function (x) x<Para_EDT_V), arr.ind = TRUE)) #Para_EDT_H
              
              if(length(which(table(Index_Dist_WS_oneBelowhh$col) > 1)) > 0){
                Index_MaxDenWS_moreOne_CloseTID <- as.numeric(names(table(Index_Dist_WS_oneBelowhh$col))[which(table(Index_Dist_WS_oneBelowhh$col) > 1)])
                Index_FINAL_Remove <- c()
                if(length(Index_MaxDenWS_moreOne_CloseTID) > 0){
                  for(RD in 1:length(Index_MaxDenWS_moreOne_CloseTID)){
                    
                    Index_Close_TID <- Index_Dist_WS_oneBelowhh$row[which(Index_Dist_WS_oneBelowhh$col == Index_MaxDenWS_moreOne_CloseTID[RD])]
                    Index_Closest_TID <- Index_Close_TID[which.min(Dist_WS_oneBelowhh[Index_Close_TID, Index_MaxDenWS_moreOne_CloseTID[RD]])]
                    Index_Remove_TID <- setdiff(Index_Close_TID, Index_Closest_TID)
                    Index_FINAL_Remove<- c(Index_FINAL_Remove,
                                           which(Index_Dist_WS_oneBelowhh$row %in% Index_Remove_TID &
                                                   Index_Dist_WS_oneBelowhh$col %in% Index_MaxDenWS_moreOne_CloseTID[RD]))
                  }
                  Index_Dist_WS_oneBelowhh <- Index_Dist_WS_oneBelowhh[-Index_FINAL_Remove,]
                } #(RD LOOP WITHIN)
              }
              
              ################################################################################################################
              # IDENTIFY THE POINTS THAT ARE NOT CLOSEST TO ANY BELOW (within Constraint) AND GIVE THEM BranchID=1 and TreeID = Next
              ################################################################################################################
              if(length(Index_Dist_WS_oneBelowhh$col) > 0){
                WS_XY_maxDen_NewTree <- WS_XY_maxDen[-Index_Dist_WS_oneBelowhh$col,]
              }else{
                WS_XY_maxDen_NewTree <- WS_XY_maxDen
              }

              if(nrow(WS_XY_maxDen_NewTree) > 0){
                WS_XY_maxDen_NewTree$BranchID <- 1
                WS_XY_maxDen_NewTree$BranchSource <- 0
                WS_XY_maxDen_NewTree$TreeID <- WS_XY_maxDen_NewTree$CentID + maxTreeID 
              }
              
              # LOOP THROUGH EACH UNIQUE "BELOW" THAT IS MATCHED WITH "WS"
              unique_Rows_below <- unique(Index_Dist_WS_oneBelowhh$row) 
              #if(hh == 9){browser()}
              WS_XY_maxDen$BranchID <- 0
              WS_XY_maxDen$BranchSource <- 0

              for(UR in 1:length(unique_Rows_below)){
                
                Index_UniqueRow <- Index_Dist_WS_oneBelowhh[which(Index_Dist_WS_oneBelowhh$row == unique_Rows_below[UR]),]
                # UPDATING TreeID
                WS_XY_maxDen$TreeID[Index_UniqueRow$col] <-  Cent_OneSl_Belowhh$TreeID[Index_UniqueRow$row] 
                
                Dist_Nearest_St_to_WSmaxDen <- c() 
                for(i in 1:nrow(Index_UniqueRow)){
                  one_Dist <- Dist_WS_oneBelowhh[Index_UniqueRow[i,1], Index_UniqueRow[i,2]]
                  Dist_Nearest_St_to_WSmaxDen <- c(Dist_Nearest_St_to_WSmaxDen, one_Dist)
                }

                Index_ClosestCol <- Index_UniqueRow[which.min(Dist_Nearest_St_to_WSmaxDen),]
                WS_XY_maxDen$BranchID[Index_ClosestCol$col] <-  Cent_OneSl_Belowhh$BranchID[Index_ClosestCol$row] #GIVEN ALL SAME BranchID BUT CHANGE BELOW JUST CHANGED THIS 
                WS_XY_maxDen$BranchSource[Index_ClosestCol$col] <- Cent_OneSl_Belowhh$BranchID[Index_ClosestCol$row]

                Index_NotClosestCol <- Index_UniqueRow[-which.min(Dist_Nearest_St_to_WSmaxDen),]
                WS_XY_maxDen$BranchSource[Index_NotClosestCol$col] <- Cent_OneSl_Belowhh$BranchID[Index_NotClosestCol$row]
                
                if(nrow(Index_NotClosestCol) > 0){
                  WS_MaxBr <- summary_TIDmaxBr$MaxBr[match(Cent_OneSl_Belowhh$TreeID[Index_NotClosestCol$row], summary_TIDmaxBr$TreeID)]
                  WS_MaxBr_tbl <- table(WS_MaxBr)
                  MoreOneBranch <- which(WS_MaxBr_tbl >  1)
                  
                  if(length(MoreOneBranch)>0){
                    for(MB in 1:length(MoreOneBranch)){
                      oneMoreOne <- as.numeric(names(MoreOneBranch)[MB])
                      Index_oneMoreOne <- which(WS_MaxBr == oneMoreOne)
                      WS_MaxBr[Index_oneMoreOne] <- seq(max(WS_MaxBr[Index_oneMoreOne]), 
                                                        max(WS_MaxBr[Index_oneMoreOne]) + length(max(WS_MaxBr[Index_oneMoreOne])), 1)+1
                    }
                    WS_XY_maxDen$BranchID[Index_NotClosestCol$col]  <- WS_MaxBr[Index_oneMoreOne]
                  }else{
                    WS_XY_maxDen$BranchID[Index_NotClosestCol$col] <-  WS_MaxBr + 1 
                  }
                }
                
              } # LOOP UR

              WS_XY_maxDen <- WS_XY_maxDen[-match(WS_XY_maxDen_NewTree$CentID,WS_XY_maxDen$CentID), ]

              New_Centroid <- rbind(WS_XY_maxDen,WS_XY_maxDen_NewTree)
              
              if(length(which(New_Centroid$TreeID == 0)) > 0){browser()}
              Cent_All <- rbind(Cent_All, New_Centroid)
              
              ############## 2h
              # UPDATING LAS
              ############## 2h
              
              #################################################################
              # MOVE CLUSTER LAS CANOPY TO St LAS (GIVEN LAS APROPRIATE TreeID)
              #################################################################
              flag <- 1

              Index_Add_Close_1 <- which(LAS_oneSl@data$CentID %in% New_Centroid$CentID )
              
              Index_Add_Close_2 <- match(LAS_oneSl@data$CentID, New_Centroid$CentID )
              LAS_oneSl@data$CentID <- New_Centroid$CentID[Index_Add_Close_2]
              LAS_oneSl@data$TreeID <- New_Centroid$TreeID[Index_Add_Close_2]
              
              # MERGING CANOPY TO St
              LAS_allSt@data <- rbind(LAS_allSt@data, LAS_oneSl@data)

            } # IF BELOW SLICE HAS CENTROID
          } # ELSE STATEMENT

          maxTreeID <- max(Cent_All$TreeID)
          
          summary_TIDmaxBr <- as.data.frame(Cent_All %>%
                                              dplyr::group_by(TreeID) %>%
                                              dplyr::summarise(MaxBr = max(BranchID), .groups = 'drop'))
        }
      } 
      print(paste("HH: ", hh, " out of ",length(Z_Inc_Slice), sep="")) 
    } # IF THERE ARE POINTS IN SLICE
  }# HH LOOP
  
  # TIMER 
  Time_Taken <- as.vector(round(((proc.time() - ptm)[3]/60),3))
  TIME_Section <- c(TIME_Section, Time_Taken)
  TIME_Section_ID <- c(TIME_Section_ID, 2)
  TIME__Count_Vox_TID <- c(TIME__Count_Vox_TID, 0)
  TIME__Count_TID <- c(TIME__Count_TID, 0)

  # REMOVE ANY DUPLICATED Cent_All
  Index_Duplicates <- which(duplicated(Cent_All[,match(c("X","Y","Z"), colnames(Cent_All))]))
  if(length(Index_Duplicates) > 0){
    Cent_All <- Cent_All[-Index_Duplicates,]
    }

  # OUTPUTTING FINAL LAS FILES
  LAS_Centroids_All <- lidR::LAS(Cent_All)
  writeLAS(LAS_Centroids_All, paste(Output_LAS,"/F", Flights[f] ,"_Centroids_All.laz",sep=''))

  writeLAS(LAS_allSt, paste(Output_LAS,"/F", Flights[f] ,"_LAS_allSt.laz",sep=''))

  # TABLE OUTPUT TID COUNT
  Code_Number <- c(Code_Number, 2)
  Count_LAS_TIDs <- c(Count_LAS_TIDs , length(unique(LAS_allSt@data$TreeID)))
  Count_Vox_TIDs <- c(Count_Vox_TIDs , 0)
  Output_Table_TID_CodeSummary <- data.frame(Code_Number, Count_LAS_TIDs, Count_Vox_TIDs)
  write.csv(Output_Table_TID_CodeSummary, paste(Output_CSV, "/F", Flight_Number ,"_2_Code_Summary.csv",sep=''), row.names = FALSE)

  write.csv(Cent_All, paste(Output_CSV, "/F", Flight_Number ,"_Cent_All.csv",sep=''), row.names = FALSE)
  
  Cent_All <- as.data.frame(Cent_All)

  ###############
  # BACKUP OUTPUT
  ###############

  Centroids_All_Backup2 <- Cent_All
  LAS_Centroids_All_Backup2 <- LAS_Centroids_All
  LAS_allSt_Backup2 <- LAS_allSt

  Cent_All <- Centroids_All_Backup2
  Cent_All <- as.data.frame(Cent_All)
  LAS_Centroids_All <- LAS_Centroids_All_Backup2
  LAS_allSt <- LAS_allSt_Backup2
  print("Backup2")

  ########################################################################################################## 3
  ########################################################################################################## 3
  # LOOP hhh: MERGING CENTROID CLUSTERS BY PROJECTING PCA BELOW MIN AND ABOV MAX
  ########################################################################################################## 3
  ########################################################################################################## 3
  ptm <- proc.time()

  # GET MIN AND MAX OF EACH Cent_All TID
  SmallestMax_Branch_Cent <- as.data.frame(Cent_All %>%
                                             dplyr::group_by(TreeID, BranchID) %>%
                                             dplyr::summarise(Z_min = min(Z),
                                                              Z_max = max(Z), .groups = 'drop'))
  
  ######################################################################################################################################################## 3
  ############################################################################## 3a
  # LOOPS THROUGH SLICED Cent STARTING WITH THE SMALLEST MAX VALUE OF ALL St
  ############################################################################## 3a
  Class_Cent_All <- class(as.data.frame(Cent_All))
  
  Slice_hhh <- seq(floor(min(SmallestMax_Branch_Cent$Z_max)),ceiling(max(SmallestMax_Branch_Cent$Z_min)) , Para_Sl_Z)
  
  for(hhh in 1:length(Slice_hhh)){
    
    # SUMMARY OF MIN AND MAX FOR EACH BRANCH WITHIN EACH St
    Range_Branch_Cent <- as.data.frame(Cent_All %>%
                                         dplyr::group_by(TreeID, BranchID) %>%
                                         dplyr::summarise(X_Zmin = X[which.min(Z)],
                                                          Y_Zmin = Y[which.min(Z)],
                                                          Z_min = min(Z),
                                                          X_Zmax = X[which.max(Z)],
                                                          Y_Zmax = Y[which.max(Z)],
                                                          Z_max = max(Z),
                                                          Point_Count  = length(Z),
                                                          Min_BrSourceID = min(BranchSource), .groups = 'drop'))
    
    Range_Branch_Cent_StMin <- Range_Branch_Cent[which(Range_Branch_Cent$BranchID ==1 & Range_Branch_Cent$Min_BrSourceID ==0),]
    
    ######################################################################################################################################################## 3
    #################################### 3b
    # PROJECT MIN BRANCH BELOW USING PCA
    #################################### 3b
    
    # GET ALL BRANCHES WITH St ABSOLUTE MIN Z WITHIN Slice_hhh[hhh] + Para_PCA_UniqueZ*Para_Sl_Z
    Orig_Sl_Min_Cent <- Range_Branch_Cent_StMin[which(Range_Branch_Cent_StMin$Z_min > Slice_hhh[hhh] &
                                                        Range_Branch_Cent_StMin$Z_min <= Slice_hhh[hhh] +  Para_PCA_UniqueZ*Para_Sl_Z ),]
    
    # IF THERE ARE NO BRANCH MINs GO TO NEXT hhh LOOP
    if(nrow(Orig_Sl_Min_Cent) > 0){
      
      # CREATE EMPTY DATA.FRAME FOR STORING XYZ OF CLUSTER BOTTOM
      All_PCA_Brch_Min_Below <- data.frame(TreeID=integer(), X=numeric(), Y=numeric(), Z=numeric(),  stringsAsFactors=FALSE) 
      
      # LOOP THOROUGH ALL St/BRANCHES WITH MIN BETWEEN Slice_hhh[hhh] < Z > Slice_hhh[hhh] + Para_PCA_UniqueZ*Para_Sl_Z
      Changed_Min_with_PCA <- 0
      for(SMB in 1:nrow(Orig_Sl_Min_Cent)){ # BELOW LOOP (B)
        
        # GET A MAX OF FOUR CENTROIDS FROM THE 6 BOTTOM SLICES (only for St/Branches TID with at least one centroids in bottom three slices)
        one_PCA_Min_Cent <- Cent_All[which(Cent_All$TreeID %in% Orig_Sl_Min_Cent$TreeID[SMB] &
                                             Cent_All$Z <= Orig_Sl_Min_Cent$Z_min[SMB] +  Para_PCA_Height_Inc *Para_Sl_Z),]  
        one_PCA_Min_Cent <- one_PCA_Min_Cent[which(one_PCA_Min_Cent$Z %in% sort(unique(one_PCA_Min_Cent$Z))[1:min(length(unique(one_PCA_Min_Cent$Z)),Para_PCA_UniqueZ)]),]
        Index_ColNames <- match(c("TreeID", "X", "Y", "Z"), colnames(one_PCA_Min_Cent))
        if(class(one_PCA_Min_Cent)[1] != Class_Cent_All){ one_PCA_Min_Cent <- data.frame(one_PCA_Min_Cent)} # CONVERTING TABLE TO DATAFRAME
        
        # IF THERE ARE AT LEAST Para_PCA_UniqueZ CENTROIDS TO GET PCA BELOW
        if(nrow(one_PCA_Min_Cent) >= Para_PCA_UniqueZ){
          
          # CALCULATE PCA FOR BOTTOM OF TID
          PCA_OneSt <- STEM_PCA_FUN(one_PCA_Min_Cent[,Index_ColNames], Para_Slice_Height =  Para_Sl_Z, Up_Down = "DOWN")  
          St_Pos_Below <- PCA_OneSt$STEM_PCA[1,]
          
          if(!is.na(St_Pos_Below$Z)){
            Dist_PCA_Orig <- rdist(St_Pos_Below[,2:3],
                                   one_PCA_Min_Cent[which.min(one_PCA_Min_Cent$Z),2:3])
            if(Dist_PCA_Orig < Para_EDT_V){
              All_PCA_Brch_Min_Below <- rbind(All_PCA_Brch_Min_Below, St_Pos_Below)
              Changed_Min_with_PCA <- Changed_Min_with_PCA + 1
            }else{ 
              St_Pos_Below <- one_PCA_Min_Cent[which.min(one_PCA_Min_Cent$Z),Index_ColNames]
              colnames(St_Pos_Below) <- c("TreeID", "X", "Y", "Z")
              St_Pos_Below$Z <- round(St_Pos_Below$Z/Para_Sl_Z )*Para_Sl_Z
              All_PCA_Brch_Min_Below <- rbind(All_PCA_Brch_Min_Below, St_Pos_Below)
            }
          }else{
            St_Pos_Below <- one_PCA_Min_Cent[which.min(one_PCA_Min_Cent$Z),Index_ColNames]
            colnames(St_Pos_Below) <- c("TreeID", "X", "Y", "Z")
            St_Pos_Below$Z <- round(St_Pos_Below$Z/Para_Sl_Z )*Para_Sl_Z
            All_PCA_Brch_Min_Below <- rbind(All_PCA_Brch_Min_Below, St_Pos_Below)
          }
        }else{
          St_Pos_Below <- one_PCA_Min_Cent[which.min(one_PCA_Min_Cent$Z),Index_ColNames]
          St_Pos_Below$Z <- round(St_Pos_Below$Z/Para_Sl_Z )*Para_Sl_Z
          All_PCA_Brch_Min_Below <- rbind(All_PCA_Brch_Min_Below, St_Pos_Below)
          }
        }
      
      ######################################################################################################################################################## 3
      ################################## 3c
      # PROJECT MAX St ABOVE USING PCA
      ################################## 3c
      
      Orig_Sl_Max_Cent <- Range_Branch_Cent[which(Range_Branch_Cent$Z_max > Slice_hhh[hhh] &
                                                    Range_Branch_Cent$Z_max <= Slice_hhh[hhh] + Para_Sl_Z),]
      All_PCA_St_Max_Above <- data.frame(TreeID=integer(), X=numeric(), Y=numeric(), Z=numeric(),  stringsAsFactors=FALSE) 
      
      Changed_Max_with_PCA <- 0
      for(SMA in 1:nrow(Orig_Sl_Max_Cent)){ 
        
        one_PCA_Max_Cent <- Cent_All[which(Cent_All$TreeID %in% Orig_Sl_Max_Cent$TreeID[SMA] &
                                             Cent_All$Z >= Orig_Sl_Max_Cent$Z_max[SMA] -  Para_PCA_Height_Inc*Para_Sl_Z  &
                                             Cent_All$Z <= Orig_Sl_Max_Cent$Z_max[SMA]),]  
        one_PCA_Max_Cent <- one_PCA_Max_Cent[which(one_PCA_Max_Cent$Z %in% rev(sort(unique(one_PCA_Max_Cent$Z)))[1:min(length(unique(one_PCA_Max_Cent$Z)),Para_PCA_UniqueZ )]),]
        Index_ColNames <- match(c("TreeID", "X", "Y", "Z"), colnames(one_PCA_Max_Cent))
        if(class(one_PCA_Max_Cent)[1] != Class_Cent_All){one_PCA_Max_Cent <- data.frame(one_PCA_Max_Cent)}
        # IF THERE ARE AT LEAST 4 CENTROIDS TO GET PCA ABOVE
        if(nrow(one_PCA_Max_Cent) >= Para_PCA_UniqueZ  ){
          # CALCULATE PCA_Sl_Max_Cent
          PCA_OneSt <- STEM_PCA_FUN(one_PCA_Max_Cent[,Index_ColNames], Para_Slice_Height =  Para_Sl_Z , Up_Down = "UP") 
          St_Pos_Above <- PCA_OneSt$STEM_PCA[2,]
          if(!is.na(St_Pos_Above$Z)){
            # ONLY USE PCA IF IT IS WITHIN 1 M OF ORIGINAL
            Dist_PCA_Orig <- rdist(St_Pos_Above[,2:3], one_PCA_Max_Cent[which.max(one_PCA_Max_Cent$Z),2:3])
            if(Dist_PCA_Orig < Para_EDT_V){
              All_PCA_St_Max_Above <- rbind(All_PCA_St_Max_Above, St_Pos_Above)
              Changed_Max_with_PCA <- Changed_Max_with_PCA +1
            }else{
              St_Pos_Above <- one_PCA_Max_Cent[which.max(one_PCA_Max_Cent$Z),Index_ColNames]
              colnames(St_Pos_Above) <- c("TreeID", "X", "Y", "Z")
              St_Pos_Above$Z <- round(St_Pos_Above$Z/Para_Sl_Z )*Para_Sl_Z
              All_PCA_St_Max_Above <- rbind(All_PCA_St_Max_Above, St_Pos_Above)
            }
          }else{
            St_Pos_Above <- one_PCA_Max_Cent[which.max(one_PCA_Max_Cent$Z),Index_ColNames]
            colnames(St_Pos_Above) <- c("TreeID", "X", "Y", "Z")
            St_Pos_Above$Z <- round(St_Pos_Above$Z/Para_Sl_Z )*Para_Sl_Z
            All_PCA_St_Max_Above <- rbind(All_PCA_St_Max_Above, St_Pos_Above)
          }
        }else{
          St_Pos_Above <- one_PCA_Max_Cent[which.max(one_PCA_Max_Cent$Z),Index_ColNames]
          colnames(St_Pos_Above) <- c("TreeID", "X", "Y", "Z")
          St_Pos_Above$Z <- round(St_Pos_Above$Z/Para_Sl_Z )*Para_Sl_Z
          All_PCA_St_Max_Above <- rbind(All_PCA_St_Max_Above, St_Pos_Above)
          }
        }
      
      # MERGE (1) MAX STEN CENTROIDS WITH (2) St CENTROIDS WITHIN SLICE THAT ARE NOT MAX CENTROIDS
      All_Cent_hhh <- rbind(All_PCA_St_Max_Above, Cent_All[which(Cent_All$Z > Slice_hhh[hhh] &
                                                                   Cent_All$Z <= Slice_hhh[hhh] + Para_Sl_Z &
                                                                   !(Cent_All$TreeID %in% All_PCA_St_Max_Above$TreeID )),Index_ColNames])
      if(nrow(All_Cent_hhh) > 0){
        
        ########################################################################################################## 3d
        # DISTANCE ANALYSIS USING BOTTOM OF BRANCH AND ALL CENTROIDS IN SLICE (PCA ADJUSTED TOP FOR SOME CENTROID)
        ########################################################################################################## 3d

        Index_Col1 <- match(c("X", "Y", "Z"), colnames(All_PCA_Brch_Min_Below))
        Index_Col2 <- match(c("X", "Y", "Z"), colnames(All_Cent_hhh))
        Dist_PCA_Branch_Cl <- as.data.frame(rdist(as.data.frame(All_PCA_Brch_Min_Below)[,Index_Col1],
                                                  as.data.frame(All_Cent_hhh)[,Index_Col2]))
        
        if(all(dim(Dist_PCA_Branch_Cl) == c(1,1))){ 
          Dist_PCA_Branch_Cl <- matrix(c(as.vector(Dist_PCA_Branch_Cl), Para_EDT_V, Para_EDT_V, Para_EDT_V), nrow=2, ncol=2)
        }
        if(nrow(Dist_PCA_Branch_Cl) == 1 & ncol(Dist_PCA_Branch_Cl) > 1){ # (2)
          Dist_PCA_Branch_Cl <- rbind(Dist_PCA_Branch_Cl, rep(Para_EDT_V, ncol(Dist_PCA_Branch_Cl))) # Para_EDT_V2
        }
        
        #############################################################################
        # INDEXING CLOSE PAIRS AND CLEAN (Remove Same Tree and Paired two Orig Trees)
        #############################################################################
        
        #GET ALL MATCHES WITHIN X M OF EACH St CENTROID
        Index_Dist_PCA_Branch_Cl <- data.frame(which(apply(Dist_PCA_Branch_Cl, 2, function (x) x<Para_EDT_V & x>0), arr.ind = TRUE)) 
        
        # EMPTY CLOSE PAIRS
        Pairs_Centroids <- data.frame(Close_BotBranch_TID  = numeric(),
                                      Close_hhhCent_TopSt_TID  = numeric())
        
        if(nrow(Index_Dist_PCA_Branch_Cl) > 0){
          Close_BotBranch_TID <- All_PCA_Brch_Min_Below$TreeID[Index_Dist_PCA_Branch_Cl$row]
          Close_hhhCent_TopSt_TID <- All_Cent_hhh$TreeID[Index_Dist_PCA_Branch_Cl$col]
          
          # Remove Close TID pairs from same St (TreeID)
          Pairs_Centroids <- data.frame(Close_BotBranch_TID, Close_hhhCent_TopSt_TID)
          Same_Tree <- which(Pairs_Centroids$Close_BotBranch_TID == Pairs_Centroids$Close_hhhCent_TopSt_TID)
          
          if(length(Same_Tree) > 0){
            Pairs_Centroids <- Pairs_Centroids[-Same_Tree,]
            Index_Dist_PCA_Branch_Cl <- Index_Dist_PCA_Branch_Cl[-Same_Tree,]
          }
        } # IF INITIAL INDEX HAS PAIRS INDEX
        
        Pairs_Centroids <- data.frame(Pairs_Centroids)
        colnames(Pairs_Centroids) <- c("TID_MinBrch_Below", "TID_MaxSt_Above")
        
        # IF THERE ARE PAIRS
        if(nrow(Pairs_Centroids) > 0){
          
          ################################################### 3e
          # UPDATING EVERYTHING USING FINAL PAIRS FOR MERGING
          ################################################### 3e
          
          # UPDATE TreeID FOR LAS FILE AND Cent_All
          Index_1 <- which(LAS_allSt@data$TreeID %in%  Pairs_Centroids$TID_MinBrch_Below )
          Index_2 <- match(LAS_allSt@data$TreeID[Index_1],  Pairs_Centroids$TID_MinBrch_Below )
          LAS_allSt@data$TreeID[Index_1] <- as.integer(Pairs_Centroids$TID_MaxSt_Above[Index_2])
          
          Index_1 <- which(Cent_All$TreeID %in%  Pairs_Centroids$TID_MinBrch_Below )
          Index_2 <- match(Cent_All$TreeID[Index_1],  Pairs_Centroids$TID_MinBrch_Below )
          Cent_All$TreeID[Index_1] <- as.integer(Pairs_Centroids$TID_MaxSt_Above[Index_2])
          
          # UPDATE BRANCHID, BRANCHSOURCEID for Cent_All
          Index_St_MergeToo <- which(Cent_All$TreeID %in%  Pairs_Centroids$TID_MaxSt_Above)
          Branch_ID_Max <- max(Cent_All$BranchID[Index_St_MergeToo])
          Cent_All$BranchID[Index_1] <- Cent_All$BranchID[Index_1] + Branch_ID_Max
        } # UPDATE IF CLOSE PAIRS
      }
    } # IF THERE ARE BRANCHES WITH MIN HEIGHT WITHIN THE hhh LOOP + 3*hhh
    print(paste(hhh, " of ", length(Slice_hhh)))
  } #  LOOP SLICE  hhh

  LAS_Centroids_All <- LAS(Cent_All)
  
  Time_Taken <- as.vector(round(((proc.time() - ptm)[3]/60),3))
  TIME_Section <- c(TIME_Section, Time_Taken)
  TIME_Section_ID <- c(TIME_Section_ID, 3)
  TIME__Count_Vox_TID <- c(TIME__Count_Vox_TID, 0)
  TIME__Count_TID <- c(TIME__Count_TID, 0)
  
  Output_Table_TIME_Section_ID <- data.frame(TIME_Section_ID, TIME_Section, TIME__Count_Vox_TID, TIME__Count_TID)
  write.csv(Output_Table_TIME_Section_ID, paste(Output_CSV, "/F", Flight_Number ,"_3_TIME_TAKEN.csv",sep=''), row.names = FALSE)
  
  # TABLE OUTPUT TID COUNT
  Code_Number <- c(Code_Number, 3)
  Count_LAS_TIDs <- c(Count_LAS_TIDs , length(unique(LAS_allSt@data$TreeID)))
  Count_Vox_TIDs <- c(Count_Vox_TIDs , 0)
  Output_Table_TID_CodeSummary <- data.frame(Code_Number, Count_LAS_TIDs, Count_Vox_TIDs)
  write.csv(Output_Table_TID_CodeSummary, paste(Output_CSV, "/F", Flight_Number ,"_3_Code_Summary.csv",sep=''), row.names = FALSE)
  
  ###############
  # BACKUP OUTPUT
  ###############
  
  Centroids_All_Backup3 <- Cent_All
  LAS_allSt_Backup3 <- LAS_allSt
  
  Cent_All <- Centroids_All_Backup3
  Cent_All <- as.data.frame(Cent_All)
  LAS_allSt <- LAS_allSt_Backup3
  print("Backup3")
  
  writeLAS(LAS_allSt, paste(Output_LAS,"/F", Flights[f] ,"_LAS_allSt_Backup3.laz",sep=''))
  write.csv(Cent_All, paste(Output_CSV, "/F", Flight_Number ,"_Cent_All_Backup3.csv",sep=''), row.names = FALSE)

  ########################################################################################################### 4
  ########################################################################################################### 4
  # V LOOP: VOXELISE AND CLUSTER UNIQUE TID (obj.rec) ... TAKES LONG TIME!!
  ########################################################################################################### 4
  ########################################################################################################### 4
  ptm <- proc.time()

  ### POINTS WITH NA TreeID ARE GIVEN VALUE Zero
  LAS_allSt@data$TreeID[is.na(LAS_allSt@data$TreeID)] <- 0
  
  LAS_Zero_NA_WS <-  filter_poi(LAS_allSt, TreeID == 0)
  writeLAS(LAS_Zero_NA_WS, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Zero_NA_WS.laz",sep=''))
  LAS_allSt <-  filter_poi(LAS_allSt, TreeID != 0)

  LAS_allSt_Orig <- LAS_allSt

  Cent_All <- data.frame(Cent_All[,match(c("TreeID", "X", "Y", "Z", "CentID", "BranchID", "BranchSource"), colnames(Cent_All))])
  Cent_All$VoxCl <- 0 
  
  LAS_Centroids_All <- LAS(Cent_All)
  Cent_All <- as.data.frame(Cent_All)
  Centroids_All_Orig <-   Cent_All 

  Unique_TID <- na.omit(unique(LAS_allSt@data$TreeID)) 
  
  max_Vox_Cl <- 0
  Max_Vox_ID <- 0
  
  if(exists("Cent_All_VoxCl")){rm(Cent_All_VoxCl)}
  
  # LOOPS THROUGH EACH TreeID 
  for (v in 1:length(Unique_TID)){

    ########################### 4a
    # VOXELISE LAS FOR ONE TID
    ########################### 4a
    LAS_oneTID <- filter_poi(LAS_allSt, TreeID ==  Unique_TID[v])
    Vx_oneTID_Output <- VOXEL_FUN(LAS_oneTID, Para_Vox_Res = Para_Vx_R_Res)
    LAS_oneTID <-Vx_oneTID_Output[[1]]
    Vx_oneTID <- Vx_oneTID_Output[[2]]
    
    # CLUSTERS Vx_oneTID 
    if(nrow(Vx_oneTID) > 1){ 
      Index_Col1 <- match(c("X", "Y", "Z"), colnames(Vx_oneTID))
      Vx_oneTID_Cl <- VoxR::distance_clustering(data = Vx_oneTID[,..Index_Col1], 
                                                d_clust=(Para_Vx_R_Res*2 + 0.01)) # REPLACE  (VoxR 1.0 from 0.5.1) During lidR 3.0.5
      colnames(Vx_oneTID_Cl)[4] <- "VoxCl"
    }else{
      Index_Col1 <- match(c("X", "Y", "Z"), colnames(Vx_oneTID))
      Vx_oneTID_Cl <- cbind(Vx_oneTID[,..Index_Col1], VoxCl = 1)
    }
    
    # UPDATE CLUSTER ID SO IT INCREMENTALLY INCREASES
    ID_Orig_VoxCl <- as.numeric(names(table(Vx_oneTID_Cl$VoxCl))) 
    if(length(ID_Orig_VoxCl) == 1){
      if(ID_Orig_VoxCl == 0){ID_Final_VoxCl <- 0}
      if(ID_Orig_VoxCl == 1){ID_Final_VoxCl <- max_Vox_Cl+1} 
    } 

    if(length(which(ID_Orig_VoxCl == 0)) == 1 & length(ID_Orig_VoxCl) > 1){
      ID_Final_VoxCl <- c(0, max_Vox_Cl + seq(1, (length(ID_Orig_VoxCl)-1),1))}
    if(length(which(ID_Orig_VoxCl == 0)) == 0 & length(ID_Orig_VoxCl) > 1){
      ID_Final_VoxCl <-  c(max_Vox_Cl + seq(1, (length(ID_Orig_VoxCl)),1))}
    
    Vx_oneTID_Cl$VoxCl <- ID_Final_VoxCl[match(Vx_oneTID_Cl$VoxCl, ID_Orig_VoxCl)]
    max_Vox_Cl <- max(Vx_oneTID_Cl$VoxCl)
    
    # GENERATE LAS USING VOXELS WITH VOXEL CLUSTERS
    Vx_oneTID_Cl <- data.frame(TreeID = Unique_TID[v],
                               Vx_oneTID_Cl,
                               Vox_PntCnt = Vx_oneTID$Point_Count_In_Vox,
                               VoxID = Vx_oneTID$VoxID +Max_Vox_ID)
    
    # UPDATING LAS_oneTID POINTS WITH Vox_ID and VoxCl
    LAS_oneTID@data$VoxID <- LAS_oneTID@data$VoxID + Max_Vox_ID
    Index_Vox_Cl<- match(LAS_oneTID@data$VoxID, Vx_oneTID_Cl$VoxID)
    LAS_oneTID@data$VoxCl <- Vx_oneTID_Cl$VoxCl[Index_Vox_Cl]
    
    # CREATING LAS_Vx_oneTID 
    LAS_Vx_oneTID <- LAS(Vx_oneTID_Cl)
    LAS_Vx_oneTID <- add_lasattribute(LAS_Vx_oneTID, x=as.numeric(Vx_oneTID_Cl$TreeID), name="TreeID", desc ="TreeID")
    LAS_Vx_oneTID <- add_lasattribute(LAS_Vx_oneTID, x=as.numeric(0), name="Under", desc ="Under")
    LAS_Vx_oneTID <- add_lasattribute(LAS_Vx_oneTID, x=as.integer(Vx_oneTID_Cl$VoxID), name="VoxID", desc ="VoxID")
    LAS_Vx_oneTID <- add_lasattribute(LAS_Vx_oneTID, x=as.integer(Vx_oneTID_Cl$Vox_PntCnt), name="Vox_PntCnt", desc ="Vox_PntCnt")
    LAS_Vx_oneTID <- add_lasattribute(LAS_Vx_oneTID, x=as.integer(Vx_oneTID_Cl$VoxCl), name="VoxCl", desc ="VoxCl")
    
    Max_Vox_ID <- max(Vx_oneTID_Cl$VoxID)
    
    # GENERATE LAS WITH ALL TreeID CLUSTERED
    if(v == 1){
      LAS_Vox_TID <- LAS_Vx_oneTID
      LAS_TID <-LAS_oneTID
    }else{
      LAS_Vox_TID@data <- rbind(LAS_Vox_TID@data, LAS_Vx_oneTID@data)
      LAS_TID@data <- rbind(LAS_TID@data, LAS_oneTID@data)
    }
    
    ################################## 4b
    # UPDATING CENTROIDA WITH VoxCl ID
    ################################## 4b
    
    One_TID_Cent <- Centroids_All_Orig[which(Centroids_All_Orig$TreeID == Unique_TID[v]),]
    if(nrow(One_TID_Cent) > 0){
      LAS_Vx_oneTID_ForCent <- filter_poi(LAS_Vox_TID, TreeID == Unique_TID[v] & VoxCl > 0)
      if(nrow(LAS_Vx_oneTID_ForCent@data) > 0){ 
        
        Index_Col1 <- match(c("X", "Y", "Z"), colnames(LAS_Vx_oneTID_ForCent@data))
        Index_Col2 <-match(c("X", "Y", "Z"), colnames(One_TID_Cent))
        Dist_Vox_Cent <- rdist(LAS_Vx_oneTID_ForCent@data[,..Index_Col1],
                               as.data.frame(One_TID_Cent)[,Index_Col2])
        Index_Dist_Vox_Cent<- apply(Dist_Vox_Cent,2,which.min)
        Cent_Vox_ID <- LAS_Vx_oneTID_ForCent@data$VoxCl[Index_Dist_Vox_Cent]
        One_TID_Cent$VoxCl <- Cent_Vox_ID
      }
    }
    
    if(!exists("Cent_All_VoxCl")){
      Cent_All_VoxCl <- One_TID_Cent
    }else{
      Cent_All_VoxCl <- rbind(Cent_All_VoxCl, One_TID_Cent)
    }
    print(paste(v, "_",  dim(Cent_All_VoxCl)[1], "   out of  ", length(Unique_TID)))
  } # V LOOP
  
  
  # TIMER   # ptm <- proc.time()
  Time_Taken <- as.vector(round(((proc.time() - ptm)[3]/60),3))
  TIME_Section <- c(TIME_Section, Time_Taken)
  TIME_Section_ID <- c(TIME_Section_ID, 4)
  TIME__Count_Vox_TID <- c(TIME__Count_Vox_TID, length(unique(LAS_Vox_TID@data$TreeID)))
  TIME__Count_TID <- c(TIME__Count_TID, length(unique(LAS_TID@data$TreeID)))
  Output_Table_TIME_Section_ID <- data.frame(TIME_Section_ID, TIME_Section, TIME__Count_Vox_TID, TIME__Count_TID)
  
  write.csv(Output_Table_TIME_Section_ID, paste(Output_CSV, "/F", Flight_Number ,"_4_TIME_TAKEN.csv",sep=''), row.names = FALSE)
  
  # TABLE OUTPUT TID COUNT
  Code_Number <- c(Code_Number, 4)
  Count_LAS_TIDs <- c(Count_LAS_TIDs , length(unique(LAS_TID@data$TreeID)))
  Count_Vox_TIDs <- c(Count_Vox_TIDs , length(unique(LAS_Vox_TID@data$TreeID)))
  Output_Table_TID_CodeSummary <- data.frame(Code_Number, Count_LAS_TIDs, Count_Vox_TIDs)
  write.csv(Output_Table_TID_CodeSummary, paste(Output_CSV, "/F", Flight_Number ,"_4_Code_Summary.csv",sep=''), row.names = FALSE)
  
  #####################
  # BACKUP VOXELISATION
  #####################
  LAS_TID_Backup4 <- LAS_TID
  LAS_Vox_TID_Backup4 <- LAS_Vox_TID
  Centroids_All_Backup4 <- Cent_All
  
  LAS_TID <- LAS_TID_Backup4
  LAS_Vox_TID <- LAS_Vox_TID_Backup4
  Cent_All <- Centroids_All_Backup4
  print("Backup4")

  ####################################################################################################################################### 5
  ####################################################################################################################################### 5
  # (CLEANING TID)  SEPERATE TID ZERO
  ####################################################################################################################################### 5
  ####################################################################################################################################### 5

  ptm <- proc.time()

  ############################################################################### 5a
  # SEPERATE TID>0 VOXELS FROM TID==0 VOXELS AND MOVE SMALL TID INTO ZERO TID
  ############################################################################### 5a
  
  # GET VOX and LAS WITH NO TreeID
  LAS_Zero <-  filter_poi(LAS_TID, TreeID == 0)
  LAS_Vox_Zero <-  filter_poi(LAS_Vox_TID, TreeID == 0)
  
  LAS_TID <-  filter_poi(LAS_TID, TreeID > 0)
  LAS_Vox_TID <-  filter_poi(LAS_Vox_TID, TreeID > 0)

  # MOVING SMALL TID THAT IS LOW TO GROUND (OR HAS FEW VOXELS AND SLIGHTLY HEIGHER FROM GROUND) TO ZERO
  Summary_TID_NotGnd <- as.data.frame(LAS_Vox_TID@data %>%
                                        dplyr::group_by(TreeID) %>%
                                        dplyr::summarise(Vox_Count = length(Z),
                                                         MinZ = min(Z),
                                                         MaxZ = max(Z),
                                                         Range_Z = range(Z)[2]- range(Z)[1], .groups = 'drop'))
  
  TID_Remove <- Summary_TID_NotGnd$TreeID[which((Summary_TID_NotGnd$Range_Z < Para_Rm_SmRng1 &
                                                   Summary_TID_NotGnd$MaxZ < Para_Rm_MaxZ)|
                                                  (Summary_TID_NotGnd$Vox_Count < Para_Rm_VCnt &
                                                     Summary_TID_NotGnd$Range_Z < Para_Rm_SmRng2))]
  
  # IDENTIFY SMALL CLUSTERS THAT GO INTO ZERO
  LAS_Vox_TID_Remove <- filter_poi(LAS_Vox_TID, (TreeID %in% TID_Remove))
  LAS_Vox_TID_Remove@data$TreeID <- as.integer(0)
  LAS_TID_Remove <- filter_poi(LAS_TID, (TreeID %in% TID_Remove))
  LAS_TID_Remove@data$TreeID <- as.integer(0)
  
  # UPDATE LAS AND VOX Zero TID and NonZero TID
  LAS_Vox_Zero@data <- rbind(LAS_Vox_Zero@data,
                             LAS_Vox_TID_Remove@data)
  LAS_Zero@data <- rbind(LAS_Zero@data,
                         LAS_TID_Remove@data)
  
  LAS_Vox_TID <- filter_poi(LAS_Vox_TID, !(TreeID %in% TID_Remove))
  LAS_TID <- filter_poi(LAS_TID, !(TreeID %in% TID_Remove))
  
  # UPDATE CENTROIDS BY REMOVING THE SMALL CLUSTERED ONES
  Centroids_Zero <-Cent_All[which(Cent_All$TreeID %in% TID_Remove),]
  Cent_All <- Cent_All[-which(Cent_All$TreeID %in% TID_Remove),]
  
  # ANOTHER BACKUP AS ABOVE TAKES LONG
  LAS_Zero_Backup5 <- LAS_Zero
  LAS_Vox_Zero_Backup5 <- LAS_Vox_Zero
  LAS_TID_Backup5 <- LAS_TID
  LAS_Vox_TID_Backup5 <- LAS_Vox_TID
  
  LAS_Zero <- LAS_Zero_Backup5
  LAS_Vox_Zero <- LAS_Vox_Zero_Backup5
  LAS_TID <- LAS_TID_Backup5
  LAS_Vox_TID <- LAS_Vox_TID_Backup5
  print("Backup5")
  
  # TIMER   # ptm <- proc.time()
  Time_Taken <- as.vector(round(((proc.time() - ptm)[3]/60),3))
  TIME_Section <- c(TIME_Section, Time_Taken)
  TIME_Section_ID <- c(TIME_Section_ID, 5)
  TIME__Count_Vox_TID <- c(TIME__Count_Vox_TID, length(unique(LAS_Vox_TID@data$TreeID)))
  TIME__Count_TID <- c(TIME__Count_TID, length(unique(LAS_TID@data$TreeID)))
  Output_Table_TIME_Section_ID <- data.frame(TIME_Section_ID, TIME_Section, TIME__Count_Vox_TID, TIME__Count_TID)
  write.csv(Output_Table_TIME_Section_ID, paste(Output_CSV, "/F", Flight_Number ,"_5_TIME_TAKEN.csv",sep=''), row.names = FALSE)

  # TABLE OUTPUT TID COUNT
  Code_Number <- c(Code_Number, 5)
  Count_LAS_TIDs <- c(Count_LAS_TIDs , length(unique(LAS_TID@data$TreeID)))
  Count_Vox_TIDs <- c(Count_Vox_TIDs , length(unique(LAS_Vox_TID@data$TreeID)))
  Output_Table_TID_CodeSummary <- data.frame(Code_Number, Count_LAS_TIDs, Count_Vox_TIDs)
  write.csv(Output_Table_TID_CodeSummary, paste(Output_CSV, "/F", Flight_Number ,"_5_Code_Summary.csv",sep=''), row.names = FALSE)
  
  LAS_Centroids_All <- lidR::LAS(Cent_All)
  LAS_Centroids_All@data$TreeID <- as.integer(LAS_Centroids_All@data$TreeID)
  LAS_Centroids_All <- add_lasattribute(LAS_Centroids_All, x=LAS_Centroids_All@data$TreeID, name="BranchID", desc ="BranchID") 
  writeLAS(LAS_Centroids_All, paste(Output_LAS,"/F", Flights[f] ,"_Centroids_Backup5.laz",sep=''))
  
  LAS_Zero@data$TreeID <- as.integer(LAS_Zero@data$TreeID)
  writeLAS(LAS_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Zero_Backup5.laz",sep=''))
  LAS_Vox_Zero@data$TreeID <- as.integer(LAS_Vox_Zero@data$TreeID)
  writeLAS(LAS_Vox_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_Zero_Backup5.laz",sep=''))
  
  LAS_TID@data$TreeID <- as.integer(LAS_TID@data$TreeID)
  writeLAS(LAS_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_TID_Backup5.laz",sep=''))
  LAS_Vox_TID@data$TreeID <- as.integer(LAS_Vox_TID@data$TreeID)
  writeLAS(LAS_Vox_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_TID_Backup5.laz",sep=''))
  
  LAS_IDW_Under@data$TreeID <- as.integer(0)
  writeLAS(LAS_IDW_Under, paste(Output_LAS,"/F", Flight_Number ,"_LAS_IDW_Under_Backup5.laz",sep=''))

  ####################################################################################################################################### 7
  ####################################################################################################################################### 7
  # LOOP VC: BOUNDING BOX OF VoxCl AND SURROUND STEMS (Obj.rec) TO CLEAN UP EACH St (BY MOVING CLUSTERS) BUT DOESN'T GENERATE ANY NEW St
  ####################################################################################################################################### 7
  ####################################################################################################################################### 7

  ptm <- proc.time()
  
  # GET LOWEST HEIGHT OF EACH TID
  MinZ_all_TID <- as.data.frame(LAS_Vox_TID@data %>%
                                  dplyr::group_by(TreeID) %>%
                                  dplyr::summarise(MinZ = min(Z),
                                                   .groups = 'drop'))
  
  # START WITH THE TID THAT HAS LOWEST MINZ
  Unique_TID <- MinZ_all_TID$TreeID[order(MinZ_all_TID$MinZ)]
  HalfHeight_LAS <- max(LAS_Vox_TID$Z)/Para_Half_StZ
  
  Max_TID <- max(Unique_TID)
  Max_Vox_Cl <- max(LAS_Vox_TID@data$VoxCl)
  
  Empty_Vx <- c()
  LAS_Empty <- filter_poi(LAS_Vox_TID, TreeID < -10000) 
  
  ###############################
  # LOOP EACH TID, CLUSTER VOXELS
  ############################### 7a
  
  for(VC in 1:length(Unique_TID)){
    if(nrow(LAS_Vox_TID@data) == 0){browser()}
    
    LAS_Update_Changes <- LAS_Empty
    
    print(paste("VoxCnt:",nrow(LAS_Vox_TID@data),"..TID:", Unique_TID[VC],"........................................................Doing:",VC, "....Total:", length(Unique_TID), sep=""))
    
    # ONE St VOX CLUSTERS 
    LAS_oneTID_Vx <- filter_poi(LAS_Vox_TID, TreeID == Unique_TID[VC])
    
    # CHECKING 
    if(nrow(LAS_oneTID_Vx@data) == 0){
      Empty_Vx <- c(Empty_Vx ,Unique_TID[VC] )
    }
    
    # IF VOX ID STILL EXISTS
    if(nrow(LAS_oneTID_Vx@data) > 0){
      
      # UPDATE BOUNDARY OF ALL TID 
      Bound_all_TID <- as.data.frame(LAS_Vox_TID@data %>%
                                       dplyr::group_by(TreeID) %>%
                                       dplyr::summarise(Vox_Count = length(Z),
                                                        MinX = min(X)-Para_Vx_R_Res,
                                                        MaxX = max(X)+Para_Vx_R_Res,
                                                        MinY = min(Y)-Para_Vx_R_Res,
                                                        MaxY = max(Y)+Para_Vx_R_Res,
                                                        MinZ = min(Z)-Para_Vx_R_Res,
                                                        MaxZ = max(Z)+Para_Vx_R_Res, 
                                                        .groups = 'drop'))
      
      Index_oneTID <- which(Bound_all_TID$TreeID == Unique_TID[VC])
      
      # GET Zero_Vox IN BOUNDING BOX OF oneTID
      LAS_oneTID_VoxZero <- filter_poi(LAS_Vox_Zero, X >= Bound_all_TID$MinX[Index_oneTID] &
                                        X <= Bound_all_TID$MaxX[Index_oneTID] &
                                        Y >= Bound_all_TID$MinY[Index_oneTID] &
                                        Y <= Bound_all_TID$MaxY[Index_oneTID] &
                                        Z >= Bound_all_TID$MinZ[Index_oneTID] &
                                        Z <= Bound_all_TID$MaxZ[Index_oneTID])
      
      #  MOVE BOUNDED Zero_Vox TO LAS_oneTID_Vx (UPDATE LAS_Vox_Zero)
      if(nrow(LAS_oneTID_VoxZero@data) > 0){ 
        
        LAS_oneTID_VoxZero@data$TreeID <- as.integer(Unique_TID[VC])
        LAS_oneTID_Vx@data <- rbind(LAS_oneTID_Vx@data, LAS_oneTID_VoxZero@data)
        
        Index_1 <- which(LAS_Vox_Zero@data$VoxID %in% LAS_oneTID_VoxZero@data$VoxID)
        LAS_Vox_Zero@data <- LAS_Vox_Zero@data[-Index_1,]
      }
      
      ###########################
      # CLUSTERING LAS_oneTID_Vx (VoxR::distance_clustering)
      ###########################
      Index_Col1 <- match(c("X", "Y", "Z"), colnames(LAS_oneTID_Vx@data))
      
      if(nrow(LAS_oneTID_Vx@data) > 1){
        Vx_TID_Cl <- VoxR::distance_clustering(data = LAS_oneTID_Vx@data[,..Index_Col1], d_clust=(Para_Vx_R_Res*2 + 0.01)) 
      }else{
        Vx_TID_Cl <-LAS_oneTID_Vx@data
        Vx_TID_Cl <- cbind(Vx_TID_Cl[,..Index_Col1], Cluster = 1)
      }
      colnames(Vx_TID_Cl)[4] <- "VoxCl"
      
      # UPDATE LAS_oneTID_Vx
      Vx_TID_Cl$VoxCl[which(Vx_TID_Cl$VoxCl >0)] <-  Max_Vox_Cl +  Vx_TID_Cl$VoxCl[which(Vx_TID_Cl$VoxCl >0)]
      Vx_TID_Cl <- as.data.frame(Vx_TID_Cl)
      LAS_oneTID_Vx@data$VoxCl <- as.integer(Vx_TID_Cl$VoxCl)
      
      Max_Vox_Cl<- max(Vx_TID_Cl$VoxCl)
      
      # BOUNDARY OF ALL CLUSTERS (LAS_oneTID_Vx)
      Bound_oneTID_allCl <- as.data.frame(LAS_oneTID_Vx@data %>%
                                            dplyr::group_by(VoxCl) %>%
                                            dplyr::summarise(TreeID = unique(TreeID),
                                                             Vox_Count = length(Z),
                                                             MinX = min(X)-2*Para_Vx_R_Res,
                                                             MaxX = max(X)+2*Para_Vx_R_Res,
                                                             MinY = min(Y)-2*Para_Vx_R_Res,
                                                             MaxY = max(Y)+2*Para_Vx_R_Res,
                                                             MinZ = min(Z)-2*Para_Vx_R_Res,
                                                             MaxZ = max(Z)+2*Para_Vx_R_Res,
                                                             Range = range(Z)[2] - range(Z)[1],
                                                             .groups = 'drop'))
      
      ########################################### 7b
      # REMOVE SMALL CLUSTERS (PUT INTO LAS_ZERO)
      ########################################### 7b
      Zero_Sm_VoxCl <- which(Bound_oneTID_allCl$VoxCl == 0 | Bound_oneTID_allCl$Vox_Count < Para_Rm_VCnt)
      
      if(length(Zero_Sm_VoxCl) > 0){

        LAS_oneTID_Vx_Zero <- filter_poi(LAS_oneTID_Vx, VoxCl %in% Bound_oneTID_allCl$VoxCl[Zero_Sm_VoxCl])
        LAS_oneTID_Vx_Zero@data$TreeID <- as.integer(0)
        LAS_Vox_Zero@data <- rbind(LAS_Vox_Zero@data, LAS_oneTID_Vx_Zero@data)  # UPDATE LAS_Vox_Zero
        
        # UPDATING LAS_TID 
        Index_1 <- which(LAS_TID@data$VoxID %in% LAS_oneTID_Vx_Zero@data$VoxID)
        if(length(Index_1) > 0){  
          LAS_TID@data$TreeID[Index_1] <- as.integer(0)
          LAS_Zero@data<- rbind(LAS_Zero@data, LAS_TID@data[Index_1,])
          LAS_TID@data<- LAS_TID@data[-Index_1,]
        }
  
        # UPDATING LAS_VOX_TID 
        Index_1 <- which(LAS_Vox_TID@data$VoxID %in% LAS_oneTID_Vx_Zero@data$VoxID)
        if(length(Index_1) > 0){ # ADDED THIS WHEN MOVING FROM V1.6 to 3.0.5
          LAS_Vox_TID@data$TreeID[Index_1] <- as.integer(0)
          LAS_Vox_TID@data<- LAS_Vox_TID@data[-Index_1,]
        }
        
        # UPDATING LAS_oneTID_Vx and Bound_oneTID_allCl
        LAS_oneTID_Vx <- filter_poi(LAS_oneTID_Vx, !(VoxCl %in% Bound_oneTID_allCl$VoxCl[Zero_Sm_VoxCl])) #VoxCl > 0
        Bound_oneTID_allCl <- Bound_oneTID_allCl[-Zero_Sm_VoxCl,]
        }
      
      #################################################
      # ASSESS LARGE CLUSTERS FOR MOVING INTO OTHER TID
      #################################################
      
      # IF THERE ARE LARGE CLUSTERS LEFT FOR TID
      if(nrow(Bound_oneTID_allCl) > 0){
        # ORDER RETAINED CLUSTERS USING RANGE (LARGEST RANGE GETS SPECIAL PROCESS: ONLY ASSESS NEIGHBOURS BELOW IT....)
        Bound_oneTID_allCl<- Bound_oneTID_allCl[with(Bound_oneTID_allCl, order(-Range)), ]
        
        # REMOVING CLUSTER WITH LARGEST RANGE IF IT BEGINS BELOW THE MID HEIGHT OF WHOLE SITE
        if(Bound_oneTID_allCl$MinZ[1]< HalfHeight_LAS){
          Index_ExcludeMoving <- 1
          Bound_oneTID_allCl_TODO <- Bound_oneTID_allCl[-Index_ExcludeMoving,]
        }else{
          Index_ExcludeMoving <- 0
          Bound_oneTID_allCl_TODO <- Bound_oneTID_allCl
        }
        
        # IF THERE ARE CLUSTERS ASSESSED FOR MOVING (NOT INCLUDING LARGEST ONE THAT HAS POINTS BELOW HalfHeight_LAS)
        if(nrow(Bound_oneTID_allCl_TODO)>0){
          
          ####################################################
          # BOUNDING BOX FOR LAS_oneTID_Vx AND SUUROUNDING TID
          ####################################################
          if(length(unique(LAS_oneTID_Vx@data$TreeID))>1){browser()}
          # BOUNDARY OF ONE TID (UPDATED WITH ZERO LAS INTEGRATED)
          Bound_oneTID <- as.data.frame(LAS_oneTID_Vx@data %>%
                                          dplyr::summarise(TreeID = unique(TreeID),
                                                           Vox_Count = length(Z),
                                                           MinX = min(X)-2*Para_Vx_R_Res,
                                                           MaxX = max(X)+2*Para_Vx_R_Res,
                                                           MinY = min(Y)-2*Para_Vx_R_Res,
                                                           MaxY = max(Y)+2*Para_Vx_R_Res,
                                                           MinZ = min(Z)-2*Para_Vx_R_Res,
                                                           MaxZ = max(Z)+2*Para_Vx_R_Res, 
                                                           .groups = 'drop'))
          
          # GET NEIGHBOUR LAS SURROUNDING "Subject" TID IN BOX (+- 0.2 +- 0.01).
          LAS_surrOneTID_Vx <- filter_poi(LAS_Vox_TID, X >= Bound_oneTID$MinX - 0.01 &
                                           X <= Bound_oneTID$MaxX + 0.01 &
                                           Y >= Bound_oneTID$MinY - 0.01 &
                                           Y <= Bound_oneTID$MaxY + 0.01 &
                                           Z >= Bound_oneTID$MinZ - 0.01 &
                                           Z <= Bound_oneTID$MaxZ + 0.01 &
                                           !(TreeID == Unique_TID[VC]))
          
          ################################################### 7c
          # LOOP THROUGH LARGE VoxCl FOR MOVING TO OTHER TID. 
          ################################################### 7c
          
          Binary_EndWhile <- 1
          while(Binary_EndWhile == 1){
            if(nrow(Bound_oneTID_allCl_TODO) > 0){
              
              #############################
              # LOOP THROUGH EACH TID VoxCl
              #############################
              Cluster_Update <- c()
              for(FP in 1:nrow(Bound_oneTID_allCl_TODO)){
                
                ####################
                # GENERATE LAS FILES
                ####################
                
                # CREATE EMPTY LAS FILES (NEEDS TO BE EMPTY AT START OF EACH LOOP SO IT DOESN'T USE PREVIOUS ONES)
                LAS_surrOneCl <- LAS_Empty                                                                             # POINTS WITHIN BOUNDING BOX OF ONE CLUSTER FROM ONE TID
                LAS_oneCl_oneTID_Vx <- filter_poi(LAS_oneTID_Vx, VoxCl %in% Bound_oneTID_allCl_TODO$VoxCl[FP])          # ONE CLUSTER FROM ONE TID
                LAS_rmOneCl_oneTID_Vx <- filter_poi(LAS_oneTID_Vx, !(VoxCl %in% Bound_oneTID_allCl_TODO$VoxCl[FP]))     # ONE TID WITHOUT THE ONE SUBJECT CLUSTER
                
                # LAS SURROUNDING ONE CLUSTER 
                if(nrow(LAS_surrOneTID_Vx@data) > 0){ 
                  LAS_surrOneCl <- filter_poi(LAS_surrOneTID_Vx, X >= Bound_oneTID_allCl_TODO$MinX[FP] - 0.01 &
                                               X <= Bound_oneTID_allCl_TODO$MaxX[FP] + 0.01 &
                                               Y >= Bound_oneTID_allCl_TODO$MinY[FP] - 0.01 &
                                               Y <= Bound_oneTID_allCl_TODO$MaxY[FP] + 0.01 &
                                               Z >= Bound_oneTID_allCl_TODO$MinZ[FP] - 0.01 &
                                               Z <= Bound_oneTID_allCl_TODO$MaxZ[FP] + 0.01 &
                                               !(TreeID == Unique_TID[VC]))
                  
                  # MERGING rmOneCl AND surrOneCl LAS DATA
                  if(nrow(LAS_rmOneCl_oneTID_Vx@data) > 0){
                    LAS_surrOneCl_rmOneCl <- LAS_rmOneCl_oneTID_Vx # ONLY TID (WITHOUT SUBJECT CLUSTER)
                    if(nrow(LAS_surrOneCl@data) > 0){ 
                      LAS_surrOneCl_rmOneCl@data <- rbind(LAS_surrOneCl_rmOneCl@data, LAS_surrOneCl@data) # BOTH SURROUND AND TID (WITHOUT SUBJECT CLUSTER)
                    }
                  }else{
                    if(nrow(LAS_surrOneCl@data) > 0){ 
                      LAS_surrOneCl_rmOneCl <- LAS_surrOneCl # ONLY SURROUND
                    }else{
                      LAS_surrOneCl_rmOneCl <- LAS_Empty # NONE
                    }
                  }
                  
                  ############################################################################## 7d
                  # VOXEL DISTANCE CALCULATION FOR LAS_oneCl_oneTID_Vx AND LAS_surrOneCl_rmOneCl
                  ############################################################################## 7d
                  
                  if(nrow(LAS_surrOneCl_rmOneCl@data) > 0){ 
                    
                    # GET CLOSE VOXEL SUMMARY BETWEEN oneCl and surrOneCl_rmOneCl
                    Index_Col1 <- match(c("X", "Y", "Z", "TreeID", "VoxID"), colnames(LAS_oneCl_oneTID_Vx@data))
                    XYZ_oneCl <- as.data.frame(LAS_oneCl_oneTID_Vx@data)[,Index_Col1]
                    Index_Col2 <- match(c("X", "Y", "Z", "TreeID", "VoxID"), colnames(LAS_surrOneCl_rmOneCl@data))
                    XYZ_surrOneCl <- as.data.frame(LAS_surrOneCl_rmOneCl@data)[,Index_Col2]
                    
                    # DISTANCE BETWEEN ALL VOXELS IN TWO LAS FILES (THE VOXELS ARE ONLY WITHIN BBox). MAX DISTANCE SPECIFIED
                    Dist_oneCl_surrOneCl <- VOX_SLICE_DIST_FUN(XYZ_oneCl,
                                                               XYZ_surrOneCl,
                                                               Max_Dist = Para_EDT_V, 
                                                               Para_Vx_R_Res= Para_Vx_R_Res,
                                                               Slice_Thickness= Para_Sl_Z)
                    
                    Summary_closeSurr <- Dist_oneCl_surrOneCl$Summary_closeSurr
                    
                    # THERE ARE CLOSE SURR, WORK OUT WHICH ONE THEY SHOULD MERGE TO
                    if(!is.null(Summary_closeSurr)){
                      Summary_closeSurr <- Summary_closeSurr[order(Summary_closeSurr$MinZ),]
                      
                      ############################################
                      # GET WHOLE TIDs WITHIN HEIGHT RANGE OF BBOX (i.e. IN surrOneCl AND SUBJECT TID)
                      ############################################
                      LAS_surrOneCl_rmOneCl_allTID <- NULL
                      for(SW in 1:nrow(Summary_closeSurr)){
                        if(is.null(LAS_surrOneCl_rmOneCl_allTID)){
                          LAS_surrOneCl_rmOneCl_allTID <- filter_poi(LAS_Vox_TID, TreeID == Summary_closeSurr$TID_Surr[SW] &
                                                                      Z >= Summary_closeSurr$MinZ[SW] - Para_Sl_Z &
                                                                      Z <= Summary_closeSurr$MaxZ[SW] + Para_Sl_Z)
                        }else{
                          LAS_surrOneCl_rmOneCl_oneTID <- filter_poi(LAS_Vox_TID, TreeID == Summary_closeSurr$TID_Surr[SW] &
                                                                      Z >= Summary_closeSurr$MinZ[SW] - Para_Sl_Z &
                                                                      Z <= Summary_closeSurr$MaxZ[SW] + Para_Sl_Z)
                          if(nrow(LAS_surrOneCl_rmOneCl_oneTID@data) > 0){ 
                            LAS_surrOneCl_rmOneCl_allTID@data <- rbind(LAS_surrOneCl_rmOneCl_allTID@data , LAS_surrOneCl_rmOneCl_oneTID@data)
                          }
                        }
                      }
                      
                      ############################################################################################
                      # IF THERE IS MORE THAN ONE CLOSEST SURROUNDING TID, UNDERTAKE SLICING PROCEDURE TO ALLOCATE
                      ############################################################################################
                      
                      if(nrow(Summary_closeSurr) > 1){
                        
                        # MIN Z FOR SURROUNDING TID
                        minZ_oneCl <- Summary_closeSurr$MinZ[1]
                        
                        # SLICE INCREMENTS
                        SeqMinMax <- c(((floor(minZ_oneCl/Para_Vx_R_Res)*Para_Vx_R_Res)-Para_Sl_Z),
                                       ceiling(max(LAS_oneCl_oneTID_Vx$Z)/Para_Vx_R_Res)*Para_Vx_R_Res + Para_Sl_Z)
                        Z_Inc_Slice <- seq(min(SeqMinMax), max(SeqMinMax), Para_Sl_Z)
                        
                        # TRIM SURROUNDS SO ITS WITHIN oneCl SLICE INCREMENTS
                        LAS_surrOneCl_rmOneCl_allTID <- filter_poi(LAS_surrOneCl_rmOneCl_allTID, Z>= min(Z_Inc_Slice) & Z<max(Z_Inc_Slice))
                        
                        if(nrow(LAS_surrOneCl_rmOneCl_allTID@data) > 0){
                          
                          #############################################
                          ############################################# 7e
                          # LOOP: SLICES OF SURROUNDING TID WITHIN BBOX
                          ############################################# 7e
                          #############################################
                          
                          Slice_SurrBBox <- SLICE_SURROUND_BBOX_FUN(LAS_surrOneCl_rmOneCl_allTID, LAS_Empty, Unique_TID = Unique_TID[VC],  
                                                                    Para_EDT_V, Z_Inc_Slice, Para_Vx_R_Res, Para_Sl_Z, 
                                                                    LAS_or_Vox = "Vox")
                          LAS_surrOneCl_rmOneCl_allTID <- Slice_SurrBBox$LAS_surrOneCl_rmOneCl_allTID
                          
                          # STORE UPDATED CHANGES IN LAS 
                          if(nrow(LAS_Update_Changes@data) > 0){ 
                            LAS_Update_Changes@data <- rbind(LAS_Update_Changes@data , LAS_surrOneCl_rmOneCl_allTID@data)
                          }else{
                            LAS_Update_Changes <- LAS_surrOneCl_rmOneCl_allTID
                          }
                          
                          ########################################################
                          # IF CLUSTER HAS NO VOXELS TO SLICE (RARE TROUBLE SHOOT)
                          ########################################################  
                        }else{
                          browser() 
                          if(nrow(LAS_Update_Changes@data) > 0){ 
                            LAS_Update_Changes@data <- rbind(LAS_Update_Changes@data , LAS_oneCl_oneTID_Vx@data)
                          }else{
                            LAS_Update_Changes <- LAS_oneCl_oneTID_Vx
                          }
                        }
                        
                        ########################################################
                        # IF ONLY ONE SURR THEN GIVE TID THE NEAREST STEM VALUE
                        ########################################################
                      }else{
                        LAS_oneCl_oneTID_Vx@data$TreeID  <- as.integer(Summary_closeSurr$TID_Surr[1])
                        
                        #  STORE UPDATED CHANGES IN LAS 
                        if(nrow(LAS_Update_Changes@data) > 0){ 
                          LAS_Update_Changes@data <- rbind(LAS_Update_Changes@data , LAS_oneCl_oneTID_Vx@data)
                        }else{
                          LAS_Update_Changes <- LAS_oneCl_oneTID_Vx
                        }
                      }
                      
                      ############
                      # IF NO SURR 
                      ############
                    }else{
                      if(nrow(LAS_Update_Changes@data) > 0){
                        LAS_Update_Changes@data <- rbind(LAS_Update_Changes@data , LAS_oneCl_oneTID_Vx@data)
                      }else{
                        LAS_Update_Changes <- LAS_oneCl_oneTID_Vx
                      }
                    }
                    
                    ########################################
                    # IF NO SURR: (NO LAS_surrOneCl_rmOneCl) 
                    ########################################    
                  }else{
                    if(nrow(LAS_Update_Changes@data) > 0){ 
                      LAS_Update_Changes@data <- rbind(LAS_Update_Changes@data , LAS_oneCl_oneTID_Vx@data)
                    }else{
                      LAS_Update_Changes <- LAS_oneCl_oneTID_Vx
                    }
                  }
                  
                  ###################################################
                  # IF NO SURR AND NO rmOneCl: (NO LAS_surrOneTID_Vx) 
                  ###################################################                   
                }else{
                  # CREATE NEW TID
                  if(is.null(LAS_rmOneCl_oneTID_Vx)){
                    LAS_oneCl_oneTID_Vx@data$TreeID <- as.integer(Max_TID + 1)
                    Max_TID <- Max_TID + 1
                    if(nrow(LAS_Update_Changes@data) > 0){ 
                      LAS_Update_Changes@data <- rbind(LAS_Update_Changes@data , LAS_oneCl_oneTID_Vx@data)
                    }else{
                      LAS_Update_Changes <- LAS_oneCl_oneTID_Vx
                    }
                  }
                }
                
              } # LOOP FP (Loop through each VoxCl)
              
              ########################
              # UPDATE Binary_EndWhile
              ########################
              # UPDATE TO DO VoxCl
              if(length(Cluster_Update) > 0){
                Bound_oneTID_allCl_TODO <- Bound_oneTID_allCl_TODO[-Cluster_Update,]
                Binary_EndWhile <- 1 # REPEAT WHILE LOOP
              }else{
                Binary_EndWhile <- 0  # END WHILE LOOP
              }
            }else{ # NO MORE VoxCl TO DO
              Binary_EndWhile <- 0 # END WHILE LOOP
            }
          } # WHILE LOOP
          
          ######################################### 7
          # MAIN UPDATE (LAS_Vox_TID and LAS_TID)
          ######################################### 7
          
          # # UPDATING CHANGES FOR VC LOOP (LAS_Vox_TID & LAS_TID) 
          if(nrow(LAS_Update_Changes@data) > 0){ 
            Index_1 <- which(LAS_Vox_TID@data$VoxID %in% c(LAS_Update_Changes@data$VoxID))
            Index_2 <- match(LAS_Vox_TID@data$VoxID[Index_1], LAS_Update_Changes@data$VoxID)
            LAS_Vox_TID@data$TreeID[Index_1] <- LAS_Update_Changes@data$TreeID[Index_2]
            
            Index_1 <- which(LAS_TID@data$VoxID %in% c(LAS_Update_Changes@data$VoxID))
            Index_2 <- match(LAS_TID@data$VoxID[Index_1], LAS_Update_Changes@data$VoxID)
            LAS_TID@data$TreeID[Index_1] <- LAS_Update_Changes@data$TreeID[Index_2]
          }
        }# IF THERE ARE CLUSTERS ASSESSED FOR MOVING
      } # IF THERE ARE LARGE VOXELS if(nrow(Bound_oneTID_allCl) > 0)
    } # IF VOX ID STILL EXISTS (SOME GET MOVED IN PREVIOUS ITERATIONS)
  } # LOOP VC THROUGH EACH St
  
  # TIMER   # ptm <- proc.time()
  Time_Taken <- as.vector(round(((proc.time() - ptm)[3]/60),3))
  TIME_Section <- c(TIME_Section, Time_Taken)
  TIME_Section_ID <- c(TIME_Section_ID, 7)
  TIME__Count_Vox_TID <- c(TIME__Count_Vox_TID, length(unique(LAS_Vox_TID@data$TreeID)))
  TIME__Count_TID <- c(TIME__Count_TID, length(unique(LAS_TID@data$TreeID)))
  Output_Table_TIME_Section_ID <- data.frame(TIME_Section_ID, TIME_Section, TIME__Count_Vox_TID, TIME__Count_TID)
  write.csv(Output_Table_TIME_Section_ID, paste(Output_CSV, "/F", Flight_Number ,"_6_TIME_TAKEN.csv",sep=''), row.names = FALSE)
  
  # TABLE OUTPUT TID COUNT
  Code_Number <- c(Code_Number, 7)
  Count_LAS_TIDs <- c(Count_LAS_TIDs , length(unique(LAS_TID@data$TreeID)))
  Count_Vox_TIDs <- c(Count_Vox_TIDs , length(unique(LAS_Vox_TID@data$TreeID)))
  Output_Table_TID_CodeSummary <- data.frame(Code_Number, Count_LAS_TIDs, Count_Vox_TIDs)
  write.csv(Output_Table_TID_CodeSummary, paste(Output_CSV, "/F", Flight_Number ,"_6_Code_Summary.csv",sep=''), row.names = FALSE)
  
  ########
  # BACKUP
  ########
  LAS_Zero_Backup6 <- LAS_Zero
  LAS_Vox_Zero_Backup6 <- LAS_Vox_Zero
  LAS_TID_Backup6 <- LAS_TID
  LAS_Vox_TID_Backup6 <- LAS_Vox_TID
  
  LAS_Zero <- LAS_Zero_Backup6
  LAS_Vox_Zero <- LAS_Vox_Zero_Backup6
  LAS_TID <- LAS_TID_Backup6
  LAS_Vox_TID <- LAS_Vox_TID_Backup6
  print("Backup6")
  
  # QUICK FIX .... this small error came in from previous step
  LAS_Vox_TID <- filter_poi(LAS_Vox_TID, !duplicated(LAS_Vox_TID@data$VoxID))
  
  LAS_Zero@data$TreeID <- as.integer(LAS_Zero@data$TreeID)
  writeLAS(LAS_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Zero_Backup6.laz",sep=''))
  LAS_Vox_Zero@data$TreeID <- as.integer(LAS_Vox_Zero@data$TreeID)
  writeLAS(LAS_Vox_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_Zero_Backup6.laz",sep=''))
  
  LAS_TID@data$TreeID <- as.integer(LAS_TID@data$TreeID)
  writeLAS(LAS_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_TID_Backup6.laz",sep=''))
  
  LAS_Vox_TID@data$TreeID <- as.integer(LAS_Vox_TID@data$TreeID)
  writeLAS(LAS_Vox_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_TID_Backup6.laz",sep=''))
  
  LAS_IDW_Under@data$TreeID <- as.integer(LAS_IDW_Under@data$TreeID)
  writeLAS(LAS_IDW_Under, paste(Output_LAS,"/F", Flight_Number ,"_LAS_IDW_Under_Backup6.laz",sep=''))

  ################################################################################################################################################################ 8
  ################################################################################################################################################################ 8
  # GROUND AND NON_GROUND MERGE
  ################################################################################################################################################################ 8
  ################################################################################################################################################################ 8
  
  LAS_Vox_TID <- merge_spatial(LAS_Vox_TID, R_Under, attribute = "Under")  
  LAS_TID <- merge_spatial(LAS_TID, R_Under, attribute = "Under") 
  LAS_IDW_Under <- merge_spatial(LAS_IDW_Under, R_Under, attribute = "Under")  
  
  for(RR in 1:2){
    ptm <- proc.time()
  ################################################################################# 8a
  # DIFFERENTIATE TID CLOSE TO GROUND (Gnd) AND THOSE NOT CLOSE TO GROUND (Non_Gnd)
  ################################################################################# 8a
  
  # TID BOUNDARY: ALSO RangeZ AND MIN DISTANCE TO UNDER.
  Bound_TID_minDistU <- as.data.frame(LAS_Vox_TID@data %>%
                                        dplyr::group_by(TreeID) %>%
                                        dplyr::summarise(
                                          MinX = min(X),
                                          MaxX = max(X),
                                          MinY = min(Y),
                                          MaxY = max(Y),
                                          MinZ = min(Z),
                                          MaxZ = max(Z),
                                          RangeZ = MaxZ - MinZ,
                                          max_Under = max(Under, na.rm =T),
                                          Under_Closest = Under[which.min(abs(Z - Under))],
                                          Dist_Under = min(abs(Z - Under),  na.rm = TRUE), 
                                          .groups = 'drop'))
  
  Bound_TID_minDistU <- Bound_TID_minDistU[order(Bound_TID_minDistU$Dist_Under),]
  
  # CONSTRAINTS FOR DIFFERENTIATING Gnd FROM Non_Gnd
  Quarter_Forest_Height <- (max(Bound_TID_minDistU$MaxZ)/4)
  Fifth_Forest_Height <- (max(Bound_TID_minDistU$MaxZ)/5)
  Gnd_TID_LgRng <- Bound_TID_minDistU[which(Bound_TID_minDistU$RangeZ > Quarter_Forest_Height &
                                              Bound_TID_minDistU$MinZ < (Fifth_Forest_Height + mean(Bound_TID_minDistU$max_Under))),] # CLUSTERS WITH LARGE Z RANGE AND CLOSE TO GROUND
  Gnd_TID_TchUnder <- Bound_TID_minDistU[which(Bound_TID_minDistU$Dist_Under < Para_EDT_V),]
  Gnd_All <- rbind(Gnd_TID_LgRng, Gnd_TID_TchUnder)
  Gnd_TID_95PTileMinZ <- Bound_TID_minDistU[which(Bound_TID_minDistU$MinZ < quantile(Gnd_All$MinZ, 0.95)),]
  
  # CREATE BBOX DF FOR Gnd AND NonGnd
  Gnd_All_TID <- unique(c(Gnd_All$TreeID, Gnd_TID_95PTileMinZ$TreeID))
  if(length(Gnd_All_TID) > 0){
    Non_Gnd_TID <- Bound_TID_minDistU[-which(Bound_TID_minDistU$TreeID %in% Gnd_All_TID),]
    Gnd_TID <- Bound_TID_minDistU[which(Bound_TID_minDistU$TreeID %in% Gnd_All_TID),]
  }else{
    Gnd_TID <- NULL
    Non_Gnd_TID <- Bound_TID_minDistU
  }
  
  # GENERATE LAS FOR Gnd AND NonGnd
  if(!is.null(Gnd_TID$TreeID)){
    LAS_Gnd_TID <- filter_poi(LAS_Vox_TID, TreeID %in% Gnd_TID$TreeID)
  }else{
    # IF NO GROUND THEN GIVE ONE TID GND (ONE WITH LARGE RANGE)
    Gnd_TID <- Bound_TID_minDistU[which.max(Bound_TID_minDistU$RangeZ),]
    LAS_Gnd_TID <- filter_poi(LAS_Vox_TID, TreeID %in% Gnd_TID$TreeID)
    Non_Gnd_TID <- Bound_TID_minDistU[-which.max(Bound_TID_minDistU$RangeZ),]
  }
  LAS_NonGnd_TID <- filter_poi(LAS_Vox_TID, TreeID %in% Non_Gnd_TID$TreeID)
  
  ########################## 8b
  # LOOP THROUGH EACH NonGnd
  ########################## 8b
  
  # EMPTY DATAFRAMES AND LAS FOR LOOPING THROUGH GROUND POINTS.
  if(exists("All_Summary_oneNG_ZeroBtw_U_CB")){rm(All_Summary_oneNG_ZeroBtw_U_CB) }
  if(exists("All_Closest_Vox_SurrGnd_ZeroWS")){ rm(All_Closest_Vox_SurrGnd_ZeroWS)}
  LAS_allNG_ZeroBtw_U_CB <- LAS_Empty # NULL
  
  Unique_NonGndTID <- Non_Gnd_TID$TreeID  
  for(NG in 1:length(Unique_NonGndTID)){
    
    # GET LAS, CENTROIDS AND BOUND OF NonGndTID
    Bound_oneNonGnd_TID <- Bound_TID_minDistU[which(Bound_TID_minDistU$TreeID == Unique_NonGndTID[NG]),]
    LAS_oneNG <- filter_poi(LAS_NonGnd_TID, TreeID == Unique_NonGndTID[NG])
    
    ######################################### 8c
    # IDENTIFYING CLOSE GnD TO SUBJECT NonGnd
    ######################################### 8c
    
    # GET SurrGND USING NonGnd Boundary
    LAS_SurrGND <- filter_poi(LAS_Gnd_TID, X >= Bound_oneNonGnd_TID$MinX - Para_EDT_V  &
                               X <= Bound_oneNonGnd_TID$MaxX + Para_EDT_V  &
                               Y >= Bound_oneNonGnd_TID$MinY - Para_EDT_V &
                               Y <= Bound_oneNonGnd_TID$MaxY + Para_EDT_V &
                               Z >= Bound_oneNonGnd_TID$MinZ - Para_EDT_V &
                               Z <= Bound_oneNonGnd_TID$MaxZ + Para_EDT_V)
    
    if(nrow(LAS_SurrGND@data) > 0){ 
      # GET CLOSE VOXEL SUMMARY BETWEEN NonGnd_TID and Gnd_TID WITHIN BBOX
      Index_Col1 <- match(c("X", "Y", "Z", "TreeID", "VoxID"), colnames(LAS_oneNG@data))
      XYZ_oneNG <- LAS_oneNG@data[,..Index_Col1]
      Index_Col2 <- match(c("X", "Y", "Z", "TreeID", "VoxID"), colnames(LAS_SurrGND@data))
      XYZ_SurrGND <- LAS_SurrGND@data[,..Index_Col2]
      
      Dist_oneNG_SurrGND <- VOX_SLICE_DIST_FUN(XYZ_oneNG, XYZ_SurrGND,
                                               Max_Dist = Para_EDT_V, 
                                               Para_Vx_R_Res= Para_Vx_R_Res,
                                               Slice_Thickness= Para_Sl_Z)
      Summary_closeSurrGND <- Dist_oneNG_SurrGND$Summary_closeSurr
      
      # IF CLOSE GnD SURROUNDS EXIST FOR MERGING WITH subject NonGnd 
      if(!is.null(Summary_closeSurrGND)){
        
        Summary_closeSurrGND <- Summary_closeSurrGND[order(Summary_closeSurrGND$MinZ),]
        
        ######################################## 8d
        # CONSTRAINT FOR WHETHER TO MERGE TO Gnd
        ######################################## 8d
        
        # IF CLOSE VOXEL IS VERY HIGH UP THE oneNG RANGE THEN DO NOT MOVE IT
        min_oneNG <- range(LAS_oneNG$Z)[1]
        max_oneNG <- range(LAS_oneNG$Z)[2]
        PortZ_ClVox_oneNG <- (max_oneNG-Summary_closeSurrGND$MinZ[1])/(max_oneNG-min_oneNG) 
        Distance_Below_CloseVoxMin <- Summary_closeSurrGND$MinZ[1]-min_oneNG 
        
        # MOVE NonGnd TO Gnd IF CONTRAINT TRUE
        if(PortZ_ClVox_oneNG > 4/5 | Distance_Below_CloseVoxMin < 2){
          
          # ORDER SUMMARY OF CLOSE VOXELS
          Summary_closeSurrGND <- Summary_closeSurrGND[order(Summary_closeSurrGND$MinZ),]
          
          # GET SurrGND FOR ALL TID THAT ARE CLOSE AND WITHIN Z Range 
          LAS_SurrGND_allTID <- NULL  
          for(SW in 1:nrow(Summary_closeSurrGND)){
            LAS_oneSurr <- filter_poi(LAS_Gnd_TID, TreeID == Summary_closeSurrGND$TID_Surr[SW] &
                                       Z >= Summary_closeSurrGND$MinZ[SW] - Para_Sl_Z &    
                                       Z <= Summary_closeSurrGND$MaxZ[SW] + Para_Sl_Z)
            if(is.null(LAS_SurrGND_allTID)){
              LAS_SurrGND_allTID <- LAS_oneSurr
            }else{
              if(nrow(LAS_oneSurr@data) > 0){
                LAS_SurrGND_allTID@data <- rbind(LAS_SurrGND_allTID@data , LAS_oneSurr@data)
              }
            }
          } # SW LOOP
          
          ######################## 8e
          # MERGING Non_Gnd TO Gnd
          ######################## 8e
          
          # IF THERE IS MORE THAN ONE SurrGND, oneNG GETS SPLIT AMONGST SurrGND TID
          if(nrow(Summary_closeSurrGND) > 1){ # MERGE
            
            # UPDATE TID for oneNG UPTO MIN HEIGHT (Above Para_EDT_V) OF SECOND SurrGND
            Top_Z_oneNG <- Summary_closeSurrGND$MinZ[1]

            # MERGE BOTH TOGETHER BEFORE SLICING THROUGH THEM
            LAS_oneNG_closeSurrGND <- LAS_oneNG
            LAS_oneNG_closeSurrGND@data <- rbind(LAS_oneNG_closeSurrGND@data, LAS_SurrGND_allTID@data)
            
            ############################################################################################################################
            ############################################################################################################################
            
            ################### 8f
            # SLICING PROCEDURE
            ################### 8f
            
            # SLICE FROM SECOND TID BotCloseVoxZ TO TOP OF oneNG (i.e. subject VoxCl)
            SeqMinMax <- c(((floor(Top_Z_oneNG/Para_Vx_R_Res)*Para_Vx_R_Res)-Para_Sl_Z),
                           ceiling(max(LAS_oneNG$Z)/Para_Vx_R_Res)*Para_Vx_R_Res + Para_Sl_Z)
            Z_Inc_Slice <- seq(min(SeqMinMax), max(SeqMinMax), Para_Sl_Z)
            
            # TEST_IF SAMVE
            pre_LAS_oneNG_closeSurrGND <- LAS_oneNG_closeSurrGND
            
            Slice_SurrBBox <- SLICE_SURROUND_BBOX_FUN(LAS_oneNG_closeSurrGND, LAS_Empty, Unique_TID = Unique_NonGndTID[NG], 
                                                      Para_EDT_V, Z_Inc_Slice, Para_Vx_R_Res, Para_Sl_Z, 
                                                      LAS_or_Vox = "Vox")
            LAS_oneNG_closeSurrGND <- Slice_SurrBBox$LAS_surrOneCl_rmOneCl_allTID 
            
          }else{ # MERGE
            
            ######################################
            # IF ONLY ONE SURR, oneNG GETS GND TID
            ######################################
            LAS_oneNG@data$TreeID  <- as.integer(Summary_closeSurrGND$TID_Surr[1])
            LAS_oneNG_closeSurrGND <- LAS_oneNG
            LAS_oneNG_closeSurrGND@data <- rbind(LAS_oneNG_closeSurrGND@data, LAS_SurrGND_allTID@data)
            
          }
          
        }else{ # Not MERGED
          # CLOSE VOXELS TOO HIGH UP NonGnd
          LAS_oneNG_closeSurrGND <- LAS_oneNG
        }
      }else{ # Not MERGED
        # SUMMARY GENERATES NO CLOSE SurrGND
        LAS_oneNG_closeSurrGND <- LAS_oneNG
      }
    }else{ # Not MERGED
      # LAS_SurrGND DOESN'T EXISTS
      LAS_oneNG_closeSurrGND <- LAS_oneNG
    }
    
    #############################################################################################################################
    #############################################################################################################################
    
    ##########################
    # UPDATE OF Gnd and NonGnd
    ##########################
    
    if(nrow(LAS_oneNG_closeSurrGND@data) > 0){
      # UPDATING LAS_Gnd_TID   
      Index_1 <- which(LAS_Gnd_TID@data$VoxID %in%  LAS_oneNG_closeSurrGND@data$VoxID )
      if(length(Index_1) > 0){
        Index_2 <- match(LAS_Gnd_TID@data$VoxID[Index_1],  LAS_oneNG_closeSurrGND@data$VoxID )
        LAS_Gnd_TID@data$TreeID[Index_1] <- as.integer(LAS_oneNG_closeSurrGND@data$TreeID[Index_2])
      }
      
      #REMOVE LAS_oneNG from LAS_NonGnd_TID
      LAS_NonGnd_TID@data <- as.data.table(anti_join(LAS_NonGnd_TID@data, LAS_oneNG_closeSurrGND@data, by=c("VoxID")))#c("X", "Y", "Z")))
      
      # MOVE LAS_oneNG INTO LAS_Gnd_TID
      LAS_Move_NonGnd_2_Gnd <- filter_poi(LAS_oneNG_closeSurrGND, VoxID %in% LAS_oneNG@data$VoxID )
      if(nrow(LAS_Move_NonGnd_2_Gnd@data) > 0){
        LAS_Gnd_TID@data <- rbind(LAS_Gnd_TID@data, LAS_Move_NonGnd_2_Gnd@data)
      }
    }
    print(paste(NG, " out of ", length(Unique_NonGndTID)))
  } # END NG LOOP
  
  # UPDATING  LAS_Vox_TID
  LAS_Vox_TID <- LAS_Gnd_TID
  
  # UPDATING TreeID field in LAS_TID
  Index_1 <- which(LAS_TID@data$VoxID %in%  LAS_Vox_TID@data$VoxID )
  Index_2 <- match(LAS_TID@data$VoxID[Index_1],  LAS_Vox_TID@data$VoxID )
  LAS_TID@data$TreeID[Index_1] <- as.integer(LAS_Vox_TID@data$TreeID[Index_2])
  LAS_TID@data$TreeID <- as.integer(LAS_TID@data$TreeID)
  
  if(RR == 1){
    LAS_Zero_Backup7 <- LAS_Zero
    LAS_Vox_Zero_Backup7 <- LAS_Vox_Zero
    LAS_TID_Backup7 <- LAS_TID
    LAS_Vox_TID_Backup7 <- LAS_Vox_TID 
    
    LAS_Zero <- LAS_Zero_Backup7
    LAS_Vox_Zero <- LAS_Vox_Zero_Backup7
    LAS_TID <- LAS_TID_Backup7
    LAS_Vox_TID <- LAS_Vox_TID_Backup7
    print("Backup7")

    # OUTPUT BACKUP OF THIS STAGE ...
    LAS_Zero@data$TreeID <- as.integer(LAS_Zero@data$TreeID)
    writeLAS(LAS_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Zero_Backup7.laz",sep=''))
    LAS_Vox_Zero@data$TreeID <- as.integer(LAS_Vox_Zero@data$TreeID)
    writeLAS(LAS_Vox_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_Zero_Backup7.laz",sep=''))
    LAS_TID@data$TreeID <- as.integer(LAS_TID@data$TreeID)
    writeLAS(LAS_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_TID_Backup7.laz",sep=''))
    LAS_Vox_TID@data$TreeID <- as.integer(LAS_Vox_TID@data$TreeID)
    writeLAS(LAS_Vox_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_TID_Backup7.laz",sep=''))
    LAS_IDW_Under@data$TreeID <- as.integer(LAS_IDW_Under@data$TreeID)
    writeLAS(LAS_IDW_Under, paste(Output_LAS,"/F", Flight_Number ,"_LAS_IDW_Under_Backup7.laz",sep=''))
    
    # TIMER   # ptm <- proc.time()
    Time_Taken <- as.vector(round(((proc.time() - ptm)[3]/60),3))
    TIME_Section <- c(TIME_Section, Time_Taken)
    TIME_Section_ID <- c(TIME_Section_ID, 8)
    TIME__Count_Vox_TID <- c(TIME__Count_Vox_TID, length(unique(LAS_Vox_TID@data$TreeID)))
    TIME__Count_TID <- c(TIME__Count_TID, length(unique(LAS_TID@data$TreeID)))
    Output_Table_TIME_Section_ID <- data.frame(TIME_Section_ID, TIME_Section, TIME__Count_Vox_TID, TIME__Count_TID)
    write.csv(Output_Table_TIME_Section_ID, paste(Output_CSV, "/F", Flight_Number ,"_8_TIME_TAKEN.csv",sep=''), row.names = FALSE)
    
    # TABLE OUTPUT TID COUNT
    Code_Number <- c(Code_Number, 8)
    Count_LAS_TIDs <- c(Count_LAS_TIDs , length(unique(LAS_TID@data$TreeID)))
    Count_Vox_TIDs <- c(Count_Vox_TIDs , length(unique(LAS_Vox_TID@data$TreeID)))
    Output_Table_TID_CodeSummary <- data.frame(Code_Number, Count_LAS_TIDs, Count_Vox_TIDs)
    write.csv(Output_Table_TID_CodeSummary, paste(Output_CSV, "/F", Flight_Number ,"_8_Code_Summary.csv",sep=''), row.names = FALSE)
    
  }else{ 
    
    LAS_Zero_Backup88 <- LAS_Zero
    LAS_Vox_Zero_Backup88 <- LAS_Vox_Zero
    LAS_TID_Backup88 <- LAS_TID
    LAS_Vox_TID_Backup88 <- LAS_Vox_TID 
    
    LAS_Zero <- LAS_Zero_Backup88
    LAS_Vox_Zero <- LAS_Vox_Zero_Backup88
    LAS_TID <- LAS_TID_Backup88
    LAS_Vox_TID <- LAS_Vox_TID_Backup88
    print("Backup88")
    
    # OUTPUT BACKUP OF THIS STAGE ...
    LAS_Zero@data$TreeID <- as.integer(LAS_Zero@data$TreeID)
    writeLAS(LAS_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Zero_Backup88.laz",sep=''))
    LAS_Vox_Zero@data$TreeID <- as.integer(LAS_Vox_Zero@data$TreeID)
    writeLAS(LAS_Vox_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_Zero_Backup88.laz",sep=''))
    LAS_TID@data$TreeID <- as.integer(LAS_TID@data$TreeID)
    writeLAS(LAS_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_TID_Backup88.laz",sep=''))
    LAS_Vox_TID@data$TreeID <- as.integer(LAS_Vox_TID@data$TreeID)
    writeLAS(LAS_Vox_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_TID_Backup88.laz",sep=''))
    LAS_IDW_Under@data$TreeID <- as.integer(LAS_IDW_Under@data$TreeID)
    writeLAS(LAS_IDW_Under, paste(Output_LAS,"/F", Flight_Number ,"_LAS_IDW_Under_Backup88.laz",sep=''))
    
    # TIMER   # ptm <- proc.time()
    Time_Taken <- as.vector(round(((proc.time() - ptm)[3]/60),3))
    TIME_Section <- c(TIME_Section, Time_Taken)
    TIME_Section_ID <- c(TIME_Section_ID, 88)
    TIME__Count_Vox_TID <- c(TIME__Count_Vox_TID, length(unique(LAS_Vox_TID@data$TreeID)))
    TIME__Count_TID <- c(TIME__Count_TID, length(unique(LAS_TID@data$TreeID)))
    Output_Table_TIME_Section_ID <- data.frame(TIME_Section_ID, TIME_Section, TIME__Count_Vox_TID, TIME__Count_TID)
    write.csv(Output_Table_TIME_Section_ID, paste(Output_CSV, "/F", Flight_Number ,"_88_TIME_TAKEN.csv",sep=''), row.names = FALSE)
    
    # TABLE OUTPUT TID COUNT
    Code_Number <- c(Code_Number, 88)
    Count_LAS_TIDs <- c(Count_LAS_TIDs , length(unique(LAS_TID@data$TreeID)))
    Count_Vox_TIDs <- c(Count_Vox_TIDs , length(unique(LAS_Vox_TID@data$TreeID)))
    Output_Table_TID_CodeSummary <- data.frame(Code_Number, Count_LAS_TIDs, Count_Vox_TIDs)
    write.csv(Output_Table_TID_CodeSummary, paste(Output_CSV, "/F", Flight_Number ,"_88_Code_Summary.csv",sep=''), row.names = FALSE)
    }
    } 
    
  ################################################################################################################################ 9
  ################################################################################################################################ 9
  # REMOVE UNDERSTOREY
  ################################################################################################################################ 9
  ################################################################################################################################ 9

  ptm <- proc.time()
  
  # RASTER OF UNDERSTOREY TOP
  R_Under<-   grid_metrics(LAS_IDW_Under, max(Z), res = Para_Vx_R_Res)
  R_Under[is.na(R_Under[])] <- 0
  
  # BOUNDARY AND SUMMARY OF TID 
  LAS_Vox_TID <- merge_spatial(LAS_Vox_TID, R_Under, attribute = "Under")  
  Bound_allTID <- as.data.frame(LAS_Vox_TID@data %>%
                                  dplyr::group_by(TreeID) %>%
                                  dplyr::summarise(MinX = min(X),
                                                   MaxX = max(X),
                                                   WideX = max(X)-min(X),
                                                   MinY = min(Y),
                                                   MaxY = max(Y),
                                                   WideY = max(Y)-min(Y),
                                                   maxWIDE = max(WideX, WideY),
                                                   MinZ = min(Z),
                                                   MaxZ = max(Z),
                                                   HeightZ = MaxZ -MinZ,
                                                   Min_ZminusUnder = min(Z-Under, na.rm = TRUE),
                                                   CountTouchUnder = length(which((Z-Under) < Para_Vx_R_Res)), 
                                                   .groups = 'drop'))
  
  # TID TOUCHING UNDERSTOREY AND NEEDS TO BE LARGE STEM.
  uniqueTID_CloseUnd_LgZ <- Bound_allTID$TreeID[which(Bound_allTID$CountTouchUnder > 0 &
                                                        Bound_allTID$MaxZ > median(Bound_allTID$MaxZ) )]
  
  
  LAS_All_Move<- NULL
  if(length(uniqueTID_CloseUnd_LgZ) > 0){
    for(TU in 1:length(uniqueTID_CloseUnd_LgZ)){
      LAS_oneTID <- filter_poi(LAS_Vox_TID , TreeID ==  uniqueTID_CloseUnd_LgZ[TU])
      
      ########################## 8
      # PROCESS STAND ATTRIBUTES
      ########################## 8
      Att_oneTID <- LAS_oneTID@data[, GET_TREE_ATTRIBUTES_FUN(X, Y, Z, TreeID, Vox_PntCnt, Para_Sl_Z, Unique_Poly_ID = 0), by = "TreeID"] 
      # GETTING OUTPUT FOR TREE ATTRIBUTE FUNCTION (GET_TREE_ATTRIBUTES_FUN)
      Poly_WS_Sl <-unlist(Att_oneTID$Poly_eachWS_All_SPDF)[[1]]
      Centroid_XYZ <-unlist(Att_oneTID$Points_Centroids_XYZ)[[1]]
      
      ######################## 8
      # CANOPY BASE ASSESSMENT
      ######################## 8
      Profile_oneTID <- GAP_DENSITY_FUN(LAS_oneTID$Z,
                                        Para_BW = Para_Canopy_Base_BW,
                                        Para_Threshold_Percent = Para_Canopy_Base_Thresh ,
                                        Plot = "No",
                                        TreeID = uniqueTID_CloseUnd_LgZ[TU])  
      
      # CANOPY BASE (WORK AROUND FOR WHEN ONE OF THEM IS INF)
      OneInf <- which(c(Profile_oneTID$End_Largest_Gap, Profile_oneTID$Start_Largest_Gap) != Inf)
      if(length(OneInf) != 2){
        C_Base_Height <- c(Profile_oneTID$End_Largest_Gap, Profile_oneTID$Start_Largest_Gap)[OneInf]
      }else{
        C_Base_Height <- max(Profile_oneTID$End_Largest_Gap, Profile_oneTID$Start_Largest_Gap)
      }
      
      # GET SUMMARY OF EACH TID AND ORDER WITH MINZ AT TOP # LAS_TID
      Poly_WS_Slice_DF <- as.data.frame(Poly_WS_Sl)
      Summary_Trunk_Range <- as.data.frame(Poly_WS_Slice_DF %>%
                                             dplyr::group_by(PolyHeight) %>%
                                             dplyr::summarise(Count_Poly = length(WSGrid_Count), 
                                                              .groups = 'drop'))
      
      Summary_Below_C_Base <- Summary_Trunk_Range[which(Summary_Trunk_Range$PolyHeight < C_Base_Height), ]
      
      # GET AVERAGE POLY COUNT IN TOP HALF OF STEM AND ALL SLICES WITH MORE THAN THAT VALUE IN BOTTOM HALF
      BelowCBase_TopH_Avg <- floor(mean(Summary_Below_C_Base$Count_Poly[floor(nrow(Summary_Below_C_Base)/2):nrow(Summary_Below_C_Base)]))
      BelowCBase_AboveAvg_BotH <- which(Summary_Below_C_Base$Count_Poly[1:floor(nrow(Summary_Below_C_Base)/2)] > BelowCBase_TopH_Avg)
      
      if(length(BelowCBase_AboveAvg_BotH) > 0){
        for(RTU in rev(1:max(BelowCBase_AboveAvg_BotH))){
          
          AboveSl_DF <- Summary_Below_C_Base[RTU + 1,]
          AboveSl_all <- Poly_WS_Slice_DF[which(Poly_WS_Slice_DF$PolyHeight ==  AboveSl_DF$PolyHeight),]
          AboveSl_allXY <- AboveSl_all[,which(colnames(AboveSl_all) %in% c("X_Cent", "Y_Cent"))]
          
          oneSl_DF <- Summary_Below_C_Base[RTU,]
          oneSl_all <- Poly_WS_Slice_DF[which(Poly_WS_Slice_DF$PolyHeight ==  oneSl_DF$PolyHeight),]
          oneSl_allXY <- oneSl_all[,which(colnames(oneSl_all) %in% c("X_Cent", "Y_Cent"))]
          
          # DO rdist BETWEEN ABOVE CENTROID AND SLICE CENTROID. THOSE THAT ARE NOT ASSIGNED TO "STEM" GO INTO UNDERSTOREY
          Distance_Above_Below_Miss <- rdist(oneSl_allXY,
                                             AboveSl_allXY) 
          
          # GET INDEX & DISTANCE OF GRID CELL THAT IS CLOSEST
          Closest_Index <- apply(Distance_Above_Below_Miss, 1, which.min)
          Closest_Distance <- apply(Distance_Above_Below_Miss, 1, min)
          Index_Far_Centroids <- which(Closest_Distance > Para_EDT_V)
          if(length(Index_Far_Centroids) > 0){
            
            # Extract LiDAR BASED ON PROXIMITY TO CENTROIDS
            LAS_oneTU_oneSl <- filter_poi(LAS_oneTID, Z >= Summary_Below_C_Base$PolyHeight[RTU ]  &
                                           Z <= Summary_Below_C_Base$PolyHeight[RTU + 1])
            
            if(nrow(LAS_oneTU_oneSl@data) > 0){ 
              
              # DISTANCE CALCULATION BETWEEN SLICE CENTROIDS AND SLICE LIDAR
              Index_Col1 <-match(c("X", "Y"), colnames(LAS_oneTU_oneSl@data))
              Dist_Vox_Slice <- rdist(LAS_oneTU_oneSl@data[,..Index_Col1],
                                      oneSl_allXY) # ROW COL
              
              # IDENTIFY CLOSEST VOXELS TO EACH CENTROID
              Closest_Vox<- apply(Dist_Vox_Slice, 1, which.min)
              Closest_Vox_Distance <- apply(Dist_Vox_Slice, 1, min)
              
              Index_Vox_Remove <- which(Closest_Vox %in% Index_Far_Centroids)
              Vox_Remove <- LAS_oneTU_oneSl@data$VoxID[Index_Vox_Remove]
              
              # MOVE TO VOXELS TO UNDERSTOREY
              LAS_Move <- filter_poi(LAS_oneTU_oneSl, VoxID %in% Vox_Remove)
              
              # KEEP TRACK OF ALL REMOVED VOXELS TO REMOVE THEM FROM "LAS_TID" AND "LAS_TID_VOX"
              if(!is.null(LAS_Move)){
                LAS_All_Move<- LAS_Move
              }else{
                if(!is.null(LAS_Move)){
                  LAS_IDW_Under@data <- rbind(LAS_IDW_Under@data, LAS_Move@data)
                  LAS_All_Move@data <- rbind(LAS_All_Move@data, LAS_Move@data)
                }
              }
              
              # REMOVE FAR CENTROIDS FROM Poly_WS_Sl BY MATCHING XY OF CENTROID WITH LINE
              Removed_Centroids <- oneSl_allXY[Index_Far_Centroids,]
              Index_Remove <- match_df(Removed_Centroids, Poly_WS_Slice_DF[,which(colnames(Poly_WS_Slice_DF) %in% c("X_Cent", "Y_Cent"))])
              Index_Remove <-as.numeric(rownames(Index_Remove))
              Poly_WS_Slice_DF <- Poly_WS_Slice_DF[-Index_Remove,]
            }
          }
        }
      }
      print(paste("TU: ", TU, " out of ", length(uniqueTID_CloseUnd_LgZ)))
    }
  }
  
  # UPDATING LAS_Vox_TID AND LAS_TID
  if(!is.null(LAS_All_Move)){ 
    LAS_Vox_TID <- filter_poi(LAS_Vox_TID, !(LAS_Vox_TID@data$VoxID %in%  LAS_All_Move@data$VoxID))
    
    LAS_TID_Move <- filter_poi(LAS_TID, (LAS_TID@data$VoxID %in%  LAS_All_Move@data$VoxID))
    LAS_TID <- filter_poi(LAS_TID, !(LAS_TID@data$VoxID %in%  LAS_All_Move@data$VoxID))
  }
  
  ########
  # BACKUP
  ########
  LAS_Zero_Backup9 <- LAS_Zero
  LAS_Vox_Zero_Backup9 <- LAS_Vox_Zero
  LAS_TID_Backup9 <- LAS_TID
  LAS_Vox_TID_Backup9 <- LAS_Vox_TID
  
  # CLEANING UP ALL TH LAS DATA !!!
  LAS_Zero <- LAS_Zero_Backup9
  LAS_Vox_Zero <- LAS_Vox_Zero_Backup9
  LAS_TID <- LAS_TID_Backup9
  LAS_Vox_TID <- LAS_Vox_TID_Backup9
  print("Backup9")
  
  # OUTPUT LAS 
  LAS_Zero@data$TreeID <- as.integer(LAS_Zero@data$TreeID)
  writeLAS(LAS_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Zero_Backup9.laz",sep=''))
  LAS_Vox_Zero@data$TreeID <- as.integer(LAS_Vox_Zero@data$TreeID)
  writeLAS(LAS_Vox_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_Zero_Backup9.laz",sep=''))
  LAS_TID@data$TreeID <- as.integer(LAS_TID@data$TreeID)
  writeLAS(LAS_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_TID_Backup9.laz",sep=''))
  LAS_Vox_TID@data$TreeID <- as.integer(LAS_Vox_TID@data$TreeID)
  writeLAS(LAS_Vox_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_TID_Backup9.laz",sep=''))
  if(!is.null(LAS_All_Move)){ 
    LAS_All_Move@data$TreeID <- as.integer(LAS_All_Move@data$TreeID)
    writeLAS(LAS_All_Move, paste(Output_LAS,"/F", Flights[f] ,"_LAS_MoveUnder_Backup9.laz",sep=''))
  }
  LAS_IDW_Under@data$TreeID <- as.integer(LAS_IDW_Under@data$TreeID)
  writeLAS(LAS_IDW_Under, paste(Output_LAS,"/F", Flights[f] ,"_LAS_IDW_Under_Backup9.laz",sep=''))

  # REMOVING UNDERSTOREY
  ptm <- proc.time()
  
  # GET RASTER OF CHM UNDERSTOREY, CHM TID, and TID OF CHM
  R_CHM_Under <-   grid_metrics(LAS_IDW_Under, max(Z), res = Para_Vx_R_Res) 
  R_CHM_Under[is.na(R_CHM_Under)] <- mean(getValues(R_CHM_Under), na.rm=T)
  
  R_CHM_TID_Vx <-   grid_metrics(LAS_Vox_TID, max(Z), res = Para_Vx_R_Res) 

  R_TID_maxZ <-   grid_metrics(LAS_Vox_TID, TreeID[which.max(Z)], res = Para_Vx_R_Res)

  # WORK AROUND FOR WHEN EXTENT OF UNDER IS NOT THE SAME AS THE OTHER TWO RASTERS IN R_Stack
  if(!all(extent(R_TID_maxZ) == extent(R_CHM_Under))){
    R_CHM_Under <- crop(R_CHM_Under,extent(R_CHM_TID_Vx))
    }
  
  R_Stack <- stack(R_TID_maxZ, R_CHM_TID_Vx, R_CHM_Under)
  
  R_Stack_DF <- data.frame(Cell_ID = seq(1, ncell(R_CHM_TID_Vx),1) , getValues(R_Stack))
  colnames(R_Stack_DF)[2:4] <- c("TID", "CHM_TID", "CHM_Under")
  
  ################################################################################
  # DEALING WITH STAGS THAT GIVE WEIRD KERNAL DENSITY OUTPUTS FOR UNDERSTOREY TOPS
  ################################################################################
  Merge_Under_Over_CHM <- overlay(R_CHM_TID_Vx, R_CHM_Under,   fun = function(x, y) {
    x[is.na(x[])] <- y[is.na(x[])]
    return(x)
  })
  
  # REPLACE GRID CELLS ABOVE 95th QUANTILE WITH 95th QUANTILE (AVOIDS SITUATIONS WHERE THERE ARE STAGS)
  CHM_95_Quantile <- quantile(getValues(Merge_Under_Over_CHM), probs = c(0.95), na.rm = T)
  Merge_Under_Over_CHM[Merge_Under_Over_CHM > CHM_95_Quantile] <- CHM_95_Quantile
  
  
  # USE DENSITY FUNCTION OF GRID CELL CHM TO DIFFERENTIATE UNDERSTOREY HEIGHT FROM CANOPY
  OneGrid_Gap <- GAP_DENSITY_FUN(na.omit(Merge_Under_Over_CHM),
                                 Para_BW = Para_K_BandWidth,
                                 Para_Threshold_Percent = Para_Gap_Threshold_Percent_U2,
                                 Plot = "No",
                                 TreeID = 1)
  Top_Under_CGaps <- OneGrid_Gap$Start_Largest_Gap
  
  # IDENTIFY LOW LYING TID
  R_Stack_DF_GAPS <- R_Stack_DF[which(R_Stack_DF$CHM_TID < Top_Under_CGaps),]
  LAS_Under_TID <- filter_poi(LAS_Vox_TID, TreeID %in% unique(R_Stack_DF_GAPS$TID))
  
  if(nrow(LAS_Under_TID@data) > 0){ 
    # GET MAX Z OF ALL LOW LYING TID
    maxZ_underTID <- as.data.frame(LAS_Under_TID@data %>%
                                     dplyr::group_by(TreeID) %>%
                                     dplyr::summarise( MaxZ = max(Z), 
                                                       .groups = 'drop'))
    
    # GET ALL TID THAT IS LOW LYING AND ALSO WITH A LOW MAX VALUE
    Under_TID <- maxZ_underTID$TreeID[which(maxZ_underTID$MaxZ < Top_Under_CGaps)]
    LAS_Under_TID_VX_Move <- filter_poi(LAS_Vox_TID, TreeID %in% Under_TID)
    
    if(nrow(LAS_Under_TID_VX_Move@data) > 0){ 
      
      # RASTER OF MIN VALUE FOR LOW LYING TID
      R_Min_UnderMove <-  grid_metrics(LAS_Under_TID_VX_Move, min(Z), res = Para_Vx_R_Res) 
      
      # RASTER OF TID VALUE FOR LOW LYING TID
      R_TID_Min_UnderMove <-  grid_metrics(LAS_Under_TID_VX_Move, TreeID[which.min(Z)], res = Para_Vx_R_Res) 
      
      R_Min_UnderMove_DF <- na.omit(data.frame(getValues(R_TID_Min_UnderMove),
                                               Cell_ID = seq(1, ncell(R_Min_UnderMove),1),
                                               getValues(R_Min_UnderMove)))
      
      # XY COORDINATES OF ALL MIN VALUES
      XYCell_UnderMove<- data.frame(xyFromCell(R_Min_UnderMove, R_Min_UnderMove_DF$Cell_ID))
      XYCell_UnderMove <-  st_as_sf(XYCell_UnderMove, coords = c("x", "y")) %>% st_set_crs(Proj_Sys)
      
      # EXTRACT CHM UNDER FOR LOW LYING POINTS.
      CHM_Under_UnderXY <- raster::extract(R_CHM_Under, XYCell_UnderMove)
      
      Under_Move_DF <- data.frame(R_Min_UnderMove_DF, CHM_Under_UnderXY)
      colnames(Under_Move_DF) <- c("TreeID" ,"Cell_ID", "Min_UnderMove", "Max_Under")
      Gap_underTID <- as.data.frame(Under_Move_DF %>%
                                      dplyr::group_by(TreeID) %>%
                                      dplyr::summarise(Gap_Under <- min(Min_UnderMove-Max_Under), 
                                                       .groups = 'drop'))
      
      # IDENTIFY ALL LOW LYING TID THAT ARE WITHIN Para_EDT_V OF UNDERSTOREY CANOPY
      Move_TID <- Gap_underTID$TreeID[which(Gap_underTID$Gap_Under < Para_EDT_V)]
      LAS_Under_TID_VX_Move <- filter_poi(LAS_Vox_TID, TreeID %in% Move_TID)
      
      if(nrow(LAS_Under_TID_VX_Move@data) > 0){ 
        
        # USING MAX Z FOR LOW LYING TID (THOSE THAT WERE IDENTIFIED DUE TO GAPS IN CANOPY) 
        Max_Z_Moved <- max(LAS_Under_TID_VX_Move$Z)
        Max_Z_allTID <- as.data.frame(LAS_Vox_TID@data %>%
                                        dplyr::group_by(TreeID) %>%
                                        dplyr::summarise(Max_Z_Moved = max(Z), 
                                                         .groups = 'drop'))
        TID_Low <- Max_Z_allTID$TreeID[which(Max_Z_allTID$Max_Z_Moved < Max_Z_Moved)]
        
        Move_TID <- unique(c(Move_TID, TID_Low))
        LAS_Under_TID_VX_Move <- filter_poi(LAS_Vox_TID, TreeID %in% Move_TID)
        LAS_Under_TID_Move <- filter_poi(LAS_TID, TreeID %in% Move_TID)
        
        # UPDATING LAS_VOX_TID and LAS_TID
        LAS_Vox_TID <- filter_poi(LAS_Vox_TID, !(TreeID %in% Move_TID))
        LAS_TID <- filter_poi(LAS_TID, !(TreeID %in% Move_TID))
        
        # UPDATING VOX_ZERO (REQUIRES USING RASTER OF LOW LYING TID AND GETTING VOXELS BELOW RASTER)
        R_Under_TID_Vx_Move<-   grid_metrics(LAS_Under_TID_VX_Move, max(Z), res = Para_Vx_R_Res) 
        
        R_Under_TID_Move<-   grid_metrics(LAS_Under_TID_Move, max(Z), res = Para_Vx_R_Res) 
        
        # GET ALL ZERO VOXELS BELOW THE TID THAT WILL BE MOVED (Vox_TID)
        LAS_Vox_Zero <- merge_spatial(LAS_Vox_Zero, R_Under_TID_Vx_Move, attribute = "WSID")  
        LAS_Vox_Zero_Move <- filter_poi(LAS_Vox_Zero, (WSID-Z) >= 0)
        if(nrow(LAS_Vox_Zero_Move@data) > 0){ 
          LAS_Vox_Zero@data <- as.data.table(anti_join(LAS_Vox_Zero@data, LAS_Vox_Zero_Move@data, by=c("VoxID"))) 
          
          # GET ALL ZERO VOXELS BELOW THE TID THAT WILL BE MOVED (Vox_TID)
          LAS_Zero <- merge_spatial(LAS_Zero, R_Under_TID_Move, attribute = "WSID") 
          LAS_Zero_Move <- filter_poi(LAS_Zero, (WSID-Z) >= 0)
          LAS_Zero@data <- as.data.table(anti_join(LAS_Zero@data, LAS_Zero_Move@data, by=c("VoxID"))) 

          # REMOVING LAS ATTRIBUTES IN PREPERATION FOR MERGING VOXELS FROM ZERO AND VOXELS FROM VOX_TID THAT NEED MOVING INTO UNDERSTOREY
          Under_ColNames <- colnames(LAS_IDW_Under@data)
          Index_Col1 <-  which(colnames(LAS_Vox_Zero_Move@data) %in%   Under_ColNames  ) 
          LAS_Vox_Zero_Move@data <- LAS_Vox_Zero_Move@data[,..Index_Col1]
          
          Missing_Columns <- setdiff(colnames(LAS_IDW_Under@data), colnames(LAS_Vox_Zero_Move@data))
          for(MC in 1:length(Missing_Columns)){
            LAS_Vox_Zero_Move <- add_lasattribute(LAS_Vox_Zero_Move, x=as.numeric(0), name=Missing_Columns[MC], desc =Missing_Columns[MC])
          }
          
          Index_Col2 <-  which(colnames(LAS_Under_TID_VX_Move@data) %in%    Under_ColNames )          
          LAS_Under_TID_VX_Move@data <- LAS_Under_TID_VX_Move@data[,..Index_Col2]
          
          Missing_Columns <- setdiff(colnames(LAS_IDW_Under@data), colnames(LAS_Under_TID_VX_Move@data))
          for(MC in 1:length(Missing_Columns)){
            LAS_Under_TID_VX_Move <- add_lasattribute(LAS_Under_TID_VX_Move, x=as.numeric(0), name=Missing_Columns[MC], desc =Missing_Columns[MC])
          }
  
          LAS_Vox_Zero_Move@data <- rbind(LAS_Vox_Zero_Move@data, LAS_Under_TID_VX_Move@data)
          LAS_IDW_Under@data <- rbind(LAS_IDW_Under@data, LAS_Vox_Zero_Move@data)
        }
      }
    }
  }
  ########
  # BACKUP
  ########
  LAS_Zero_Backup10 <- LAS_Zero
  LAS_Vox_Zero_Backup10 <- LAS_Vox_Zero
  LAS_TID_Backup10 <- LAS_TID
  LAS_Vox_TID_Backup10 <- LAS_Vox_TID
  LAS_IDW_Under_Backup10 <- LAS_IDW_Under
  
  # CLEANING UP ALL TH LAS DATA !!!
  LAS_Zero <- LAS_Zero_Backup10
  LAS_Vox_Zero <- LAS_Vox_Zero_Backup10
  LAS_TID <- LAS_TID_Backup10
  LAS_Vox_TID <- LAS_Vox_TID_Backup10
  LAS_IDW_Under <- LAS_IDW_Under_Backup10
  print("Backup10")
  
  # OUTPUT LAS
  LAS_Zero@data$TreeID <- as.integer(LAS_Zero@data$TreeID)
  writeLAS(LAS_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Zero_Backup10.laz",sep=''))
  LAS_Vox_Zero@data$TreeID <- as.integer(LAS_Vox_Zero@data$TreeID)
  writeLAS(LAS_Vox_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_Zero_Backup10.laz",sep=''))
  LAS_TID@data$TreeID <- as.integer(LAS_TID@data$TreeID)
  writeLAS(LAS_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_TID_Backup10.laz",sep=''))
  LAS_Vox_TID@data$TreeID <- as.integer(LAS_Vox_TID@data$TreeID)
  writeLAS(LAS_Vox_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_TID_Backup10.laz",sep=''))
  LAS_IDW_Under@data$TreeID <- as.integer(LAS_IDW_Under@data$TreeID)
  writeLAS(LAS_IDW_Under, paste(Output_LAS,"/F", Flights[f] ,"_LAS_IDW_Under_Backup10.laz",sep=''))
  
  # TIMER   # 
  Time_Taken <- as.vector(round(((proc.time() - ptm)[3]/60),3))
  TIME_Section <- c(TIME_Section, Time_Taken)
  TIME_Section_ID <- c(TIME_Section_ID, 10)
  TIME__Count_Vox_TID <- c(TIME__Count_Vox_TID, length(unique(LAS_Vox_TID@data$TreeID)))
  TIME__Count_TID <- c(TIME__Count_TID, length(unique(LAS_TID@data$TreeID)))
  Output_Table_TIME_Section_ID <- data.frame(TIME_Section_ID, TIME_Section, TIME__Count_Vox_TID, TIME__Count_TID)
  write.csv(Output_Table_TIME_Section_ID, paste(Output_CSV, "/F", Flight_Number ,"_10_TIME_TAKEN.csv",sep=''), row.names = FALSE)
  
  ptm <- proc.time()

  ####################################################################################################################################### 13
  ####################################################################################################################################### 13
  # CHARACTERISING THE St (FOR N BA AND SA PREDICTIONS)
  ####################################################################################################################################### 13
  ####################################################################################################################################### 13
  
  LAS_TID <- add_lasattribute(LAS_TID, x=as.numeric(1), name="Vox_PntCnt", desc ="Vox_PntCnt")
  
  ptm <- proc.time()
  
  # SUMMARY LAS SO THAT I CAN IDENTIFY SURROUNDING St WITHIN SUBJECT St IN VIDEO
  Bound_all_TID_Range <- as.data.frame(LAS_TID@data %>%
                                         dplyr::group_by(TreeID) %>%
                                         dplyr::summarise(MinX = min(X)-Para_Vid_Surrt ,
                                                          MaxX = max(X)+Para_Vid_Surrt ,
                                                          MinY = min(Y)-Para_Vid_Surrt ,
                                                          MaxY = max(Y)+Para_Vid_Surrt ,
                                                          MinZ = min(Z)-Para_Vid_Surrt ,
                                                          MaxZ = max(Z)+Para_Vid_Surrt ,
                                                          Range = range(Z)[2]-range(Z)[1], 
                                                          .groups = 'drop'))
  # LOOP THROUGH EACH St
  Unique_TID <- Bound_all_TID_Range$TreeID[rev(order(Bound_all_TID_Range$Range))] # ORDERED
  for(NP in 1:length(Unique_TID)){
    
    LAS_One_TID <- filter_poi(LAS_TID, TreeID == Unique_TID[NP])
    LAS_One_TID_Plot <- LAS_One_TID
    if(nrow(LAS_One_TID@data) > 0){ 
      if(nrow(LAS_One_TID@data) > 100){

        ##########################
        # PROCESS STAND ATTRIBUTES
        ##########################
        Att_St_Sl <- LAS_One_TID@data[, GET_TREE_ATTRIBUTES_FUN(X, Y, Z, TreeID, Vox_PntCnt, Para_Sl_Z, Unique_Poly_ID = 0), by = "TreeID"] 
        
        # GETTING OUTPUT FOR TREE ATTRIBUTE FUNCTION (GET_TREE_ATTRIBUTES_FUN)
        Poly_WS_Slice <-unlist(Att_St_Sl$Poly_eachWS_All_SPDF)[[1]]
        Poly_CHull_Slice <-unlist(Att_St_Sl$Poly_eachSl_All_SPDF)[[1]]
        Centroid_XYZ <-unlist(Att_St_Sl$Points_Centroids_XYZ)[[1]]

        ########################
        # CANOPY BASE ASSESSMENT
        ########################
        
        Shift_X_Plot <- min(LAS_One_TID_Plot$X)
        Shift_Y_Plot <-min(LAS_One_TID_Plot$Y)
        
        LAS_One_TID_Plot$X <- LAS_One_TID_Plot$X - Shift_X_Plot
        LAS_One_TID_Plot$Y <- LAS_One_TID_Plot$Y - Shift_Y_Plot
        
        St_Profile <- GAP_DENSITY_FUN(LAS_One_TID_Plot$Z,
                                      Para_BW = Para_Canopy_Base_BW ,
                                      Para_Threshold_Percent = Para_Canopy_Base_Thresh ,
                                      Plot = "No",
                                      TreeID = Unique_TID[NP])
        
        # WORK AROUND FOR WHEN ONE OF THEM IS INF
        OneInf <- which(c(St_Profile$End_Largest_Gap, St_Profile$Start_Largest_Gap) != Inf)
        if(length(OneInf) != 2){
          St_Base_Height <- c(St_Profile$End_Largest_Gap, St_Profile$Start_Largest_Gap)[OneInf]
          C_Base_Height <- c(St_Profile$End_Largest_Gap, St_Profile$Start_Largest_Gap)[OneInf]
        }else{
          St_Base_Height <- min(St_Profile$End_Largest_Gap, St_Profile$Start_Largest_Gap)
          C_Base_Height <- max(St_Profile$End_Largest_Gap, St_Profile$Start_Largest_Gap)
        }
        
        C_Base_XYZ <-as.data.frame(Poly_WS_Slice[which.min(abs(Poly_WS_Slice$PolyHeight - C_Base_Height)), match(c("X_Cent", "Y_Cent", "PolyHeight"), names(Poly_WS_Slice))])
        St_Base_XYZ <- as.data.frame(Poly_WS_Slice[which.min(abs(Poly_WS_Slice$PolyHeight - St_Base_Height)), match(c("X_Cent", "Y_Cent", "PolyHeight"), names(Poly_WS_Slice))])
        
        # GET SUMMARY OF EACH TID AND ORDER WITH MINZ AT TOP # LAS_TID
        Summary_Trunk_Range <- as.data.frame(Poly_WS_Slice@data %>%
                                               dplyr::group_by(PolyHeight) %>%
                                               dplyr::summarise(WSGrid_Count_Max = max(WSGrid_Count),
                                                                WSGrid_Count_Total = sum(WSGrid_Count),
                                                                DiamMax = max(DiamMax),
                                                                DiamMin = max(DiamMin),
                                                                X_Cent =  mean(X_Cent),
                                                                Y_Cent =  mean(Y_Cent),
                                                                Count_Poly = length(WSGrid_Count), 
                                                                .groups = 'drop'))
        
        Within_St <- Summary_Trunk_Range[which(Summary_Trunk_Range$PolyHeight >= St_Base_XYZ$PolyHeight &
                                                 Summary_Trunk_Range$PolyHeight <= C_Base_XYZ$PolyHeight), ]
        
        Cutoff_St <-Within_St[c(1, nrow(Within_St)),]
        
        Below_StBase <- Summary_Trunk_Range[which(Summary_Trunk_Range$PolyHeight <= St_Base_XYZ$PolyHeight), ]
        Above_CBase <- Summary_Trunk_Range[which(Summary_Trunk_Range$PolyHeight >= C_Base_XYZ$PolyHeight), ]
        
        if(C_Base_Height != St_Base_Height){
          
          ##################
          # UPDATING St BASE
          ##################
          # CALCULATING BOTTOM OF St
          Below_Small <- which(Below_StBase$WSGrid_Count_Total < Cutoff_St$WSGrid_Count_Total[1]*Para_CutOff_Grid_Multiplier  & Below_StBase$Count_Poly <= Cutoff_St$Count_Poly [1])
          
          List_Rnge_ConsecRuns_Below <- rev(CONSECUTIVE_RUN_RANGE(Below_Small)%>% str_match_all("[0-9]+"))
          
          St_Base_Height <- Below_StBase$PolyHeight[as.numeric(List_Rnge_ConsecRuns_Below[[1]][1])]
          First_Consec_Range <- C_Base_Height - St_Base_Height
          
          if(length(List_Rnge_ConsecRuns_Below) > 1){
            for(SR in 2:length(List_Rnge_ConsecRuns_Below)){
              Below_Range <- as.numeric(List_Rnge_ConsecRuns_Below[[SR]])
              if(length(Below_Range) > 1){
                Range <- Below_StBase$PolyHeight[Below_Range[2]] - Below_StBase$PolyHeight[Below_Range[1]]
                if(Range > First_Consec_Range-Para_1st_Consec_Range_Reduce ){
                  # IF THE RANGE IS GREATER THAN THE FIRST CONSECUTIVE RANGE
                  St_Base_Height <- Below_StBase$PolyHeight[as.numeric(List_Rnge_ConsecRuns_Below[[SR]][1])]
                }
              }
            }
          }
          St_Base_XYZ$PolyHeight <- St_Base_Height
          St_Base_XYZ$X_Cent <- Below_StBase$X_Cent[which(Below_StBase$PolyHeight == St_Base_XYZ$PolyHeight)]
          St_Base_XYZ$Y_Cent <-Below_StBase$Y_Cent[which(Below_StBase$PolyHeight == St_Base_XYZ$PolyHeight)]
          
          ######################
          # UPDATING CANOPY BASE
          ######################
          
          # CALCULATING ABOVE St
          Above_Small <- which(Above_CBase$WSGrid_Count_Total < Cutoff_St$WSGrid_Count_Total[2]*Para_CutOff_Grid_Multiplier  & Above_CBase$Count_Poly <= Cutoff_St$Count_Poly [2])
          List_Rnge_ConsecRuns_Above <- CONSECUTIVE_RUN_RANGE(Above_Small)%>% str_match_all("[0-9]+")
          
          First_Consec_Range <- Above_CBase$PolyHeight[as.numeric(List_Rnge_ConsecRuns_Above[[1]])[length(as.numeric(List_Rnge_ConsecRuns_Above[[1]]))]] - St_Base_Height
          
          Upper_Range <- as.numeric(List_Rnge_ConsecRuns_Above[[1]])[length(as.numeric(List_Rnge_ConsecRuns_Above[[1]]))]
          
          C_Base_Height   <- Above_CBase$PolyHeight[max(Upper_Range)]
          
          if(length(List_Rnge_ConsecRuns_Above) > 1){
            
            for(SR in 2:length(List_Rnge_ConsecRuns_Above)){
              
              Upper_Range <- as.numeric(List_Rnge_ConsecRuns_Above[[SR]])
              if(length(Upper_Range) > 1){
                Range <- Above_CBase$PolyHeight[Upper_Range[2]] - Above_CBase$PolyHeight[Upper_Range[1]]
                
                if(Range > (First_Consec_Range-1)){
                  # IF THE RANGE IS GREATER THAN THE FIRST CONSECUTIVE RANGE
                  C_Base_Height <- Above_CBase$PolyHeight[as.numeric(List_Rnge_ConsecRuns_Above[[SR]][2])]
                }
              }
            }
          }
          
          C_Base_XYZ$PolyHeight <- C_Base_Height
          C_Base_XYZ$X_Cent <- Above_CBase$X_Cent[which(Above_CBase$PolyHeight == C_Base_XYZ$PolyHeight)]
          C_Base_XYZ$Y_Cent <-Above_CBase$Y_Cent[which(Above_CBase$PolyHeight == C_Base_XYZ$PolyHeight)]
          
          #######################################
          # REDUCE St LENGTH IF BRANCHES PRESENT
          #######################################
          
          # UPDATE St
          Within_St <- Summary_Trunk_Range[which(Summary_Trunk_Range$PolyHeight >= St_Base_Height &
                                                   Summary_Trunk_Range$PolyHeight <= C_Base_Height), ]
          Cutoff_St <-Within_St[c(1, nrow(Within_St)),]
          
          Within_Small <- which(Within_St$WSGrid_Count_Total < max(Cutoff_St$WSGrid_Count_Total)*Para_CutOff_Grid_Multiplier2  & Within_St$Count_Poly <= Cutoff_St$Count_Poly[1])
          
          # IF THERE ARE NO HEIGHTS WITH SMALL POLYGON GET THE MEDIAN POLYGON SIZE
          if(length(Within_Small) == 0){
            Within_St <- Summary_Trunk_Range[which(Summary_Trunk_Range$PolyHeight >= St_Base_Height &
                                                     Summary_Trunk_Range$PolyHeight <= C_Base_Height), ]
            Within_Small <- which(Within_St$WSGrid_Count_Total < median(Within_St$WSGrid_Count_Total) & Within_St$Count_Poly <= Cutoff_St$Count_Poly[1])
          }
          
          List_Range_ConsecRuns_WithinSt <- CONSECUTIVE_RUN_RANGE(Within_Small) %>% str_match_all("[0-9]+")
          Lowest_Range <- as.numeric(List_Range_ConsecRuns_WithinSt[[1]])
          
          # St BASE IS MIN OF THE LOWEST RUN OF CONSECTIVE POLYHEIGHTS
          St_Base_Height  <- Within_St$PolyHeight[min(Lowest_Range)]
          C_Base_Height  <- Within_St$PolyHeight[max(Lowest_Range)]
          First_Consec_Range <- C_Base_Height - St_Base_Height
          # IF THERE ARE OTHER CONSECUTIVE RUNS OF SMALL POLYGONS
          if(length(List_Range_ConsecRuns_WithinSt) > 1){
            for(SR in 2:length(List_Range_ConsecRuns_WithinSt)){
              Upper_Range <- as.numeric(List_Range_ConsecRuns_WithinSt[[SR]])
              if(length(Upper_Range) > 1){
                
                Range <- Within_St$PolyHeight[Upper_Range[2]] - Within_St$PolyHeight[Upper_Range[1]] 
                if(Range >= (First_Consec_Range-2)){ 
                  # IF THE RANGE IS GREATER THAN THE FIRST CONSECUTIVE RANGE - 1 m
                  C_Base_Height <- Within_St$PolyHeight[as.numeric(List_Range_ConsecRuns_WithinSt[[SR]][2])]
                }
              }
            }
          }
          
          St_Base_XYZ$PolyHeight <- St_Base_Height
          St_Base_XYZ$X_Cent <- Within_St$X_Cent[which(Within_St$PolyHeight == St_Base_Height)]
          St_Base_XYZ$Y_Cent <-Within_St$Y_Cent[which(Within_St$PolyHeight == St_Base_Height)]
          
          C_Base_XYZ$PolyHeight <- C_Base_Height
          C_Base_XYZ$X_Cent <- Within_St$X_Cent[which(Within_St$PolyHeight == C_Base_Height)]
          C_Base_XYZ$Y_Cent <-Within_St$Y_Cent[which(Within_St$PolyHeight == C_Base_Height)]
          
          # IF CANOPY BASE HEIGHT AND St BASE HEIGHT ARE THE SAME
        }else{
          Diam_Ls_2m <- which(Summary_Trunk_Range$DiamMax < Para_DiamMax  | Summary_Trunk_Range$Count_Poly == Para_Count_Poly)
          List_Rnge_Diam_Ls_2m <- CONSECUTIVE_RUN_RANGE(Diam_Ls_2m)%>% str_match_all("[0-9]+")
          
          if(length(List_Rnge_Diam_Ls_2m[[1]]) > 1){
            St_Base_Height <- Summary_Trunk_Range$PolyHeight[as.numeric(List_Rnge_Diam_Ls_2m[[1]][1])]
            C_Base_Height <- Summary_Trunk_Range$PolyHeight[as.numeric(List_Rnge_Diam_Ls_2m[[1]][2])]
          }
          else{
            Index_First_Range <- min(which(do.call(rbind, lapply(List_Rnge_Diam_Ls_2m, function(x) length(x))) == 2))
            
            if(Index_First_Range != Inf){
              St_Base_Height <- Summary_Trunk_Range$PolyHeight[as.numeric(List_Rnge_Diam_Ls_2m[[Index_First_Range]][1])]
              C_Base_Height <- Summary_Trunk_Range$PolyHeight[as.numeric(List_Rnge_Diam_Ls_2m[[Index_First_Range]][2])]
            }else{
              St_Base_Height <- median(Summary_Trunk_Range$PolyHeight)
              C_Base_Height <- median(Summary_Trunk_Range$PolyHeight)
            }
          }
          First_Consec_Range <- C_Base_Height - St_Base_Height
          # IF THERE ARE OTHER CONSECUTIVE RUNS OF SMALL POLYGONS
          if(length(List_Rnge_Diam_Ls_2m) > 1){
            for(SR in 2:length(List_Rnge_Diam_Ls_2m)){
              Upper_Range <- as.numeric(List_Rnge_Diam_Ls_2m[[SR]])
              if(length(Upper_Range) > 1){
                Range <-  Summary_Trunk_Range$PolyHeight[Upper_Range[2]] - Summary_Trunk_Range$PolyHeight[Upper_Range[1]] 
                Distance_Between_Ranges <- Summary_Trunk_Range$PolyHeight[Upper_Range[1]] - C_Base_Height
                if((Range > First_Consec_Range-Para_1st_Consec_Range_Reduce ) & Distance_Between_Ranges < min(c(Range,First_Consec_Range))/Para_Dist_Btw_Ranges ){
                  # IF THE RANGE IS GREATER THAN THE FIRST CONSECUTIVE RANGE
                  C_Base_Height <- Summary_Trunk_Range$PolyHeight[as.numeric(List_Rnge_Diam_Ls_2m[[SR]][2])]
                }
              }
            }
          }
          
          # UPDATE St BASE XY
          St_Base_XYZ$PolyHeight <- St_Base_Height
          Index_St_XY <- which(Summary_Trunk_Range$PolyHeight == St_Base_XYZ$PolyHeight)
          if(length(Index_St_XY) == 0){
            Index_St_XY <-which.min(abs(Summary_Trunk_Range$PolyHeight - St_Base_XYZ$PolyHeight))
          }
          St_Base_XYZ$X_Cent <- Summary_Trunk_Range$X_Cent[Index_St_XY]
          St_Base_XYZ$Y_Cent <-Summary_Trunk_Range$Y_Cent[Index_St_XY]
          
          # UPDATE CANOPY BASE XY
          C_Base_XYZ$PolyHeight <- C_Base_Height
          Index_C_XY <- which(Summary_Trunk_Range$PolyHeight == C_Base_XYZ$PolyHeight)
          if(length(Index_C_XY) == 0){
            Index_C_XY <-which.min(abs(Summary_Trunk_Range$PolyHeight - C_Base_XYZ$PolyHeight))
          }
          C_Base_XYZ$X_Cent <- Summary_Trunk_Range$X_Cent[Index_C_XY]
          C_Base_XYZ$Y_Cent <-Summary_Trunk_Range$Y_Cent[Index_C_XY]
          
        } # END ELSE LOOP (StBase and CBase same)

        ##################################
        # CANOPY VOLUME AND St POLY DIAMETER/AREA
        ##################################
        
        Canopy_Poly_Volume <- sum(Poly_WS_Slice$PolyVolume[which(Poly_WS_Slice$PolyHeight > C_Base_XYZ$PolyHeight)])
        Canopy_Pnt_Count <- length(which(LAS_One_TID$Z  >= C_Base_XYZ$PolyHeight))
        Canopy_Poly_Pnt_Den <- Canopy_Pnt_Count/Canopy_Poly_Volume
        
        Canopy_Grid_Volume <- sum(Poly_WS_Slice$WSGrid_Volume[which(Poly_WS_Slice$PolyHeight > C_Base_XYZ$PolyHeight)])
        Canopy_Grid_Pnt_Den <- Canopy_Pnt_Count/Canopy_Grid_Volume

        St_Height_Summary <- as.data.frame(Poly_WS_Slice@data %>%
                                             dplyr::group_by(PolyHeight) %>%
                                             dplyr::summarise(PolyArea = sum(PolyArea),
                                                              WSGrid_Area = sum(WSGrid_Count)*Para_Vx_R_Res *Para_Vx_R_Res ,
                                                              Count_St = length(PolyArea),
                                                              Mean_minDiam = median(DiamMin),
                                                              Mean_maxDiam = median(DiamMax),
                                                              PolyVolume = sum(PolyVolume),
                                                              WSGrid_Volume = sum(WSGrid_Volume)
                                                              , .groups = 'drop'))
        
        Index_St_Height_Range <- which(St_Height_Summary$PolyHeight <= C_Base_XYZ$PolyHeight & St_Height_Summary$PolyHeight >= St_Base_XYZ$PolyHeight)
        # IF THERE IS NO STEM THEN THIS BECOMES ZERO
        if(length(Index_St_Height_Range) == 0){
          St_PolyArea_Median <- 0
          St_StCount_Median <- 0
          St_minDiam_Median <- 0
          St_maxDiam_Median <- 0
        }else{
          St_PolyArea_Median <- median(St_Height_Summary$PolyArea[Index_St_Height_Range])
          St_StCount_Median <- median(St_Height_Summary$Count_St[Index_St_Height_Range])
          St_minDiam_Median <- median(St_Height_Summary$Mean_minDiam[Index_St_Height_Range])
          St_maxDiam_Median <-  median(St_Height_Summary$Mean_maxDiam[Index_St_Height_Range])
        }
        C_Max_Branches <- max(St_Height_Summary$Count_St)
        
        Canopy_max_WSGrid_Area <-  max(St_Height_Summary$WSGrid_Area)
        Canopy_max_Poly_Area <-  max(St_Height_Summary$PolyArea)
        
        #############################################
        # GENERATE PCA WITH BEARING AND ANGLE OF St
        #############################################
        Index_Trunk_Heights <- which(Poly_WS_Slice$PolyHeight >= St_Base_Height& Poly_WS_Slice$PolyHeight <= C_Base_Height )
        
        # IF THERE IS A TRUNK IDENTIFIED
        if(length(Index_Trunk_Heights) != 0){
          # IF TRUNK LENGTH IS SHORT
          if(length(Index_Trunk_Heights) < Para_Min_TrunkHeight ){
            Index_Trunk_Heights <- Index_Trunk_Heights[1:min(length(Poly_WS_Slice),length(Index_Trunk_Heights))]
          }
          Centroids_StBase_to_CBase <- Poly_WS_Slice[Index_Trunk_Heights,
                                                     match(c("X_Cent", "Y_Cent", "PolyHeight") ,names(Poly_WS_Slice))]
          if(nrow(Centroids_StBase_to_CBase) > 3){
            names(Centroids_StBase_to_CBase) <- c("X", "Y", "Z") 
            Centroids_StBase_to_CBase <- data.frame(TreeID = rep(Unique_TID[NP], nrow(Centroids_StBase_to_CBase)), Centroids_StBase_to_CBase)
            
            STEM_PCA <- STEM_PCA_FUN(Centroids_StBase_to_CBase, Para_Slice_Height = min(Centroids_StBase_to_CBase$Height), Up_Down = "DOWN") 
            endpts <- STEM_PCA$STEM_PCA[,2:4]
            
            Bearing_St <-  earth.bear(endpts[1,2], endpts[1,1], endpts[2,2], endpts[2,1])
            St_TopBotDist_inXY <- sqrt((endpts[2,1]- endpts[1,1])^2 + (endpts[2,2]- endpts[1,2])^2)
            SlopingSt_Length_fromGround <- sqrt(St_TopBotDist_inXY^2 + endpts[2,3]^2)
            Angel_St <- endpts[2,3]/(SlopingSt_Length_fromGround)*90 
            Min_Pnt_Height <- 0 
            plot_PCASeg <- 1
          }else{
            Bearing_St <- 0
            Angel_St <- 0
            plot_PCASeg <- 0
          }
        }else{
          Bearing_St <- 0
          Angel_St <- 0
          plot_PCASeg <- 0
        }
        # St PERCENTILE AND QUANTILES
        LAS_oneSt_C <- filter_poi(LAS_One_TID, Z >= C_Base_Height)
        
        C_Quantiles <- t(data.frame(quantile(LAS_oneSt_C$Z, seq(0.05, 0.95, 0.05)))) 
        row.names(C_Quantiles) <- ""
        colnames(C_Quantiles) <- paste("C_PTile_", seq(5, 95, 5), sep="")
        
        C_Z_SD <- sd(LAS_oneSt_C$Z, na.rm = FALSE)
        C_Z_mean <- mean(LAS_oneSt_C$Z, na.rm = FALSE)
        C_Z_median <- median(LAS_oneSt_C$Z, na.rm = FALSE)
        C_XY_SD <- mean(sd(LAS_oneSt_C$X, na.rm = FALSE), sd(LAS_oneSt_C$Y, na.rm = FALSE))
        C_Length <- max(LAS_One_TID$Z) - C_Base_Height
        Final_TID <- Unique_TID[NP]
        
        Range_z <- max(LAS_One_TID$Z) - min(LAS_One_TID$Z) 
        Hits <- length(LAS_One_TID$Z)
        min_z <- min(LAS_One_TID$Z)
        
        # IF THERE IS NO STEM THEN RatioTrunk_RangeZ_Diam BECOMES ZERO
        if(St_maxDiam_Median == 0){
          RatioTrunk_RangeZ_Diam <- 0
        }else{
          RatioTrunk_RangeZ_Diam <- (C_Base_Height-St_Base_Height)/St_maxDiam_Median
        }
        
        ################
        # OUTPUT RESULTS
        ################

        if(class(St_Base_XYZ) != "data.frame"){
          browser()
          St_Base_X <- data.frame(St_Base_XYZ)
          St_Base_X <- data.frame(St_Base_XYZ)}
        
        oneSt_Tree_Attributes <- c(  TID = Final_TID,
                                     Range_z = Range_z,
                                     Hits = Hits,
                                     min_z = min_z,
                                     RatioTrunk_RangeZ_Diam = RatioTrunk_RangeZ_Diam,
                                     St_Base_Height = St_Base_Height,
                                     St_Base_X = St_Base_XYZ[1,1],
                                     St_Base_Y = St_Base_XYZ[1,2],
                                     St_Height = max(LAS_One_TID$Z),
                                     St_TrunkCount =  St_StCount_Median,
                                     St_TrunkArea = St_PolyArea_Median,
                                     St_TrunkminDiam =St_minDiam_Median,
                                     St_TrunkmaxDiam = St_maxDiam_Median,
                                     Canopy_Base_Height = C_Base_Height,
                                     Canopy_MaxBranch = C_Max_Branches,
                                     Canopy_Poly_Volume = Canopy_Poly_Volume,
                                     Canopy_Grid_Volume = Canopy_Grid_Volume,
                                     Canopy_Pnt_Count = Canopy_Pnt_Count,
                                     Canopy_Poly_Pnt_Den = Canopy_Poly_Pnt_Den,
                                     Canopy_Grid_Pnt_Den = Canopy_Grid_Pnt_Den,
                                     Canopy_Height_maxPolyVol = max(St_Height_Summary$PolyHeight[which(St_Height_Summary$PolyVolume == max(St_Height_Summary$PolyVolume))]),
                                     Canopy_Volume_maxWSGridVol = max(St_Height_Summary$WSGrid_Volume[which(St_Height_Summary$WSGrid_Volume == max(St_Height_Summary$WSGrid_Volume))]),
                                     Canopy_max_WSGrid_Area = Canopy_max_WSGrid_Area,
                                     Canopy_max_Poly_Area = Canopy_max_Poly_Area,
                                     Bearing_St = Bearing_St,
                                     Angel_St = Angel_St,
                                     C_Quantiles = C_Quantiles,
                                     C_Z_SD = C_Z_SD,
                                     C_Z_mean = C_Z_mean,
                                     C_Z_median = C_Z_median,
                                     C_XY_SD = C_XY_SD,
                                     C_Length = C_Length
        )
        
        
        if(class(oneSt_Tree_Attributes) != "numeric"){
          browser()}
        
        ##################################
        # AGGREGATING DATA INTO PLOT LEVEL
        ##################################
        
        if(NP == 1){
          Poly_WS_Slice_All <- Poly_WS_Slice
          Poly_CHull_Slice_All <- Poly_CHull_Slice
          Centroid_XYZ_All <-  Centroid_XYZ
          AllSt_Tree_Attributes <- oneSt_Tree_Attributes
        }else{
          Poly_WS_Slice_All <- rbind(Poly_WS_Slice_All, Poly_WS_Slice)
          Poly_CHull_Slice_All <- rbind(Poly_CHull_Slice_All, Poly_CHull_Slice)
          Centroid_XYZ_All <-  rbind(Centroid_XYZ_All,Centroid_XYZ)
          AllSt_Tree_Attributes <- rbind(AllSt_Tree_Attributes, oneSt_Tree_Attributes)
        }
      }
    }
  } # NP LOOP!
  
  AllSt_Tree_Attributes_DF <- as.data.frame(AllSt_Tree_Attributes)
  write.csv(AllSt_Tree_Attributes_DF, paste(Output_CSV, "/F", Flight_Number ,"_St_Attributes.csv",sep=''), row.names = FALSE)
  

  #########################################################
  # SUBSET LAS FILE TO ONLY INCLUDED STEMS WITHIN 30 m PLOT
  #########################################################

  SUBSET_Half_Dim <- 30

  Flight_Centroid <- Flight_Centroids[which(as.numeric(Flight_Centroids$Flight) == Flight_Number),]
  Subset_Extent <- matrix(c(c((Flight_Centroid$POINT_X - SUBSET_Half_Dim),(Flight_Centroid$POINT_X + SUBSET_Half_Dim)), c((Flight_Centroid$POINT_Y - SUBSET_Half_Dim),(Flight_Centroid$POINT_Y + SUBSET_Half_Dim))), 2,2)

  AllSt_Tree_Attributes_Within_30Plot <- AllSt_Tree_Attributes_DF[which(AllSt_Tree_Attributes_DF$St_Base_X >= Subset_Extent[1,1]&
                                                                          AllSt_Tree_Attributes_DF$St_Base_X <=  Subset_Extent[2,1] &
                                                                          AllSt_Tree_Attributes_DF$St_Base_Y >= Subset_Extent[1,2]&
                                                                          AllSt_Tree_Attributes_DF$St_Base_Y <=  Subset_Extent[2,2]),]


  write.csv(AllSt_Tree_Attributes_Within_30Plot, paste(Output_CSV, "/F", Flight_Number ,"_St_Attributes_60X60.csv",sep=''), row.names = FALSE)

  # OUTPUTTING FINAL LAS (REMOVAL OF EDGE EFFECT BY ONLY SELECTING STEM WITHIN 30X30 M PLOT AFTER RUNNING OVER 40X40 m PLOT)
  LAS_TID_Within_30Plot <- filter_poi(LAS_TID, TreeID %in% AllSt_Tree_Attributes_Within_30Plot$TID)
  LAS_Vox_TID_Within_30Plot <- filter_poi(LAS_Vox_TID, TreeID %in% AllSt_Tree_Attributes_Within_30Plot$TID)
  LAS_Vox_Zero_Within_30Plot <- filter_poi(LAS_Vox_Zero, TreeID %in% AllSt_Tree_Attributes_Within_30Plot$TID)
  LAS_Zero_Within_30Plot <- filter_poi(LAS_Zero, TreeID %in% AllSt_Tree_Attributes_Within_30Plot$TID)

  LAS_TID_Within_30Plot@data$TreeID <- as.integer(LAS_TID_Within_30Plot@data$TreeID)
  writeLAS(LAS_TID_Within_30Plot, paste(Output_LAS,"/F", Flight_Number ,"_LAS_TID_Within_30Plot.laz",sep=''))
  LAS_Vox_TID_Within_30Plot@data$TreeID <- as.integer(LAS_Vox_TID_Within_30Plot@data$TreeID)
  writeLAS(LAS_Vox_TID_Within_30Plot, paste(Output_LAS,"/F", Flight_Number ,"_LAS_Vox_TID_Within_30Plot.laz",sep=''))

  if(nrow(LAS_Vox_Zero_Within_30Plot@data) > 0){ 
  LAS_Vox_Zero_Within_30Plot@data$TreeID <- as.integer(LAS_Vox_Zero_Within_30Plot@data$TreeID)
  writeLAS(LAS_Vox_Zero_Within_30Plot, paste(Output_LAS,"/F", Flight_Number ,"_LAS_Vox_Zero_Within_30Plot.laz",sep=''))
  }

  if(nrow(LAS_Zero_Within_30Plot@data) > 0){
    LAS_Zero_Within_30Plot@data$TreeID <- as.integer(LAS_Zero_Within_30Plot@data$TreeID)
    writeLAS(LAS_Zero_Within_30Plot, paste(Output_LAS,"/F", Flight_Number ,"_LAS_Zero_Within_30Plot.laz",sep=''))
    }


  ########################################################################################################################################## 11
  ########################################################################################################################################## 11
  # OUTPUTTING FINAL LAS FILES
  ########################################################################################################################################## 11
  ########################################################################################################################################## 11

  LAS_Zero@data$TreeID <- as.integer(LAS_Zero@data$TreeID)
  writeLAS(LAS_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Zero_BackupFINAL.laz",sep=''))
  LAS_Vox_Zero@data$TreeID <- as.integer(LAS_Vox_Zero@data$TreeID)
  writeLAS(LAS_Vox_Zero, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_Zero_BackupFINAL.laz",sep=''))

  LAS_TID@data$TreeID <- as.integer(LAS_TID@data$TreeID)
  writeLAS(LAS_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_TID_BackupFINAL.laz",sep=''))
  LAS_Vox_TID@data$TreeID <- as.integer(LAS_Vox_TID@data$TreeID)
  writeLAS(LAS_Vox_TID, paste(Output_LAS,"/F", Flights[f] ,"_LAS_Vox_TID_BackupFINAL.laz",sep=''))

  write.csv(AllSt_Tree_Attributes, paste(Output_CSV, "/F", Flights[f] ,"_St_Attributes.csv",sep=''), row.names = FALSE)

} #f LOOP (FLIGHT LOOP)





