
##########################################################################################################################################################################################
###FUNCTIONS IN THIS SCRIPT
##########################################################################################################################################################################################
# numextract
# numextract_all
# VOX_SLICE_DIST_FUN
# WS_CLUSTERING_FUN
# WS_MOVING_FUN
# SLICE_SURROUND_BBOX_FUN
# GAP_DENSITY_FUN
# VORONOI_FOCAL_FUN
# VOXEL_FUN
# CONVHULL_FUN
# GET_TREE_ATTRIBUTES_FUN
# STEM_PCA_FUN
# plot_UAS_function
# plot_function

##########################################################################################################################################
##########################################################################################################################################
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

numextract_all <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

##########################################################################################################################################
##########################################################################################################################################
# DISTANCE BETWEEN ALL VOXELS IN TWO FILES. MAX DISTANCE SPECIFIED
VOX_SLICE_DIST_FUN = function(XYZ_oneCl = XYZ_oneCl, 
                                         XYZ_Surr =XYZ_surrOneCl, 
                                         Max_Dist = Para_EDT_V, 
                                         Para_Vx_R_Res= Para_Vx_R_Res,
                                         Slice_Thickness= Para_Sl_Z){
  
  Slice_Seq <- seq(min(c(min(XYZ_oneCl$Z), min(XYZ_Surr$Z))), max(c(max(XYZ_oneCl$Z), max(XYZ_Surr$Z))) , Slice_Thickness)

  for(SS in 1:length(Slice_Seq)){
    
    # GET DATA FOR SLICE 
    Vx_oneCl_oneSl <- XYZ_oneCl [which(XYZ_oneCl$Z >= Slice_Seq[SS]-Slice_Thickness & XYZ_oneCl$Z <= Slice_Seq[SS]+2*Slice_Thickness ),]
    Vx_Surr_oneSl <- XYZ_Surr[which(XYZ_Surr$Z >= Slice_Seq[SS]-Slice_Thickness & XYZ_Surr$Z <= Slice_Seq[SS]+2*Slice_Thickness ),]
    
    # IF THERE IS SURROUND VOXELS
    if(nrow(Vx_Surr_oneSl) > 0){
      if(nrow(Vx_oneCl_oneSl) > 0){
        
        ######################  
        # DISTANCE CALCULATION
        ###################### 
        Index_Col1 <- match(c("X", "Y", "Z"), colnames(Vx_oneCl_oneSl))
        Index_Col2 <- match(c("X", "Y", "Z"), colnames(Vx_Surr_oneSl))
        Dist_oneCl_Surr <- rdist(as.data.frame(Vx_oneCl_oneSl)[,Index_Col1], 
                               as.data.frame(Vx_Surr_oneSl)[,Index_Col2]) # Row Col
        Dist_oneCl_Surr <- round( Dist_oneCl_Surr, 3)
        Index_Dist_oneCl_Surr <- as.data.frame(which(Dist_oneCl_Surr <= Max_Dist, arr.ind = TRUE))  
        
        # UPDATE TABLE OF ALL CLOSE MATCHES
        if(nrow(Index_Dist_oneCl_Surr) > 0){

          Closest_Vox_Surr<- Vx_Surr_oneSl[Index_Dist_oneCl_Surr$col,]
          colnames(Closest_Vox_Surr)[which(colnames(Closest_Vox_Surr) %in% c("TreeID", "VoxID"))] <- c("TID_Surr", "VoxID_Surr")
          
          Closest_Vox_St <- Vx_oneCl_oneSl[Index_Dist_oneCl_Surr$row,]
          colnames(Closest_Vox_St)[which(colnames(Closest_Vox_St) %in% c("TreeID", "VoxID"))] <- c("TID_oneCl", "VoxID_oneCl")
          Dist_St_Surr_Close <- Dist_oneCl_Surr[ cbind(Index_Dist_oneCl_Surr$row, Index_Dist_oneCl_Surr$col) ] 
          Closest_Vox_Surr$Dist <- Dist_St_Surr_Close
          Closest_Vox_Surr_oneCl <- data.frame(Closest_Vox_Surr, Closest_Vox_St)
          
          if(!exists("All_Closest_Vox_Surr_oneCl")){
            All_Closest_Vox_Surr_oneCl <- Closest_Vox_Surr_oneCl
          }else{
            All_Closest_Vox_Surr_oneCl<- rbind(All_Closest_Vox_Surr_oneCl, Closest_Vox_Surr_oneCl)
          }
        }
      }
    }
  }
  
  if(exists("All_Closest_Vox_Surr_oneCl")){
    # GENERATE SUMMARY OF THE CLOSEST SURROUNDING TID
    Index1 <- match(c("X", "Y", "Z"), colnames(All_Closest_Vox_Surr_oneCl))
    colnames(All_Closest_Vox_Surr_oneCl)[Index1] <- c("X_Surr", "Y_Surr", "Z_Surr")
    Index2 <- match(c("X.1", "Y.1", "Z.1"), colnames(All_Closest_Vox_Surr_oneCl))
    colnames(All_Closest_Vox_Surr_oneCl)[Index2] <- c("X_oneCl", "Y_oneCl", "Z_oneCl")

    All_Closest_Vox_Surr_oneCl <- All_Closest_Vox_Surr_oneCl[order(All_Closest_Vox_Surr_oneCl$VoxID_oneCl, All_Closest_Vox_Surr_oneCl$Dist), ] #sort by id and reverse of abs(value)
    All_Closest_Vox_Surr_oneCl_NoDup <-  All_Closest_Vox_Surr_oneCl[!duplicated(All_Closest_Vox_Surr_oneCl$VoxID_oneCl),] 

    Summary_closeSurr <- as.data.frame(All_Closest_Vox_Surr_oneCl_NoDup %>%
                                                 dplyr::group_by(TID_Surr) %>%
                                                 dplyr::summarise(minDist = min(Dist),
                                                                  cntVox = length(which(Dist <= sqrt(Para_Vx_R_Res^2 + Para_Vx_R_Res^2) + 0.01)),
                                                                  MinZ = min(Z_Surr),
                                                                  MaxZ = max(Z_Surr),
                                                                  .groups = 'drop'))
  }else{
    All_Closest_Vox_Surr_oneCl_NoDup <- NULL 
    Summary_closeSurr <- NULL
  }
  return(list(All_Closest_Vox_Surr_oneCl_NoDup = All_Closest_Vox_Surr_oneCl_NoDup,
              Summary_closeSurr = Summary_closeSurr))
}

##########################################################################################################################################
##########################################################################################################################################
# THIS FUNCTION COMPUTES THE CENTROIDS USING THE WS DERIVED FROM THE POINT DENSITY DISTRIBUTIONS OVER A SLICE. 
WS_CLUSTERING_FUN = function(LAS_WS_Cl,  Para_Vx_R_Res = 0.2, LAS_or_Vox = "Vox"){

  # RASTERISING THE LAS SLICE with R_Den, R_MaxZ, R_TID_MaxZ 
  if(LAS_or_Vox == "Vox"){
    R_Den <-   grid_metrics(LAS_WS_Cl, sum(Vox_PntCnt), res = Para_Vx_R_Res)
  }else{
    R_Den <-   grid_metrics(LAS_WS_Cl, length(Z), res = Para_Vx_R_Res)
  }
  R_MaxZ <-   grid_metrics(LAS_WS_Cl, max(Z), res = Para_Vx_R_Res)
  R_TID_MaxZ <-   grid_metrics(LAS_WS_Cl, TreeID[which.max(Z)], res = Para_Vx_R_Res)

  if(R_Den@ncols > 1 & R_Den@nrows > 1){
    
    # WATERSHED USING LAS DENSITY 
    algo = lidR::watershed(R_Den, tol = 1, ext=5)
    LAS_WS_Cl <- segment_trees(LAS_WS_Cl, algorithm = algo, attribute = "CentID", uniqueness = "incremental")
    R_WS <- algo()

    if(length(na.omit(R_WS@data@values)) > 0){

      Unique_CentID <- as.numeric(names(table(LAS_WS_Cl@data$CentID)))
      
      # WORK AROUND FOR WHEN THE ROWS AND COLS ARE BUGGY (FLIPPED FOR UNKNOWN REASON)
      if(!((R_WS@ncols == R_Den@ncols)&(R_WS@nrows == R_Den@nrows))){
        R_WS_Nrows <- R_WS@nrows
        R_WS_Ncols <- R_WS@ncols
        R_WS@ncols <-  R_WS_Nrows
        R_WS@nrows <-  R_WS_Ncols
        }

      # STACK RASTERS 
      R_Stack <- stack(R_Den,
                       R_WS,
                       R_MaxZ,
                       R_TID_MaxZ) # R_oneSl_ClID,
      names(R_Stack@layers[[1]]) <- "ZCount"
      names(R_Stack@layers[[2]]) <- "WS"
      names(R_Stack@layers[[3]]) <- "maxZ"
      names(R_Stack@layers[[4]]) <- "TID_MaxZ"
      
      # GENERATES DF OF RASTER VALUES (...CONVERT FACTOR TO NUMERIC.... REMOVES NA)
      R_Stack_DF <- data.frame(Cell_ID = seq(1, ncell(R_Den),1) , getValues(R_Stack))
      R_Stack_DF <- data.frame(sapply(R_Stack_DF, function(x) as.numeric(as.character(x))))
      R_Stack_DF <- na.omit(R_Stack_DF)
      
      ###################################################
      # GET LOCATION OF MAX DENSITY GRID CELL FOR EACH WS
      ###################################################

      # LOCATE MAX DENSITY CELL_ID FOR EACH WS
      Cell_ID_maxDen <- as.data.frame(R_Stack_DF %>%
                                        dplyr::select(Cell_ID, ZCount, WS) %>%
                                        dplyr::group_by(WS) %>%
                                        dplyr::summarise(Cell_ID_maxDen = Cell_ID[which.max(ZCount)][1], .groups = 'drop'))

      # CLEAN Cell_ID_maxDen
      Cell_ID_maxDen <- na.omit(Cell_ID_maxDen)
      Cell_ID_maxDen <- Cell_ID_maxDen[which(Cell_ID_maxDen$WS %in% Unique_CentID),]

      # XY OF EACH CELL (NEW CENTROID) REPRESENTING MAX DEN FOR EACH WS
      WS_XY_maxDen <- xyFromCell(R_WS, Cell_ID_maxDen$Cell_ID_maxDen, spatial=FALSE)
      colnames(WS_XY_maxDen) <- c("X", "Y")
      Max_Den <- raster::extract(R_Den, WS_XY_maxDen)

      CentID<-apply(X = data.frame(WS_XY_maxDen), MARGIN = 1, FUN = function(xy) R_WS@data@values[which.min(replace(distanceFromPoints(R_WS, xy), is.na(R_WS), NA))]) 
      WS_XY_maxDen <- data.frame(WS_XY_maxDen, Max_Den, CentID)

      # GET Z VALUE FOR NEW CENTROID (USING LIDAR WITH WS ID ASSOCIATED WITH CENTROID)
      WS_Median_Z_TID <- as.data.frame(LAS_WS_Cl@data %>%
                                          dplyr::group_by(CentID) %>%
                                          dplyr::summarise(Z = median(Z),
                                                           CentID = as.numeric(names(table(CentID)[which.max(table(CentID))]))[1], .groups = 'drop'))
      
      WS_Median_Z_TID <- na.omit(WS_Median_Z_TID)
      WS_XY_maxDen$Z <- WS_Median_Z_TID$Z[match(WS_XY_maxDen$CentID, WS_Median_Z_TID$CentID)] 

    return(list(WS_XY_maxDen = WS_XY_maxDen,
                R_Stack = R_Stack,
                LAS_WS_Cl = LAS_WS_Cl))
    }else{
      return(list(WS_XY_maxDen = NULL,
                  R_Stack = NULL,
                  LAS_WS_Cl = LAS_WS_Cl))
    }

  }else{
    return(list(WS_XY_maxDen = NULL,
                R_Stack = NULL,
                LAS_WS_Cl = LAS_WS_Cl))
  }

}

##########################################################################################################################################
##########################################################################################################################################
WS_MOVING_FUN = function(LAS_Above_Sl, LAS_Below_Sl, Para_Vx_R_Res = 0.2, Slice_Thickness = 0.4, Z_Inc_Slice_mm = Z_Inc_Slice[mm], LAS_or_Vox = "Vox"){
  
  #############################################
  # APPLY WS CLUSTERING FUNCTION ON ABOVE SLICE
  #############################################
  WS_Above_Cl <- WS_CLUSTERING_FUN(LAS_Above_Sl, Para_Vx_R_Res = Para_Vx_R_Res, LAS_or_Vox = "Vox")
  
  if(!is.null(WS_Above_Cl$WS_XY_maxDen)){
    
    WS_Class_Above <- WS_Above_Cl$R_Stack$WS
    WS_Class_Above_TID <- WS_Above_Cl$R_Stack$TID_MaxZ
    
    # DATA.FRAME OF GRID CELLS WITH WS CLUSTER ID
    XYCell_Above <- xyFromCell(WS_Class_Above, which(!is.na(WS_Class_Above@data@values)))
    XYCell_Above <- data.frame(XYCell_Above,
                               WS_ID=na.omit(WS_Class_Above@data@values),
                               TID = WS_Class_Above_TID@data@values[which(!is.na(WS_Class_Above@data@values))])
    Index_Col1 <-match(c("X", "Y"), colnames(LAS_Above_Sl@data))
    
    #############################################
    # APPLY WS CLUSTERING FUNCTION ON BELOW SLICE
    #############################################
    WS_Below_Cl <- WS_CLUSTERING_FUN(LAS_Below_Sl, Para_Vx_R_Res = Para_Vx_R_Res, LAS_or_Vox = "Vox")
    
    if(!is.null(WS_Below_Cl$WS_XY_maxDen)){
      
      Below_MaxDen_WS_XY<- WS_Below_Cl$WS_XY_maxDen
      
      WS_Class_Below_WS <- WS_Below_Cl$R_Stack$WS
      WS_Class_Below_GID  <- WS_Class_Below_WS
      
      XYCell_Below <- xyFromCell(WS_Class_Below_WS, which(!is.na(WS_Class_Below_WS@data@values)))
      XYCell_Below <- data.frame(XYCell_Below, WS_ID=na.omit(WS_Class_Below_WS@data@values))
      
      
      ##########################################################################
      # GIVE ALL LAS_BELOW_SL A GRI_ID (FOR NA USE rdist TO GET CLOSEST GRID ID)
      ##########################################################################
      
      # GET GRID ID OF ALL BELOW GRIDS WITH VALUE
      WS_Below_Values <- WS_Class_Below_WS$WS@data@values
      G_ID <- which(!is.na(WS_Class_Below_WS$WS@data@values))
      WS_Below_Values[!is.na(WS_Class_Below_WS$WS@data@values)] <- G_ID
      WS_Class_Below_GID <- setValues(WS_Class_Below_GID, WS_Below_Values)
      XYCell_Below <- data.frame(XYCell_Below,
                                 G_ID=WS_Class_Below_GID[which(!is.na(WS_Class_Below_WS@data@values))])
      #browser()
      LAS_Below_Sl <- merge_spatial(LAS_Below_Sl, WS_Class_Below_GID, attribute = "WSID") 
      index_Col1<-match(c("X", "Y"), colnames(LAS_Below_Sl@data))
      index_NA_GID <- which(is.na(LAS_Below_Sl@data$WSID))
      
      if(length(index_NA_GID)>0){
        Dist_BelowG_BeloNA <- rdist(XYCell_Below[, 1:2], LAS_Below_Sl@data[index_NA_GID, ..index_Col1]) 
        Nearest_G_Index <- apply(Dist_BelowG_BeloNA, 2, which.min)
        LAS_Below_Sl@data$WSID[index_NA_GID] <- G_ID[Nearest_G_Index]
      }
      
      ##########################################
      # GET TID (MOST COMMON) FOR EACH GRID CELL
      ##########################################
      Grid_Cell_TID <- as.data.frame(LAS_Below_Sl@data %>%
                                       dplyr::group_by(WSID) %>%
                                       dplyr::summarise(TID = TreeID[which.max(length(Z))], .groups = 'drop'))
      XYCell_Below <- XYCell_Below[which(XYCell_Below$G_ID %in% Grid_Cell_TID$WSID),]
      XYCell_Below$TID <- Grid_Cell_TID$TID
      
      ####################################
      # ID NEAREST BELOW WS FOR EACH VOXEL 
      ###################################
      
      # DISTANCE BETWEEN TOP AND BOTTOM GRID CELLS
      Dist_WS_Below_Above  <- rdist(XYCell_Below[,1:2],
                                    XYCell_Above[,1:2]) 
      Dist_WS_Below_Above  <- round(Dist_WS_Below_Above, 3)
      Dist_WS_Below_Above <- ifelse(Dist_WS_Below_Above > Slice_Thickness*2,NA,Dist_WS_Below_Above) # REMOVE PAIRS THAT ARE VERY FAR APART
      
      # NEAREST BELOW WS FOR EACH VOXEL
      if(nrow(Dist_WS_Below_Above) ==1) {
        Dist_WS_Below_Above<- as.vector(Dist_WS_Below_Above)
        
        NearestG_Below_To_Above <- which.min(Dist_WS_Below_Above)
        if(length(NearestG_Below_To_Above) == 0){
          NearestG_Below_To_Above <- NA
          NearestDist_Below_To_Above <- NA
        }else{
          NearestDist_Below_To_Above <- min(Dist_WS_Below_Above, na.rm=T)
        }
      }else{
        Dist_WS_Below_Above  <-data.frame(Dist_WS_Below_Above)
        NearestG_Below_To_Above <-apply(Dist_WS_Below_Above, 2, function(x) {if (all(is.na(x))) {NA}  else {which.min(x)} })
        NearestDist_Below_To_Above <- apply(Dist_WS_Below_Above, 2, function(x) {if (all(is.na(x))) {NA}  else {min(x, na.rm=T)} })
        
      }
      Index_Nearest <- which(NearestDist_Below_To_Above < Slice_Thickness*2 ) #Para_EDT_V == Slice_Thickness*2
      
      ################################################################
      # UPDATE ABOVE TID (USING NEAREST BELOW TID FOR EACH ABOVE GRID)
      ################################################################
      
      # UPDATE DATAFRAME OF GRID CELLS
      XYCell_Above$TID[Index_Nearest] <- XYCell_Below$TID[NearestG_Below_To_Above[Index_Nearest]]
      
      # UPDATE ABOVE RASTER
      WS_Above_Values <- WS_Class_Above$WS@data@values
      WS_Above_Values[which(!is.na(WS_Above_Values))] <- XYCell_Above$TID
      
      WS_Class_Above@data@names <- "TID"
      WS_Class_Above <- setValues(WS_Class_Above, WS_Above_Values)
      
      # UPDATE ABOVE LAS
      LAS_Above_Sl <- merge_spatial(LAS_Above_Sl, WS_Class_Above, attribute = "WSID") 
      LAS_Above_Sl@data$TreeID <- as.integer(LAS_Above_Sl@data$WSID) 
      
      
      # UPDATE LAS THAT HAD NA (FROM WS algorithm) USING rdist
      index_Col1<-match(c("X", "Y"), colnames(LAS_Above_Sl@data))
      index_NA_NewTID <- which(is.na(LAS_Above_Sl@data$WSID))
      
      if(length(index_NA_NewTID)>0){
        Dist_AboveG_AboveNA <- rdist(XYCell_Above[, 1:2], LAS_Above_Sl@data[index_NA_NewTID, ..index_Col1])
        Dist_AboveG_AboveNA <- round(Dist_AboveG_AboveNA, 2)
        Nearest_G_Index <- apply(Dist_AboveG_AboveNA, 2, which.min)
        if(length(Nearest_G_Index)> 0){
          LAS_Above_Sl@data$WSID[index_NA_NewTID] <- XYCell_Above$TID[Nearest_G_Index]
          LAS_Above_Sl@data$TreeID <- as.integer(LAS_Above_Sl@data$WSID) # ABOVE WSID AND TreeID IS SAME
        }
      }
      
      # UPDATE LAS TID (ONLY ONE SLICE NOT THE WHOLE 4)
      LAS_Above_oneSl_New <- filter_poi(LAS_Above_Sl, Z >= Z_Inc_Slice_mm &
                                         Z < Z_Inc_Slice_mm+ Slice_Thickness &
                                         WSID > 0)
      # IF WS BELOW ISN'T NULL
    }else{
      LAS_Above_oneSl_New <- filter_poi(LAS_Above_Sl, TreeID < -10000) # GENREATES EMPTY LAS 
    }
    # IF WS ABOVE ISN'T NULL
  }else{
    LAS_Above_oneSl_New <- filter_poi(LAS_Above_Sl, TreeID < -10000) # GENREATES EMPTY LAS 
  }
  return(list(LAS_Above_oneSl_New = LAS_Above_oneSl_New))
}

##########################################################################################################################################
##########################################################################################################################################
# SLICE_SURROUND_BBOX_FUN:
  # GETS ALL THE TID THAT FALL WITHIN THE BBoX Z RANGE
  # INCREMENTALLY MOVES UP THE SLICE AND PERFORMS WS ON ABOVE SLICE TO CHANGE ITS TID BASED ON PROXIMITY TO BELOW SLICE
  # THE WHILE LOOP USES THE NEW ABOVE SLICE VOXEL TID  TO REPEATIDLY UNDERTAKING DISTANCE CALCULATION BETWEEN ORIGINAL AND NEW TID VALUES.
    # THE LOOP STOPS WHEN NO FURTHER VOXELS HAVE BEEN CHANGED DUE TO CONSTRAINTS IN VOX_SLICE_DIST_FUN
SLICE_SURROUND_BBOX_FUN = function(LAS_surrOneCl_rmOneCl_allTID, LAS_Empty, Unique_TID,
                                   Para_EDT_V, Z_Inc_Slice, Para_Vx_R_Res, Para_Sl_Z, 
                                   LAS_or_Vox = "Vox"){
  
  LAS_bboxTID_allSl <- LAS_Empty
  Move_Above_Sl_NonGnd_All <- "No"
  for(mm in 1:length(Z_Inc_Slice)){

    # GET BELOW AND ABOVE SLICE
    LAS_Above_oneSl <- filter_poi(LAS_surrOneCl_rmOneCl_allTID, Z >= Z_Inc_Slice[mm] & Z < Z_Inc_Slice[mm]+(Para_EDT_V))
    LAS_Below_oneSl <- filter_poi(LAS_surrOneCl_rmOneCl_allTID, Z >= (Z_Inc_Slice[mm]-Para_EDT_V) & Z < Z_Inc_Slice[mm])

    if(nrow(LAS_Above_oneSl@data) > 0 & nrow(LAS_Below_oneSl@data) > 0){ 
      
      #####################################################
      # WS CLUSTER:  FOR MOVING ABOVE SLICE and BELOW SLICE 
      #####################################################
      
      WS_Move_Slice<- WS_MOVING_FUN(LAS_Above_oneSl, LAS_Below_oneSl,
                                    Para_Vx_R_Res = Para_Vx_R_Res,
                                    Slice_Thickness = Para_EDT_V,
                                    Z_Inc_Slice_mm = Z_Inc_Slice[mm],
                                    LAS_or_Vox = "Vox")
      LAS_Above_oneSl_New <- WS_Move_Slice$LAS_Above_oneSl_New
      
      if(nrow(LAS_Above_oneSl_New@data) > 0){ 
        
        # UPDATE WS_MOVING_FUN RESULTS  (ASSIGN TreeID AS WSID AND REMOVE WSID ATTRIBUTE)    
        LAS_Above_oneSl_New@data$TreeID <- as.integer(LAS_Above_oneSl_New@data$WSID) 
        Index_Remove <- which(colnames(LAS_Above_oneSl_New@data) == "WSID")
        LAS_Above_oneSl_New@data <- LAS_Above_oneSl_New@data[,-Index_Remove, with=FALSE]
        
        Index_Above_Sl_OrigTID <- which(unique(LAS_Above_oneSl_New@data$TreeID) %in% Unique_TID)
        if(length(Index_Above_Sl_OrigTID) > 0){
          
          # WHILE LOOP ONLY UPDATES VOXELS WITHIN Para_EDT_V AND REPEATS UNTIL ALL UPDATED WITHIN SLICE
          Check_Orig_Count <- 1
          Change_Orig_Count <- 0
          while(Check_Orig_Count != 0){

            # USE VOX PROXIMITY TO REPLACE ANY NonGnd TID with Neighbours
            LAS_Above_oneSl_New_Orig <-  filter_poi(LAS_Above_oneSl_New, TreeID %in% Unique_TID)
            LAS_Above_oneSl_New_Changed <-  filter_poi(LAS_Above_oneSl_New, !(TreeID %in% Unique_TID)) # CHANGED DUE TO "WS_MOVING_FUN"
 
            if(nrow(LAS_Above_oneSl_New_Changed@data) > 0){ 
              # GET CLOSE VOXEL SUMMARY BETWEEN NonGnd_TID and Gnd_TID
              Index_Col1 <- match(c("X", "Y", "Z", "TreeID", "VoxID"), colnames(LAS_Above_oneSl_New_Orig@data))
              XYZ_New_Orig <- LAS_Above_oneSl_New_Orig@data[,..Index_Col1]
              Index_Col2 <- match(c("X", "Y", "Z", "TreeID", "VoxID"), colnames(LAS_Above_oneSl_New_Changed@data))
              XYZ_New_Changed <- LAS_Above_oneSl_New_Changed@data[,..Index_Col2]
              
              ################################################################################################ 12
              # DIST CALC FOR ABOVE SLICE (BETWEEN VOXELS WITH ORIG TID AND THOSE THAT CHANGED  WITH "WS_MOVING_FUN")
              ################################################################################################ 12
              
              Dist_NewOrig_NewChange <- VOX_SLICE_DIST_FUN(XYZ_New_Orig,
                                                           XYZ_New_Changed, # THESE ARE LIKE THE SURROUNDING TID
                                                           Max_Dist = Para_EDT_V, 
                                                           Para_Vx_R_Res= Para_Vx_R_Res,
                                                           Slice_Thickness= Para_Sl_Z)
              
              Table_CloseVox <- Dist_NewOrig_NewChange$All_Closest_Vox_Surr_oneCl_NoDup
              
              # UPDATING TID IN ABOVE SLICE
              Index <- match(Table_CloseVox$VoxID_oneTID, LAS_Above_oneSl_New@data$VoxID )
              LAS_Above_oneSl_New@data$TreeID[Index] <- Table_CloseVox$TID_Surr
    
              # CHECKING IF WHILE LOOP SHOULD BE REPEATED
              Check_Orig_Count <- length(which(LAS_Above_oneSl_New@data$TreeID %in% Unique_TID))
              if(Change_Orig_Count == Check_Orig_Count){
                Check_Orig_Count <- 0 # BREAK OUT OF WHILE LOOP
              }
            }else{
              Check_Orig_Count <- 0 # BREAK OUT OF WHILE LOOP
            }
            Change_Orig_Count <- Check_Orig_Count
          } # WHILE LOOP
          
        } # IF ABOVE SLICE STILL HAS NonGnd TID

        if(nrow(LAS_Above_oneSl_New@data) > 0){ 
          # UPDATING LAS WITH ABOVE SLICE TID
          Index_1 <- which(LAS_surrOneCl_rmOneCl_allTID@data$VoxID %in%  LAS_Above_oneSl_New@data$VoxID )
          Index_2 <- match(LAS_surrOneCl_rmOneCl_allTID@data$VoxID[Index_1],  LAS_Above_oneSl_New@data$VoxID )
          LAS_surrOneCl_rmOneCl_allTID@data$TreeID[Index_1] <- as.integer(LAS_Above_oneSl_New@data$TreeID[Index_2])

          if(length(which(LAS_Above_oneSl_New@data$TreeID == Unique_TID)) > 0) {
            LAS_Above_oneSl_NonGnd <- filter_poi(LAS_Above_oneSl_New, TreeID %in% Unique_TID)
            if(nrow(LAS_bboxTID_allSl@data) == 0){
              LAS_bboxTID_allSl <- filter_poi(LAS_surrOneCl_rmOneCl_allTID, Z < Z_Inc_Slice[mm]+(Para_EDT_V) & TreeID == Unique_TID)
            }else{
              LAS_bboxTID_allSl@data <- rbind(LAS_bboxTID_allSl@data, LAS_Above_oneSl_NonGnd@data)
            }
            Move_Above_Sl_NonGnd_All <- "No"
          }else{
            Move_Above_Sl_NonGnd_All <- "Yes"
          }

          if(nrow(LAS_bboxTID_allSl@data) > 0 & Move_Above_Sl_NonGnd_All == "Yes"){ 

            # GET CLOSE VOXEL SUMMARY BETWEEN NonGnd_TID and Gnd_TID
            Index_Col1 <- match(c("X", "Y", "Z", "TreeID", "VoxID"), colnames(LAS_bboxTID_allSl@data))
            XYZ_Above_Sl_NonGnd_All <- LAS_bboxTID_allSl@data[,..Index_Col1]
            Index_Col2 <- match(c("X", "Y", "Z", "TreeID", "VoxID"), colnames(LAS_Above_oneSl_New@data))
            XYZ_Above_Sl_New <- LAS_Above_oneSl_New@data[,..Index_Col2]
            
            ################################################################################################
            # DIST CALC FOR ABOVE SLICE (BETWEEN VOXELS WITH ORIG TID AND THOSE THAT HAVE NOT BEEN CHANGED)
            ################################################################################################ 
            Dist_Above_Sl_NonGnd_All <- VOX_SLICE_DIST_FUN(XYZ_Above_Sl_NonGnd_All,
                                                           XYZ_Above_Sl_New,
                                                           Max_Dist = 40, # Para_EDT_H,
                                                           Para_Vx_R_Res= Para_Vx_R_Res,
                                                           Slice_Thickness= Para_Sl_Z)

            Dist_Above_Sl_NonGnd_All_Summary <- Dist_Above_Sl_NonGnd_All$Summary_closeSurr
            if(!is.null(Dist_Above_Sl_NonGnd_All_Summary)){
              if(nrow(Dist_Above_Sl_NonGnd_All_Summary) == 1){
                LAS_bboxTID_allSl@data$TreeID <- as.integer(Dist_Above_Sl_NonGnd_All_Summary$TID_Surr)
              }else{
                #browser()
                Dist_Above_Sl_NonGnd_All_VoxMatch <- Dist_Above_Sl_NonGnd_All$All_Closest_Vox_Surr_oneCl_NoDup
                
                Index_1 <- which(LAS_bboxTID_allSl@data$VoxID %in%  Dist_Above_Sl_NonGnd_All_VoxMatch$VoxID_oneCl)
                Index_2 <- match(LAS_bboxTID_allSl@data$VoxID[Index_1],  Dist_Above_Sl_NonGnd_All_VoxMatch$VoxID_oneCl)
                LAS_bboxTID_allSl@data$TreeID[Index_1] <- as.integer(Dist_Above_Sl_NonGnd_All_VoxMatch$TID_Surr[Index_2])
              }
              Index_1 <- which(LAS_surrOneCl_rmOneCl_allTID@data$VoxID %in%  LAS_bboxTID_allSl@data$VoxID )
              Index_2 <- match(LAS_surrOneCl_rmOneCl_allTID@data$VoxID[Index_1],  LAS_bboxTID_allSl@data$VoxID )
              LAS_surrOneCl_rmOneCl_allTID@data$TreeID[Index_1] <- as.integer(LAS_bboxTID_allSl@data$TreeID[Index_2])
            }
            LAS_bboxTID_allSl <- LAS_Empty
            Move_Above_Sl_NonGnd_All <- "No"
          }
        }
      }
    } # IF THERE ARE VOXELS ABOVE SLICE AND VOXELS BELOW SLICE
  } # END mm Slice LOOP
  return(list(LAS_surrOneCl_rmOneCl_allTID = LAS_surrOneCl_rmOneCl_allTID))
}

##########################################################################################################################################
##########################################################################################################################################
GAP_DENSITY_FUN = function(Z_Values,
                           Para_BW = 0.4,
                           Para_Threshold_Percent = 0.2,
                           Plot = "No",
                           TreeID = 1) {

  pdens <- density(Z_Values, bw=Para_BW)

  # FINDING PEAKS WHEN THERE ARE NO ZEROS
  PEAK_DIP_FUN = function(pdens)
  {
    TurnPnts <- turnpoints(pdens$y)
    Index <- c(1, TurnPnts$tppos, TurnPnts$n)
    Z <- c(0, pdens$x [TurnPnts$tppos], max(pdens$x)) # Adding one extra turning point at start and end
    Density <- c(0, pdens$y[TurnPnts$tppos], 0)    # Adding one extra turning point at start and end
    Peak_Dip_Summary <- t(data.frame(rbind(Z, Density, Index)))
    colnames(Peak_Dip_Summary) <- c("Z", "Density", "Index")
    Peak_Dip_Summary<- data.frame( Peak_Dip_Summary, Peak_Dip = rep("Dip", nrow(Peak_Dip_Summary)), stringsAsFactors = FALSE)
    Peak_Dip_Summary$Peak_Dip [c(1, nrow(Peak_Dip_Summary))] <- c("Start_End", "Start_End")

    #################
    # PEAKS AND DIPS
    #################

    Y_Peak <- TurnPnts$points[TurnPnts$peaks]
    X_Peak <- pdens$x[TurnPnts$pos[TurnPnts$peaks]]

    Peak_Dip_Summary$Peak_Dip[match(X_Peak,Peak_Dip_Summary$Z)] <- "Peak"

    Dip_DF <- Peak_Dip_Summary[which(Peak_Dip_Summary$Peak_Dip == "Dip"),]
    minDip_DF <- Dip_DF[which.min(Dip_DF$Density),]

    Peak_DF <- Peak_Dip_Summary[which(Peak_Dip_Summary$Peak_Dip == "Peak"),]
    maxPeak_DF <- Peak_DF[which.max(Peak_DF$Density),]

    Range_Den_minDip_maxPeak <- maxPeak_DF$Density - minDip_DF$Density

    return(list(Peak_Dip_Summary=Peak_Dip_Summary,
                minDip_DF = minDip_DF,
                maxPeak_DF = maxPeak_DF,
                Range_Den_minDip_maxPeak = Range_Den_minDip_maxPeak))
  }

  # GET Peak_Dip_Summary
  Peak_Dip_Summary <- PEAK_DIP_FUN(pdens)

  # GET ALL SECTIONS THAT ARE WITHIN (Para_Threshold_Percent*100)% ABOVE THE MINIMUM DIP
  Den_Threshold <- Peak_Dip_Summary$minDip_DF$Density + Peak_Dip_Summary$Range_Den_minDip_maxPeak * Para_Threshold_Percent

  # IF THERE ARE NO MINIMUM DIPS
  if(length(Den_Threshold) == 0){
    Start_Largest_Gap <- max(Peak_Dip_Summary$Peak_Dip_Summary$Z)
    End_Largest_Gap <- max(Peak_Dip_Summary$Peak_Dip_Summary$Z)

  }else{

    Near_min_density <- which(pdens$y < Den_Threshold)
    
    # REMOVE STRING OF LOW DENSITIES NEAR BOTTOM HEIGHT AND TOP CANOPY HEIGHT
    Index_PeakDip_inMinDen <- Peak_Dip_Summary$Peak_Dip_Summary$Index[which(Peak_Dip_Summary$Peak_Dip_Summary$Index %in% Near_min_density)]

    # IDENTIFY THE SECTION THAT HAS THE LARGEST STRETCH AND USE THOSE BOUNDS TO DETERMINE
    # UNDER TOP AND CANOPY BASE AS FIRST PASS

    # VEG_PROGILE_GAPS _WITH_LOW_LIDAR_DENSITY FUNCTION
    NEAR_MIN_DENSITY_FUN = function(Near_min_density)
    {
      # CALCULATING WHERE THE GAP IS
      Zero_Density <- rle(diff(Near_min_density))
      myZero_Density <- which(Zero_Density$values == TRUE & Zero_Density$lengths > 0)
      Zero_Density.lengths.cumsum <- cumsum(Zero_Density$lengths)
      ends <- Zero_Density.lengths.cumsum[myZero_Density]
      newindex <- ifelse(myZero_Density>1, myZero_Density-1, 0)
      starts <- Zero_Density.lengths.cumsum[newindex] + 1
      starts_2 <- starts
      if(length(which(newindex == 0))>0){starts_2 <-  c(1,starts)}

      starts <- starts_2
      Start_Height <- pdens$x[Near_min_density[starts]]
      End_Height <- pdens$x[Near_min_density[ends]]

      # REMOVING THE GAP BELOW THE GROUND
      if(End_Height[1] < 0){
        Start_Height <- Start_Height[-1]
        End_Height <- End_Height[-1]
      }

      return(list(Start_Height=Start_Height,
                  End_Height=End_Height))
    } # END NEAR_MIN_DENSITY_FUN
    
    # GET OUTPUT
    Near_Zero_Density_List <-NEAR_MIN_DENSITY_FUN(Near_min_density)
    Start_Height <- Near_Zero_Density_List$Start_Height
    End_Height <- Near_Zero_Density_List$End_Height

    # GET LARGEST GAP WITH DENSITY CLOSE TO ZERO AND DETERMINE THE UNDER AND CANOPY BASE
    Gaps <- End_Height -Start_Height
    Largest_Gap <- max(Gaps)
    Start_Largest_Gap <- Start_Height[which(Gaps == Largest_Gap)]
    End_Largest_Gap <- End_Height[which(Gaps == Largest_Gap)]

    # GET ALL PEAKS AND DIPS WITHIN FIRST PASS OF START AND END
    PeakDip_Within_minDen <- Peak_Dip_Summary$Peak_Dip_Summary[which(Peak_Dip_Summary$Peak_Dip_Summary$Index %in% Index_PeakDip_inMinDen),]

    # CALCULATE START OF GAP AS "FIRST MIN" WITHIN GAP
    Density_Dips_Start <- PeakDip_Within_minDen$Density[ which(PeakDip_Within_minDen$Z > Start_Largest_Gap &
                                                                 PeakDip_Within_minDen$Peak_Dip == "Dip")]
    Start_Largest_Gap <- min(PeakDip_Within_minDen$Z[ which(PeakDip_Within_minDen$Z > Start_Largest_Gap &
                                                              PeakDip_Within_minDen$Peak_Dip == "Dip" &
                                                              PeakDip_Within_minDen$Density <= median(Density_Dips_Start))])

    # CALCULATE END OF GAP AS "LAST MIN" WITHIN  GAP
    Density_Dips_Ends <- PeakDip_Within_minDen$Density[ which(PeakDip_Within_minDen$Z < End_Largest_Gap &
                                                                PeakDip_Within_minDen$Peak_Dip == "Dip")]
    End_Largest_Gap <- max(PeakDip_Within_minDen$Z[which(PeakDip_Within_minDen$Z < End_Largest_Gap &
                                                           PeakDip_Within_minDen$Peak_Dip == "Dip" &
                                                           PeakDip_Within_minDen$Density <= median(Density_Dips_Ends))])

  }

  if(Plot == "Yes"){

    plot(pdens$x, pdens$y, type="l", ylim=c(0, 0.6),
         main= paste( "F:", TreeID),
         ylab = "Density", 
         xlab = "Height (m)",
         cex.lab = 1.5,
         cex.axis = 2)

    abline(v = Start_Largest_Gap, col="blue", lwd = 3)
    abline(v = End_Largest_Gap, col="dark green", lwd = 3)
    abline(h = Den_Threshold, col= "red", lwd = 2)

    legend("topright", legend = c("Threshold Height", "Bottom Minima", "Top Minima"),
           lty= 1, col=c("red", "blue", "dark green"), lwd=3,
           cex= 1.7)
  }

  return(list(Start_Largest_Gap=Start_Largest_Gap,
              End_Largest_Gap=End_Largest_Gap,
              Peak_Dip_Summary = Peak_Dip_Summary))
}

##########################################################################################################################################
##########################################################################################################################################
# LOCAL MAX Focal function With Voronoi Assessment around maximas
VORONOI_FOCAL_FUN = function(Raster_Layer, Win_row=15, Win_col=15, min_or_max)
{
  # MAX VALUE WITHIN SEARCH WINDOW
  if(min_or_max == "max"){
    f <- function(X) max(X, na.rm=TRUE)
  }else{
    f <- function(X) min(X, na.rm=TRUE)
  }

  # MAX HEIGHT (CHM) IN SEARCH WINDOW
  ww <- matrix(1, nrow=Win_row, ncol=Win_col) 
  CHM_localmax <- focal(Raster_Layer, fun=f, w=ww, pad=TRUE, padValue=NA)

  # LOCAL MAXIMAS (CELL XY COORDS) 
  Local_Max <- Raster_Layer==CHM_localmax
  maxXY <- xyFromCell(Local_Max, Which(Local_Max==1, cells=TRUE))
  maxZvalue <- raster::extract(Raster_Layer, maxXY)

  # GETTING VORONOI OF EACH MAXIMA AND GETTING MEDIAN HEIGHT WITHIN VORONOI 
  maxXY_ID <- data.frame(maxXY, ID = seq(1, nrow(maxXY), 1))
  coordinates(maxXY_ID) = ~x+y
  # browser()
  Raster_Layer_BoundPoly <- as(extent(Raster_Layer), 'SpatialPolygons') 
  Maxima_Vori_Poly <- voronoi.polygons(maxXY_ID, Raster_Layer_BoundPoly)

  Median_Height_MaximLoc_VoriPoly <- raster::extract(Raster_Layer, Maxima_Vori_Poly,
                                                     fun =  median,
                                                     na.rm = TRUE)
  
  # SMALL WORK AROUND FOR CORNER CASE WITH VALUE Inf ...
  Gap_VoriPoly <- raster::extract(Raster_Layer, Maxima_Vori_Poly,
                                  fun =  function(x,...)GAP_DENSITY_FUN(na.omit(x), Para_BW = 0.4),
                                  na.rm = TRUE) 

  # OUTPUT DF
  max_XYZ_MedianVoriPoly <- data.frame(maxXY, maxZvalue, Median_Height_MaximLoc_VoriPoly, data.frame(Gap_VoriPoly[,1:2]))

  return(list(max_XYZ_MedianVoriPoly=max_XYZ_MedianVoriPoly, Maxima_Vori_Poly=Maxima_Vori_Poly))
}

##########################################################################################################################################
##########################################################################################################################################
VOXEL_FUN <- function(LAS_XYZ_ID, Para_Vox_Res = 1)
{
  
  # DETERMINING WHICH VOXEL EACH TID GOES INTO   
  Vox_XYZ <- voxel_metrics(LAS_XYZ_ID, list(PID), res = Para_Vox_Res) # Assign XYZ of Voxel to each Voxel
  colnames(Vox_XYZ) <- c("X", "Y", "Z", "PID")
  Vox_XYZ_Collapse <-Vox_XYZ[,list(PID = paste(unique(as.numeric(PID)),collapse = ',')),by = c("X", "Y", "Z")] # Gives each voxel one row and Points clisted in PID
  # browser()
  Vox_XYZ_Collapse <- cbind(VoxID=seq(1, nrow(Vox_XYZ_Collapse), 1), Vox_XYZ_Collapse) # GIVES EACH VOXEL A UNIQUE IDENTIFYER

  # Create vector of points ordered in same order as VoxID 
  Vox_XYZ_Collapse_PointID <-paste(Vox_XYZ_Collapse$PID, ",")
  Vox_XYZ_Collapse_PointID <-unlist(strsplit(Vox_XYZ_Collapse_PointID, ","))
  Vox_XYZ_Collapse_PointID <-as.numeric(Vox_XYZ_Collapse_PointID)

  # CALCULATING HOW MANY POINTS IN EACH VOXEL
  Vox_XYZ_Point_Count <- voxel_metrics(LAS_XYZ_ID, length(Z), res = Para_Vox_Res) 
  colnames(Vox_XYZ_Point_Count)[ncol(Vox_XYZ_Point_Count)] <- "Point_Count_In_Vox" 
  # browser()
  Vox_XYZ_Point_Count <- cbind(Vox_XYZ_Point_Count, PID= Vox_XYZ_Collapse$PID, VoxID=Vox_XYZ_Collapse$VoxID)

  # CREATING VECTOR THE LENGTH OF POINTS WITH VoxID FOR EACH POINT
  Vox_XYZ_Expand <- unlist(mapply(rep, Vox_XYZ_Collapse$VoxID,  Vox_XYZ_Point_Count$Point_Count_In_Vox))
  if(class(Vox_XYZ_Expand)[1] == "matrix"){Vox_XYZ_Expand <-as.vector(Vox_XYZ_Expand)}

  #ASSIGNING EACH POINT TO VOXEL and PUT VOX ID INTO STEM LAS FILE
  VoxID_TID_df <- data.frame(PID=Vox_XYZ_Collapse_PointID, VoxID=Vox_XYZ_Expand)
  LAS_XYZ_ID@data$VoxID[match(VoxID_TID_df$PID,LAS_XYZ_ID@data$PID)] <- VoxID_TID_df$VoxID
  Output <- list(LAS_XYZ_ID, Vox_XYZ_Point_Count, VoxID_TID_df)
}

##########################################################################################################################################
##########################################################################################################################################
# CONVEX HULL FUNCTION (GET TREE POLYGON)
CONVHULL_FUN = function(X, Y, TID)
{
  CHull = chull(X, Y)
  x_hul = X[CHull]
  y_hul = Y[CHull]
  Poly = Polygons(list(Polygon(cbind(x_hul,y_hul))), TID[1])
  return(list(X_Hull = list(x_hul), Y_Hull = list(y_hul), Poly = list(Poly)))
}

##########################################################################################################################################################################################
###################################################################################################################################################################################
GET_TREE_ATTRIBUTES_FUN = function(X, Y, Z, Classify,  Vox_PntCnt, Slice_Size_Zaxis, Unique_Poly_ID = 0)
{
  # BINNING TREE
  if((length(seq(min(Z), max(Z),Slice_Size_Zaxis))-1) >1) { # THIS IF STATEMENT GETS AROUND TREE HAVING ONE BIN
    Slice_Height <- cut(Z, seq(min(Z), max(Z),Slice_Size_Zaxis),labels=seq(min(Z), max(Z),Slice_Size_Zaxis)[-length(seq(min(Z), max(Z),Slice_Size_Zaxis))]) # labels=1:(length(seq(min(Z), max(Z),Slice_Size_Zaxis))-1))
  }else{
    Slice_Height <- rep(1, length(Z))
  }

  Slices_ID_DF <-na.omit(data.frame(X,Y,Z,Vox_PntCnt,Slice_Height))
  Unique_BIN_ID <- sort(as.numeric(unique(Slices_ID_DF$Slice_Height)))
  
  # Splitting data.frame into lists for each bin
  splits <- split(Slices_ID_DF, Slices_ID_DF$Slice_Height)
  Count_Splits <- length(sapply(splits, NROW))
  if(length(splits) != length(Unique_BIN_ID)){ # WORK AROUND WHEN THERE IS AN EMPTY LIST
    splits <- splits[- which(as.vector(sapply(splits, NROW)) ==0)]
  }
  
  # GET POLYGON OF EACH BIN AND AREA IT REPRESENTS
  Raster_Res <- 0.2
  Polygon_Slices <- lapply(splits, function(x) {
    # LAS SLICE
    LAS_oneSt_oneSl <- LAS(x)
  
    # RASTERISING THE SLICE LAS COUNT, SLICE 
    if(nrow(LAS_oneSt_oneSl@data)>1){
      R_oneSl_ZCount <-   grid_metrics(LAS_oneSt_oneSl, ~sum(Vox_PntCnt), res = Raster_Res)
      if(ncell(R_oneSl_ZCount) > 1){
        algo = lidR::watershed(R_oneSl_ZCount, tol = 1, ext=5) 
        LAS_oneSt_oneSl <- segment_trees(LAS_oneSt_oneSl, algorithm = algo, attribute = "WSID", uniqueness = "incremental")
        R_WS <- algo() 
        Unique_WS <- as.numeric(names(table(LAS_oneSt_oneSl@data$WSID)))
        Unique_WS <- Unique_WS[which(Unique_WS > 0)]
        
        if(length(Unique_WS) > 0){
          WS_Pnt_Count <- as.vector(table(LAS_oneSt_oneSl@data$WSID))
          WS_Grid_Count <-as.vector(table(getValues(R_WS)))
          
          # IF SOME WS GRIDS DON'T HAVE LAS POINTS (SLIGHT OFFSET OR ON EDGE DOESN'T REGISTER)
          if(length(WS_Pnt_Count) != length(WS_Grid_Count)){
            Remove_Length <- length(WS_Grid_Count) - length(WS_Pnt_Count)
            Remove_Value <- sort(WS_Grid_Count)[Remove_Length]
            Remove_Index <- which(WS_Grid_Count %in% Remove_Value)[Remove_Length]
            WS_Grid_Count <- WS_Grid_Count[-Remove_Index]
          }
  
          # STACK THE RASTERS
          # WORK AROUND FOR WHEN THE ROWS AND COLS ARE BUGGY (FLIPPED FOR UNKNOWN REASON)
          if((R_WS@ncols == R_oneSl_ZCount@ncols)&(R_WS@nrows == R_oneSl_ZCount@nrows)){
            R_Stack <- stack(R_oneSl_ZCount,  R_WS) # R_oneSl_ClID,
          }else{
            R_WS_Nrows <- R_WS@nrows
            R_WS_Ncols <- R_WS@ncols
            R_WS@ncols <-  R_WS_Nrows
            R_WS@nrows <-  R_WS_Ncols
            R_Stack <- stack(R_oneSl_ZCount,  R_WS) # R_oneSl_ClID,
          }
  
          R_Stack_DF <- data.frame(Cell_ID = seq(1, ncell(R_oneSl_ZCount),1) , getValues(R_Stack))
          colnames(R_Stack_DF) <- c("Cell_ID",
                                    "ZCount",
                                    "WS")
          
          R_Stack_DF <- data.frame(sapply(R_Stack_DF, function(x) as.numeric(as.character(x))))
          
          # LOCATE MAX DENSITY CELL_ID FOR EACH WS
          Cell_ID_maxDen <- as.data.frame(R_Stack_DF %>%
                                            dplyr::select(Cell_ID, ZCount, WS) %>%
                                            dplyr::group_by(WS) %>%
                                            dplyr::summarise(Cell_ID_maxDen = Cell_ID[which.max(ZCount)][1], .groups = 'drop'))
          
          Cell_ID_maxDen <- na.omit(Cell_ID_maxDen)
          
          if(nrow(Cell_ID_maxDen) > 0){
            Cell_ID_maxDen <- Cell_ID_maxDen[which(Cell_ID_maxDen$WS %in% Unique_WS),]
            
            XY_WS_MaxDen <- xyFromCell(R_WS, Cell_ID_maxDen$Cell_ID_maxDen, spatial=FALSE)
            
            # GET POLYGON OF EACH WS
            Poly_oneSl_allPoly <- SpatialPolygons(list(), proj4string = CRS(Proj_Sys))
            PolyDiamMax_oneSt_oneSl <- c()
            PolyDiamMin_oneSt_oneSl <- c()
            for(UWS in 1:length(Unique_WS)){
              
              Unique_Poly_ID <- Unique_Poly_ID + 1
              LAS_oneSt_oneSl_onePoly <- filter_poi(LAS_oneSt_oneSl, WSID == Unique_WS[UWS])
              
              # WHOLE SLICE WS
              Slice_CHull = CONVHULL_FUN(LAS_oneSt_oneSl_onePoly@data$X,
                                         LAS_oneSt_oneSl_onePoly@data$Y,
                                         Unique_Poly_ID)
              PolyDiamMax_oneSt_oneSl = c(PolyDiamMax_oneSt_oneSl, max(dist(cbind(unlist(Slice_CHull$X_Hull),unlist(Slice_CHull$Y_Hull)))))
              Poly_SP_oneSl_onePoly <- SpatialPolygons(Slice_CHull$Poly, proj4string = CRS(Proj_Sys))
              CoordPoly_oneSl_onePoly  <- Poly_SP_oneSl_onePoly@polygons[[1]]@Polygons[[1]]@coords
              Poly_oneSl_allPoly <- SpatialPolygons(c(slot(Poly_oneSl_allPoly,
                                                           "polygons"), Slice_CHull$Poly))
            }
            
            Poly_oneSl_allPolyArea <- unlist(lapply(Poly_oneSl_allPoly@polygons, function(x) slot(x, "area")))

            # MIN DIAMETER IS CALCULATED USING MAX DIAMETERS AND AREA IN AN ELIPSOID CALCULATION
            PolyDiamMin_oneSt_oneSl =Poly_oneSl_allPolyArea/as_units((pi * PolyDiamMax_oneSt_oneSl/2)*2, "m^2") 
            
            Poly_oneSl_allPolyArea_Total <- sum(Poly_oneSl_allPolyArea)
  
          }else{
            Poly_oneSl_allPoly <- SpatialPolygons(list())
            Poly_oneSl_allPolyArea <- 0
            Poly_oneSl_allPolyArea_Total <- 0
            PolyDiamMax_oneSt_oneSl <- 0
            PolyDiamMin_oneSt_oneSl <- 0
            XY_WS_MaxDen <- 0
            WS_Pnt_Count <- 0
            WS_Grid_Count <- 0
          }
        }else{
          Poly_oneSl_allPoly <- SpatialPolygons(list())
          Poly_oneSl_allPolyArea <- 0
          Poly_oneSl_allPolyArea_Total <- 0
          PolyDiamMax_oneSt_oneSl <- 0
          PolyDiamMin_oneSt_oneSl <- 0
          XY_WS_MaxDen <- 0
          WS_Pnt_Count <- 0
          WS_Grid_Count <- 0
        }
        
      }else{
        Poly_oneSl_allPoly <- SpatialPolygons(list())
        Poly_oneSl_allPolyArea <- 0
        Poly_oneSl_allPolyArea_Total <- 0
        PolyDiamMax_oneSt_oneSl <- 0
        PolyDiamMin_oneSt_oneSl <- 0
        XY_WS_MaxDen <- 0
        WS_Pnt_Count <- 0
        WS_Grid_Count <- 0
      }
      
    }else{
      Poly_oneSl_allPoly <- SpatialPolygons(list())
      Poly_oneSl_allPolyArea <- 0
      Poly_oneSl_allPolyArea_Total <- 0
      PolyDiamMax_oneSt_oneSl <- 0
      PolyDiamMin_oneSt_oneSl <- 0
      XY_WS_MaxDen <- 0
      WS_Pnt_Count <- 0
      WS_Grid_Count <- 0
      }
    
    list(Poly_oneSl_allPoly = Poly_oneSl_allPoly,
         Poly_oneSl_allPolyArea = Poly_oneSl_allPolyArea,
         Poly_oneSl_allPolyArea_Total = Poly_oneSl_allPolyArea_Total,
         PolyDiamMax_oneSt_oneSl = PolyDiamMax_oneSt_oneSl,
         PolyDiamMin_oneSt_oneSl = PolyDiamMin_oneSt_oneSl,
         XY_WS_MaxDen = XY_WS_MaxDen,
         WS_Pnt_Count = WS_Pnt_Count,
         WS_Grid_Count = WS_Grid_Count) 
  }) # END POLYGON FUNCTION
  
  ########################
  # EACH WS POLYGON OUTPUT
  ########################
  
  Centroids_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[6]), use.names = FALSE)
  
  # ID SLICES THAT WERE NOT COMPUTED DUE TO RASTER COUNT IN WS CALCULATION
  Centroids_allSl_allPoly_X <- Centroids_allSl_allPoly[which(Centroids_allSl_allPoly < 5000000)]
  Zero_Centroid<- which(Centroids_allSl_allPoly_X == 0)
  #browser()
  if(length(Zero_Centroid) >0){
    Centroids_allSl_allPoly_X <- Centroids_allSl_allPoly_X[-Zero_Centroid]
    Centroids_allSl_allPoly_Y <- Centroids_allSl_allPoly[which(Centroids_allSl_allPoly > 5000000)]
    Centroids_allSl_allPoly <-Centroids_allSl_allPoly [which(Centroids_allSl_allPoly == 0)]
    
    AreaPoly_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[2])   , use.names = TRUE)
    AreaPoly_allSl_allPoly_Height <-str_extract(names(AreaPoly_allSl_allPoly), "\\-*\\d+\\.*\\d*")[-Zero_Centroid]
    AreaPoly_allSl_allPoly <- as.vector(AreaPoly_allSl_allPoly)[-Zero_Centroid]
    
    DiamMax_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[4]), use.names = FALSE)[-Zero_Centroid]
    DiamMin_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[5]), use.names = FALSE)[-Zero_Centroid]
    WSPntCount_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[7]), use.names = FALSE) [-Zero_Centroid]
    WSGridCount_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[8]), use.names = FALSE) [-Zero_Centroid]
    Poly_eachWS <- unlist(sapply(Polygon_Slices,function(x) x[1]), use.names = TRUE)[-Zero_Centroid]

  }else { # ALL WS RASTERS HAVE BEEN COMPUTED SO NO SLICES NEED REMOVAL
    Centroids_allSl_allPoly_X <- Centroids_allSl_allPoly_X
    Centroids_allSl_allPoly_Y <- Centroids_allSl_allPoly[which(Centroids_allSl_allPoly > 5000000)]
    Centroids_allSl_allPoly <-Centroids_allSl_allPoly
    
    AreaPoly_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[2])   , use.names = TRUE)
    AreaPoly_allSl_allPoly_Height <-str_extract(names(AreaPoly_allSl_allPoly), "\\-*\\d+\\.*\\d*")
    AreaPoly_allSl_allPoly <- as.vector(AreaPoly_allSl_allPoly)
    
    DiamMax_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[4]), use.names = FALSE)
    DiamMin_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[5]), use.names = FALSE)
    WSPntCount_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[7]), use.names = FALSE)
    WSGridCount_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[8]), use.names = FALSE)
    Poly_eachWS <- unlist(sapply(Polygon_Slices,function(x) x[1]), use.names = TRUE)
    
    AreaSl_allSl_allPoly <- unlist(sapply(Polygon_Slices,function(x) x[3]) , use.names = TRUE)
    AreaSl_allSl_allPoly_Height <-str_extract(names(AreaSl_allSl_allPoly), "\\-*\\d+\\.*\\d*")
    AreaSl_allSl_allPoly <- as.vector(AreaSl_allSl_allPoly)
  }

  Output_eachWS  <- data.frame( TreeID = rep(Classify, length(AreaPoly_allSl_allPoly_Height)),
                                PolyHeight = as.numeric(as.character(AreaPoly_allSl_allPoly_Height)),
                                PolyArea = AreaPoly_allSl_allPoly,
                                PolyVolume = AreaPoly_allSl_allPoly * Slice_Size_Zaxis,
                                DiamMax = DiamMax_allSl_allPoly,
                                DiamMin = DiamMin_allSl_allPoly,
                                WS_PntCount = WSPntCount_allSl_allPoly,
                                WSGrid_Count = WSGridCount_allSl_allPoly,
                                WSGrid_Volume = WSGridCount_allSl_allPoly * (Raster_Res*Raster_Res*Slice_Size_Zaxis),
                                X_Cent = Centroids_allSl_allPoly_X,
                                Y_Cent = Centroids_allSl_allPoly_Y)
  
  
  Points_Centroids_XYZ <- Output_eachWS[, match(c("X_Cent", "Y_Cent", "PolyHeight", "TreeID"), colnames(Output_eachWS))]
  coordinates(Points_Centroids_XYZ) = ~X_Cent + Y_Cent
  crs(Points_Centroids_XYZ) <- crs(Proj_Sys)
  
  for(a in 1:length(Poly_eachWS)){
    if(a == 1){
      Poly_eachWS_All <- Poly_eachWS[[a]]
      PolyID_Max <- max(as.numeric(sapply(slot(Poly_eachWS_All, "polygons"), function(x) slot(x, "ID"))))
    }else{
      Poly_oneWS <-Poly_eachWS[[a]]
      new_IDs = as.numeric(sapply(slot(Poly_oneWS, "polygons"), function(x) slot(x, "ID"))) + PolyID_Max
      if(length(new_IDs) > 0){
        for (i in 1:length(slot(Poly_oneWS, "polygons"))){
          slot(slot(Poly_oneWS, "polygons")[[i]], "ID") = as.character(new_IDs[i])
        }
        Poly_eachWS_All <- rbind(Poly_eachWS_All, Poly_oneWS)
        PolyID_Max <- max(as.numeric(sapply(slot(Poly_eachWS_All, "polygons"), function(x) slot(x, "ID"))))
      }else{
        
        # WORK AROUND FOR WHEN A Poly_eachWS IS EMPTY AND YOU NEED TO REMOVE A ROW FROM THE OUTPUT DF TO MATCH IT
        Remove_Output_Row <- Output_eachWS[which(Output_eachWS$PolyHeight == unique(Output_eachWS$PolyHeight)[a]),]
        Remove_Output_Row_Area <- Remove_Output_Row$PolyArea[which(Remove_Output_Row$PolyArea == min(Remove_Output_Row$PolyArea))]
        Index_Remove_Output_Row <- which(Output_eachWS$PolyHeight == unique(Output_eachWS$PolyHeight)[a] &
                                           Output_eachWS$PolyArea == Remove_Output_Row_Area)
        if(length(Index_Remove_Output_Row) > 0){
          Output_eachWS <- Output_eachWS[-Index_Remove_Output_Row,]
        }
      }
    }
  }

  rname_eachWS <- row.names(Output_eachWS)
  rname_Poly <- sapply(slot(Poly_eachWS_All, "polygons"), function(x) slot(x, "ID"))
  
  setdiff_RowNames_Output_eachWS_missing <- setdiff(rname_eachWS, rname_Poly)
  if(length(setdiff_RowNames_Output_eachWS_missing) > 0 ){
    Index_Remove1 <- which( rname_eachWS %in% setdiff_RowNames_Output_eachWS_missing)
    Output_eachWS <- Output_eachWS[-Index_Remove1,]
  }
  
  setdiff_RowNames_Poly_eachWS_All_missing <- setdiff(rname_Poly, rname_eachWS)
  if(length(setdiff_RowNames_Poly_eachWS_All_missing) > 0 ){
    Index_Remove2 <- which( rname_Poly %in% setdiff_RowNames_Poly_eachWS_All_missing)
    Poly_eachWS_All <- Poly_eachWS_All[-Index_Remove2,]
  }

  Poly_eachWS_All_SPDF <- SpatialPolygonsDataFrame(Poly_eachWS_All, Output_eachWS)
  
  # OUTPUT
  return(list(Points_Centroids_XYZ = list(Points_Centroids_XYZ),
              Poly_eachWS_All_SPDF = list(Poly_eachWS_All_SPDF))) 
}

##########################################################################################################################################################################################
##########################################################################################################################################################################################
# Stem_PCA_Funciton
STEM_PCA_FUN = function(One_Br, Para_Slice_Height, Up_Down = "UP"){ #, Type_Tree
 
  PCA_Output <- data.frame(TreeID=numeric(),
                           X=numeric(),
                           Y=numeric(),
                           Z=numeric())

  STEM_PCA <- data.frame(TreeID=numeric(),
                           X=numeric(),
                           Y=numeric(),
                           Z=numeric())

    if(Up_Down == "UP"){
      Min_Tree <- min(One_Br$Z)
      Max_Tree <- max(One_Br$Z) + Para_Slice_Height#[S]
      TreeID <- One_Br[1,1]
      }
    if(Up_Down == "DOWN"){
      Min_Tree <- min(One_Br$Z) - Para_Slice_Height#[S]
      Max_Tree <- max(One_Br$Z)
      TreeID <- One_Br[1,1]
    }

    mean_xyz <- apply(One_Br[,match(c("X", "Y", "Z") ,colnames(One_Br))], 2, mean) #
    xyz_pca   <- princomp(One_Br[,match(c("X", "Y", "Z") ,colnames(One_Br))])
    dirVector <- xyz_pca$loadings[, 1]   # PC1
    #Comp_1 <- xyz_pca$score[,1]

    # Calculating coordinates of the top and bottom of each tree (endpts)
    t_ends <- c((Min_Tree-mean_xyz[3])/dirVector[3], (Max_Tree-mean_xyz[3])/dirVector[3])
    endpts <- rbind(mean_xyz + t_ends[[1]]*dirVector, mean_xyz + t_ends[[2]]*dirVector)
    colnames(endpts) <- c("X", "Y", "Z")

    # OUTPUTTING PCA DIRECTIONAL VECTOR AND MEAN XYZ
    Output_MeanXYZ_dirVector <- rbind(mean_xyz, dirVector)
    Output_MeanXYZ_dirVector <- cbind(TreeID = rep(TreeID,2),Output_MeanXYZ_dirVector)
    PCA_Output <- rbind(PCA_Output, Output_MeanXYZ_dirVector)

    # FINAL OUTPUT
    STEM_PCA_one <- data.frame(TreeID = rep(TreeID,2), endpts)
    STEM_PCA <- rbind(STEM_PCA, STEM_PCA_one)
    #}
return(list(STEM_PCA=STEM_PCA, PCA_Output=PCA_Output))
} # END STEM_PCA_FUN_Simple

##########################################################################################################################################################################################
##########################################################################################################################################################################################
# GET RANGE OF ALL RUNS WITH CONSECUTIVE NUMBERS WITHIN A VECTOR
CONSECUTIVE_RUN_RANGE <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x) %in% 1:2) as.character(x) else paste0(x[1], "-", x[length(x)])
  }), use.names=FALSE)
}

##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
# UAS_PLOTTING FUNCTION

plot_UAS_function <- function(LIDAR,
                          Tree_Attribute_Table,
                          PCA_DBSCAN,
                          Vol_per_minDen,
                          Type)
{

plot(LIDAR, color="DBSCAN_ID")
rgl.bg(color = "white")
rgl.bbox(color=c("#333377","white"), emission="#333377",
         specular="#3333FF", shininess=5, alpha=0.8 )

  # TEXT
text3d((mean(LIDAR@data$X)-2),(mean(LIDAR@data$Y)-2),(max(LIDAR@data$Z)+4),
       paste("T:", unique(LIDAR@data$DBSCAN_ID)[tt]," ESP:",ESP[r]," DEN:",Den[d], sep=""), adj=c(0, 0), col="black", cex=1)
text3d((mean(LIDAR@data$X)-2),(mean(LIDAR@data$Y)-2),(max(LIDAR@data$Z+2)),
       paste(length(unique(PCA_DBSCAN$All_STEM_Slice_PCA$Classify)), " ", Type, " Clusters: ", sep=""), adj=c(0, 0), col="black", cex=2)
text3d((mean(LIDAR@data$X)-2),(mean(LIDAR@data$Y)-2),(min(LIDAR@data$Z)-2),
       paste(Vol_per_minDen, " Pnt per cm^3: ", sep=""), adj=c(0, 0), col="black", cex=1)
text3d(PCA_DBSCAN$All_STEM_Slice_PCA$llx[seq(2,nrow(PCA_DBSCAN$All_STEM_Slice_PCA),2)],
       PCA_DBSCAN$All_STEM_Slice_PCA$lly[seq(2,nrow(PCA_DBSCAN$All_STEM_Slice_PCA),2)],
       PCA_DBSCAN$All_STEM_Slice_PCA$llz[seq(2,nrow(PCA_DBSCAN$All_STEM_Slice_PCA),2)],
       PCA_DBSCAN$All_STEM_Slice_PCA$Classify[seq(2,nrow(PCA_DBSCAN$All_STEM_Slice_PCA),2)], adj=c(0, 0), col="black", cex=2, lwd=3)

# POLYGONS
One_Tree_Polygon_SP <- SpatialPolygons(Tree_Attribute_Table$TreePoly, proj4string = CRS(Projection_System) )
for(p in 1:length(One_Tree_Polygon_SP)){
  One_Tree_Polygon_SP_Vertices_XY <- st_coordinates(One_Tree_Polygon_SP) 
  polygon3d(One_Tree_Polygon_SP_Vertices_XY[,1],One_Tree_Polygon_SP_Vertices_XY[,2], rep(max(Tree_Attribute_Table$TreeTop_Z[p]), nrow(One_Tree_Polygon_SP_Vertices_XY)), fill = FALSE, col=Tree_Attribute_Table[1,p+1], lwd=6)
  }
segments3d(PCA_DBSCAN$All_STEM_Slice_PCA[,2:4], lwd=6, col= "red")

}

##########################################################################################################################################################################################
# PLOT FUNCTION
###############

plot_function <- function(LIDAR = LiDAR_Over_Mid_All ,
                          ID_Over = ID_Poly_Over_All[Plot_Trees],
                          ID_AllTree = unique(LiDAR_Over_Mid_All$WSID),
                          M_Title = LiDAR_Over_Mid_All,
                          LiDAR_Type = "UAV")
{

  # GETTING COLOUR SCHEME
  Greens <- colorRampPalette(c("chartreuse","darkgreen","green", "darkseagreen4","darkolivegreen1", "green4"))
  Greens <- Greens(30)
  RedYellowBlue <- colorRampPalette(c("red", "yellow", "purple", "pink", "yellow", "purple", "red", "yellow", "purple", "pink", "yellow", "purple"))
  RedYellowBlue <- RedYellowBlue(20)
  GreyScale_Light <- colorRampPalette(c("grey40","grey80"))
  GreyScale_Light <- GreyScale_Light(8)

  GreyScale_Dark<- colorRampPalette(c("grey0","grey40"))
  GreyScale_Dark <- GreyScale_Dark(8)

  length_All <- length(ID_AllTree)
  if(length_All == 1){length_All <- 2 }
  AllTree_Colours <- rep(RedYellowBlue, ceiling(length_All/length(RedYellowBlue)))[1:length_All]
  if(ID_Over[1] != ""){
    AllTree_Colours[match(ID_Over, ID_AllTree)] <- rep(Greens, ceiling(length(ID_Over)/length(Greens)))[1:length(ID_Over)]
  }

  NoID_Points <- which(substr(as.character(ID_AllTree),1,5) == "NoID_" )

  if(length(NoID_Points > 0)){
    AllTree_Colours[NoID_Points] <- rep(GreyScale_Dark, ceiling(length(NoID_Points )/length(GreyScale_Dark)))[1:length(NoID_Points )]
  }

  #palette(AllTree_Colours)
  cols = AllTree_Colours[match(LIDAR$WSID,ID_AllTree)]
  open3d()

  # xlim, ylim and zlim
  if(LiDAR_Type != "UAV"){
    X_plotCoord <- LIDAR@data$X 
    Y_plotCoord <- LIDAR@data$Y 
    Z_plotCoord <- LIDAR@data$Z
  }else{
    X_plotCoord <- LIDAR@coords[,1] 
    Y_plotCoord <- LIDAR@coords[,2] 
    Z_plotCoord <- LIDAR@coords[,3]
  }
  LiDAR_Coords <-data.frame(X_plotCoord, Y_plotCoord, Z_plotCoord)
  colnames(LiDAR_Coords) <- c("X", "Y", "Z")

  range_X <- max(X_plotCoord)-min(X_plotCoord)
  range_Y <- max(Y_plotCoord)-min(Y_plotCoord)

  min_Range <- min(range_X, range_Y)
  max_Range <- max(range_X, range_Y)
  Range_dif <- max_Range-min_Range

  if(range_X == max_Range){
    New_X_Min <-  min(X_plotCoord)
    New_X_Max <-  max(X_plotCoord)
    New_Y_Min <-  min(Y_plotCoord)- Range_dif/2
    New_Y_Max <-  max(Y_plotCoord)+ Range_dif/2
  }else{
    New_X_Min <-  min(X_plotCoord)- Range_dif/2
    New_X_Max <-  max(X_plotCoord)+ Range_dif/2
    New_Y_Min <-  min(Y_plotCoord)
    New_Y_Max <-  max(Y_plotCoord)
  }
  plot3d(LiDAR_Coords, col=cols, main= M_Title, ylim = c(New_Y_Min,New_Y_Max),xlim=c(New_X_Min, New_X_Max))
}
