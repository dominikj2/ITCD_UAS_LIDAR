rm(list=ls())
options(digits = 12)
library(gstat)
library(sp)
library(rgdal)
library(stats)
library(lidR)
library (pastecs)
library(lattice)
library (maptools)
library(boot)
library(plotrix)
library(fields)
library(VoxR)
library(raster)
library(spatstat)
library(stringr)
library(dplyr)
library(e1071)
library(rgeos)
library(rgl)
library(clusterCrit)
library(dbscan)
library(data.table)
library(SDraw)
library(MASS)
library(stringr)
library(scatterplot3d)
library(fossil)
library(deldir)
library(EBImage)
library(animation)
library(magick)
library(stringr)
library(spatialEco)
library(DescTools)
library(randomcoloR)
library(RColorBrewer)
library(rtop)

##########################################################################################################################################
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
##########################################################################################################################################
# SETTING UP DIRECTORIES
########################

FOLDER_ITCD <- # FOLDER THAT HAS ALL THE OUTPUT FROM RUNNING ITCD_ALGO.R 

Field_Data_Files <- list.files(FOLDER_ITCD, pattern = "Flight")
Field_Data_ID <- numextract_all(Field_Data_Files)

SCE_Folder <- "" # NAME OF FOLDER FOR OPTIMISERS OUTPUT
dir.create(file.path(FOLDER_ITCD, SCE_Folder), showWarnings = FALSE)

FOLDER_FIELD_DATA_SHAPEFILE <- ""# FOLDER CONTAINING FIELD DATA
FILE_FIELD_DATA_PREFIX <- "" # PREFIX FOR THE FIELD DATA (File names much be prefix_1.shp where 1 represents the Flight_ID)
    # NOTE: FIELD DATA NEEDS TO BE SHAPEFILE THAT CAN BE READ IN AS readOGR
        # SHAPEFILE NEEDS TO HAVE AN ATTRIBUTE "TID" WHICH REPRESENTS A UNIQUE TREE IDENTIFIER FOR EACH STEM IN PLOT

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
#  OPTIMISER FUNCTION FOR OMISSION/COMISSION RATES
Minimise_Function <- function(x){
  # COEFFICIENTS TO OPTIMISE
  Para_Range_z <- x[1]
  Para_min_z <- x[2]                                                                                                                        
  
  # CONSTRAINTS
  Index_TID_Minised <- which(TID_Stats$Range_z >=  Para_Range_z &
                               (TID_Stats$min_z)  <=  Para_min_z )   

  # IF NO TREES IDENTIFIED WITH THE CONSTRAINT GIVE HIGH MINIMISATION VALUE
  if(length(Index_TID_Minised) == 0){
    Minimise_Function <- nrow(TID_Stats)}
  
  # IF TREES SELECTED WITH CONSTRANT
  else{
    # SUBSET TREE STAT AND PREDICTED TREES USING SCE CONSTRAINT
    TID_Stats_Minimised <- TID_Stats[Index_TID_Minised,]
    TID_XY_Minimised  <- TID_Stats_Minimised[, match(c("St_Base_X", "St_Base_Y"), colnames(TID_Stats_Minimised))] 

    # CALCULATE NEAREST MANUALLY DETECT TREE TO PREDICTED TREE
    colnames(TID_XY_Minimised) <- c("X", "Y")

    Index_XY <- match(c("X", "Y"), colnames(Field_Loc_Data))
    nearest <- RANN::nn2(Field_Loc_Data[,Index_XY],TID_XY_Minimised, k=1, searchtype="standard", treetype="bd")

    TID_Minimised_DF<- data.frame(TID_minimised = TID_Stats_Minimised$TID, 
                                               TID_XY_Minimised, 
                                               FieldID_Nearest = Field_Loc_Data$TID[nearest$nn.idx], 
                                               Dist_TID_FieldID=nearest$nn.dists)
    
    # OMISSION AND COMMISSION RATES
    Commission <- as.data.frame(table(TID_Minimised_DF$FieldID_Nearest))
    colnames(Commission) <- c("Obs_Tree", "Pred_Count")
    Commission$Obs_Tree <- as.integer(as.character(Commission$Obs_Tree))
    Omitted_Values <- setdiff(Field_Loc_Data$TID, Commission$Obs_Tree)

    if(length(Omitted_Values) >0){
      Omitted <- data.frame(Obs_Tree=Omitted_Values, Pred_Count=rep(0,length(Omitted_Values)))
      Omit_Commission_df <- rbind(Commission, Omitted)
      if(nrow(Omit_Commission_df)>1){
        Omit_Commission_df<- Omit_Commission_df[order(Omit_Commission_df$Obs_Tree),]
        }
    }else{
      Omit_Commission_df <- Commission
      if(nrow(Omit_Commission_df)>1){
        Omit_Commission_df<- Omit_Commission_df[order(Omit_Commission_df$Obs_Tree),]
        }
    }

    ########################################################
    # CALCULATE F SCORES
    ########################################################

    
    TreeID_1to1  <- Omit_Commission_df$Obs_Tree[which(Omit_Commission_df$Pred_Count == 1)]
    TreeID_Omitted  <- Omit_Commission_df$Obs_Tree[which(Omit_Commission_df$Pred_Count == 0)]
    TreeID_Comm  <- Omit_Commission_df$Obs_Tree[which(Omit_Commission_df$Pred_Count > 1)]
    
    TP <- length(TreeID_1to1)
    FN  <- length(TreeID_Omitted)
    FP <- sum( Omit_Commission_df$Pred_Count[which(Omit_Commission_df$Obs_Tree %in% TreeID_Comm)])-length(TreeID_Comm) 
    r_score <- TP/(TP+FN)
    p_score<- TP/(TP+FP)
   
    # MINIMISE FUNCTION
    Minimise_Function <-  eval(parse(text=(Minimisation_Function_List[m]))) 
    
    assign("Omit_Commission_df", Omit_Commission_df, envir = .GlobalEnv) 
    assign("TID_Minimised_DF", TID_Minimised_DF, envir = .GlobalEnv) 
    assign("TID_Stats_Minimised", TID_Stats_Minimised, envir = .GlobalEnv) 
    return(Minimise_Function)
  }
} # END Minimise_Function


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
# LOOP THROUGH FLIGHTS/PLOTS
First_Loop <- 1
X_NA_Length <- c()
Y_NA_Length <- c()
for(v in 1:length(Field_Data_ID)){ 

  Flight_ID <- Field_Data_ID[v]
  Done <- which(Field_Data_ID == Flight_ID)
  
  LAS_Dir <- paste(FOLDER_ITCD, "/Flight_", Flight_ID, "/LAS", sep="")
  Files_Done <- list.files(LAS_Dir, pattern = paste("F", Flight_ID, "_LAS_TID_Within_30Plot.laz", sep=""))
  if(length(Files_Done) == 1){  
    Field_GISCorr_Shp <- readOGR(dsn = paste(FOLDER_ITCD, "/" ,FOLDER_FIELD_DATA_SHAPEFILE,  "/Flight_",Flight_ID , sep=""),
                                 layer = paste(FILE_FIELD_DATA_PREFIX, Flight_ID, sep=""))
    Field_DF <-as.data.frame(cbind(Field_GISCorr_Shp@coords, Field_GISCorr_Shp@data))
    colnames(Field_DF)[1:2] <- c("X","Y")

    # LAS FILES
    Input_Folder <- paste(FOLDER_ITCD, "/Flight_", Flight_ID, sep="")
    LAS_Folder_Input <- paste(Input_Folder, "/LAS", sep="")
    LAS_File <-paste(LAS_Folder_Input, "/F", Flight_ID, "_LAS_TID_Within_30Plot.laz", sep="")
    LAS_TID <- lidR::readLAS(LAS_File, select = "xyzp0")
    X_Shift <- min(LAS_TID@data$X)
    Y_Shift <- min(LAS_TID@data$Y)
    LAS_TID@data$X <- LAS_TID@data$X - X_Shift
    LAS_TID@data$Y <- LAS_TID@data$Y - Y_Shift
    Range_Z_LAS <- range(LAS_TID@data$Z)

    CSV_Folder <-paste(Input_Folder, "/CSV", sep="")
    TID_Stats_File <- list.files(CSV_Folder, pattern = "_St_Attributes_60X60.csv")
    TID_Stats <- read.csv(paste(CSV_Folder, "/", TID_Stats_File, sep=""))
    TID_Stats$St_Base_X <- TID_Stats$St_Base_X - X_Shift
    TID_Stats$St_Base_Y <- TID_Stats$St_Base_Y - Y_Shift

    Field_DF$X <- Field_DF$X- X_Shift
    Field_DF$Y <- Field_DF$Y- Y_Shift
    Field_Loc_Data <- Field_DF

    dir.create(file.path(paste(FOLDER_ITCD, "/",SCE_Folder, sep=""), paste("Flight_", Flight_ID, sep="")), showWarnings = FALSE)
    SCE_Optim_Folder <-paste(FOLDER_ITCD, "/", SCE_Folder , "/Flight_", Flight_ID, sep="")

    dir.create(file.path(SCE_Optim_Folder, "Plots"), showWarnings = FALSE)
    SCE_Optim_Folder_PDF <- paste(FOLDER_ITCD, "/", SCE_Folder , "/Flight_", Flight_ID,"/Plots", sep="")

    # CREATE OUTPUT FOLDERS
    dir.create(file.path(paste(SCE_Optim_Folder, "/", sep=""), "LAS"), showWarnings = FALSE)
    Output_LAS_Folder <- paste(SCE_Optim_Folder, "/LAS", sep="")

    # LIST OF MINIMISATION FUNCTIONS
    Minimisation_Function_List <- list()
    Minimisation_Function_List$FScore = " 1 - 2*((r_score*p_score)/(r_score+p_score)) "
    Minimiser_Name <- c("FScore")

    for(m in 1:length(Minimisation_Function_List)){

      SCE_Outputs <- data.frame(matrix(ncol = 13, nrow = 0))
      colnames(SCE_Outputs) <- c("Flight","Para_Set", "Minimiser",
                                 "Parameter_1", "Parameter_2",
                                 "Value", "ParConvergence", "Iterations", "Timeout")         
      
      All_Summary_Omit_Commission <- data.frame(matrix(ncol = 8, nrow = 0))
      colnames(All_Summary_Omit_Commission) <- c("Flight", "Para_Set", "Minimiser",
                                                 "Total_Omitted", "Total_Commission", "Total_1_to_1",
                                                 "Total_Field_St", "Percent_1_1")

      #########################################################################################################################################
      # SCE OPTIMISER
      ###############
     
      Output_SCEUA <-sceua(OFUN = Minimise_Function,
                           pars = c(median(TID_Stats$Range_z),
                                    median(TID_Stats$min_z)),
                           lower=c(min(TID_Stats$Range_z),
                                   min(TID_Stats$min_z)),
                           upper=c(max(TID_Stats$Range_z),
                                   max(TID_Stats$min_z)))

      #########################################################################################################################################

      # OUTPUTTING OPIMISED PREDICTION COUNT FOR EACH OBS_TREE
      write.table(Omit_Commission_df, file=paste(SCE_Optim_Folder, "/Omit_Commiss_",Flight_ID, "_",  Minimiser_Name[m],"_GISCorr_FSCORE.csv", sep=''), row.names=F, col.names=T, sep=",")
      TID_Stats_Minimised$St_Base_X <- TID_Stats_Minimised$St_Base_X + X_Shift
      TID_Stats_Minimised$St_Base_Y <- TID_Stats_Minimised$St_Base_Y + Y_Shift
      write.table(TID_Stats_Minimised, file=paste(SCE_Optim_Folder, "/TID_Stats_Minimised_",Flight_ID, "_",  Minimiser_Name[m],"_GISCorr_FSCORE.csv", sep=''), row.names=F, col.names=T, sep=",")
      TID_Minimised_DF$X <- TID_Minimised_DF$X + X_Shift
      TID_Minimised_DF$Y <- TID_Minimised_DF$Y + Y_Shift
      write.table(TID_Minimised_DF, file=paste(SCE_Optim_Folder, "/TID_Minimised_DF_",Flight_ID, "_",  Minimiser_Name[m],"_GISCorr_FSCORE.csv", sep=''), row.names=F, col.names=T, sep=",")

      # SUMMARISING COMMISSIONS AND OMISSIONS FOR OUTPUTTING
      Total_Omitted <- length(which(Omit_Commission_df$Pred_Count == 0))
      Total_Commission <- sum(Omit_Commission_df$Pred_Count[which(Omit_Commission_df$Pred_Count > 1)]-1)
      Total_1_to_1 <- length(which(Omit_Commission_df$Pred_Count == 1))
      Total_Field_St <- nrow(Field_Loc_Data)
      Percent_1_1 <- Total_1_to_1/Total_Field_St
      
      TreeID_Comm  <- Omit_Commission_df$Obs_Tree[which(Omit_Commission_df$Pred_Count > 1)]
      
      TP <- Total_1_to_1
      FN  <- Total_Omitted
      FP <- sum( Omit_Commission_df$Pred_Count[which(Omit_Commission_df$Obs_Tree %in% TreeID_Comm)])-length(TreeID_Comm) #length(TreeID_Comm) #
      r_score <- TP/(TP+FN)
      p_score<- TP/(TP+FP)
      F_score<- 2*((r_score*p_score)/(r_score+p_score)) 
      
      Summary_Omit_Commission<- data.frame(Flight= Flight_ID, Minimiser= Minimiser_Name[m], Total_Omitted, Total_Commission, Total_1_to_1, 
                                           Total_Field_St, Percent_1_1, F_score)

      #########################################################################################################################################
      # APPENDING ONE SCE OUTPUT TO DATA.FRAME.... OUTPUTTING AS ONE FILES (ALL PARAMETER SETS) FOR WHOLE FLIGHT
      Obs_Tree_Count <-  nrow(Field_Loc_Data)

      Parameter_Outputs <- c(Flight_ID, Obs_Tree_Count,  Minimiser_Name[m], Output_SCEUA$par, Output_SCEUA$value,
                             Output_SCEUA$convergence$parConvergence, Output_SCEUA$iterations,  Output_SCEUA$timeout)
      One_SCE_Outputs <- as.data.frame(t(Parameter_Outputs))
      colnames(One_SCE_Outputs) <- c("Flight", "Obs_Tree_Count", "Minimiser",
                                     "Parameter_1", "Parameter_2",
                                     "Value", "ParConvergence", "Iterations", "Timeout")                                                      
      SCE_Outputs <- rbind(SCE_Outputs, One_SCE_Outputs)

      # OUTPUTTING RESULT FOR ONE FLIGHT
      colnames(SCE_Outputs) <- c("Flight", "Obs_Tree_Count",  "Minimiser",
                                 "Parameter_1", "Parameter_2",
                                 "Value", "ParConvergence", "Iterations", "Timeout")                                                           

      write.table(SCE_Outputs, file=paste(SCE_Optim_Folder, "/SCE_Flight",Flight_ID, "_", Minimiser_Name[m],"_GISCorr_FSCORE.csv", sep=''), row.names=F, col.names=T, sep=",")


      write.table(Summary_Omit_Commission, file=paste(SCE_Optim_Folder, "/Summary_Omit_Commiss_",Flight_ID, "_", Minimiser_Name[m],"_GISCorr_FSCORE.csv", sep=''), row.names=F, col.names=T, sep=",")

      if(v == First_Loop){
        All_Output <- Summary_Omit_Commission
      }else{
        All_Output <- rbind(All_Output, Summary_Omit_Commission)
      }

      write.table(All_Output, file=paste(FOLDER_ITCD, "/", SCE_Folder, "/All_Summary_Omit_Commission_",  Minimiser_Name[m],"_GISCorr_FSCORE.csv", sep=''), row.names=F, col.names=T, sep=",")
    }  # LOOP MINIMISER

  #OUTPUT LAS REMOVED AND KEPT
  TID_Keep <- TID_Stats[which(TID_Stats$min_z < as.numeric(as.character(SCE_Outputs$Parameter_2)) & TID_Stats$Range_z > as.numeric(as.character(SCE_Outputs$Parameter_1))),]
  TID_Removed  <- TID_Stats[-(which(TID_Stats$min_z < as.numeric(as.character(SCE_Outputs$Parameter_2)) & TID_Stats$Range_z > as.numeric(as.character(SCE_Outputs$Parameter_1)))),]

  LAS_TID_Keep <- filter_poi(LAS_TID, TreeID  %in% TID_Keep$TID)
  LAS_TID_Removed  <- filter_poi(LAS_TID, TreeID  %in% TID_Removed$TID)

  writeLAS(LAS_TID_Keep, paste(SCE_Optim_Folder, "/LAS/F", Flight_ID ,"_LAS_TID_KEEP_OPTIM_GISCorr_FSCORE.laz",sep=''))
  writeLAS(LAS_TID_Removed, paste(SCE_Optim_Folder, "/LAS/F", Flight_ID ,"_LAS_TID_REMOVE_OPTIM_GISCorr_FSCORE.laz",sep=''))

  }
} # LOOP THROUGH FLIGHTS
