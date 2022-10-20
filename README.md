# Individual tree detection and crown delineation from Unmanned Aircraft System (UAS) LiDAR in structurally complex mixed species eucalypt forests

The Algorithm represents the following published work (https://www.sciencedirect.com/science/article/pii/S0924271620302963):

Jaskierniak, D., A. Lucieer, G. Kuczera, D. Turner, P. N. J. Lane, R. G. Benyon, and S. Haydon. 2021. 'Individual tree detection and crown delineation from Unmanned Aircraft System (UAS) LiDAR in structurally complex mixed species eucalypt forests', ISPRS Journal of Photogrammetry and Remote Sensing, 171: 171-87.

![image](https://user-images.githubusercontent.com/22090635/137672008-502176d3-b208-4028-ae20-0cfb4bac3298.png)

https://user-images.githubusercontent.com/22090635/196833925-58ff1e34-fa53-4f21-a112-2a219b7eb309.mp4

**Figure 1:** Example ITCD with typical forest conditions using red to represent unclassified low density points, brown-to-green colour scale to represent understorey points,  and all other colours to represent detected trees, where: (a) shows intertwining canopy (IC) as well as suppressed (S) and mid-storey (MS) trees, (b) shows leaning trees (L), overstorey canopy base touching understorey (CB) and MS, (c) shows IC, S and MS, (d) shows IC and dead stags (D), (e) shows MS, L and S, and (f-h) shows deteriorating point densities down the vegetation profile (PD) under circumstances with CB, S, D and IC.

# Procedure for running algorithm

The Individual tree and crown delineation (ITCD) procedure requires the execution of two files: 

**1. ITCD_ALGO.R** performs a point cloud segmentation procedure. The algorithm first uses kernel densities to stratify the vegetation profile and differentiate understorey from the rest of the vegetation. For vegetation above understorey, the ITCD algorithm adopted a novel watershed clustering procedure on point density measures within horizontal slices. A Principal Component Analysis (PCA) procedure is then applied to merge the slice-specific clusters into trunks, branches, and canopy clumps, before a voxel connectivity procedure clustered these biomass segments into overstorey trees. The algorithm makes use of all the customised functions found in FUNCTIONS_ITCD.R 
To run this algorithm the following information and data is required:
* **Output_Dir:** DIRECTORY FOR OUTPUTTING THE SEGMENTATION POINT CLOUD
* **Data_Dir:**  DIRECTORY FOR PLOT CENTROIDS FILE (CSV), WHICH CONTAIN THE FOLLOWING THREE COLUMNS   c("Flight", "POINT_X", "POINT_Y")
* **LAS_Dir:** DIRECTORY FOR LAS FILES REPRESENTING EACH UAS FLIGHT
* **FUN_Dir:** DIRECTORY FOR FUNCTIONS_ITCD.R FILE 

**2. OTIMISER_UAS_SCE.R** is used to calibrate two site-specific parameters that capture the localised conditions using a Shuffled Complex Evolution (SCE) optimiser. For the purpose of applying ITCD_ALGO.R to sites without field measurements, Jaskierniak et al. (2021) provides regression relationships for predicting the two site-specific parameters using structural characteristics of vegetation clusters derived from ITCD_ALGO.R. This approach generalises the algorithm across new UAS LiDAR data without undertaking time-consuming ground measurements within tall eucalypt forests with complex vegetation structure. 
To run this algorithm, the following information and data is required:
* **FOLDER_ITCD:** DIRECTORY THAT HAS ALL OUTPUT FROM ITCD_ALGO.R (i.e. Output_Dir)
* **FOLDER_FIELD_DATA_SHAPEFILE:** DIRECTORY CONTAINING THE FIELD  DATA
* **FILE_FIELD_DATA_PREFIX:** PREFIX FOR FIELD DATA FILE NAMES (i.e. prefix_1.shp WHERE 1 IS THE FLIGHT_ID CORRESPONDING TO THE PLOT CENTROID FILE’S “Flight”
* **FIELD DATA** IS A SHAPEFILE WITH ATTRIBUTE “TID” TO REPRESENT THE UNIQUE TREE IDENTIFIER FOR EACH STEM IN PLOT.

Example data is provided in DATA_EXAMPLE. Note that LAZ file is too big to upload on GITHUB. 

