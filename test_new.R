library(plyr)
library(reshape2)
library(moddicom)
#library("misc3d", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
source("./functions/featureExtraction_doseVoxel.R")

# -----------------------------------------------------------------------------
# calcolo features

obj <- geoLet()
obj$openDICOMFolder(pathToOpen = "test/dose")

obj$getROIList()
roi <- obj$getROIVoxels("PTV_RA")
pix.sp <- c(roi$geometricalInformationOfImages$pixelSpacing,as.numeric(roi$geometricalInformationOfImages$SliceThickness))
dose <- obj$extractDoseVoxels(ROIName = "PTV_RA")
features <- dose.feature.extraction(obj.dose = dose, feature.family = c("stat","glcm","rlm","szm"), pixelSpacing =  pix.sp)