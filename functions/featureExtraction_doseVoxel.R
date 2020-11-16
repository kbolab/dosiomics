source("./functions/01_statisticalFeatures_dose.R")
source("./functions/03_GLCMtexturalFeatures_dose.R")
source("./functions/03_GLCMtexturalFeatures2-5D.R")
source("./functions/03_GLCMtexturalFeatures25Dmerged.R")
source("./functions/03_GLCMtexturalFeatures_dose.R")
source("./functions/03_GLCMtexturalFeaturesMerged.R")
source("./functions/05_GLSZMtexturalFeatures_dose.R")
source("./functions/04_GLRLMtexturalFeatures2-5D.R")
source("./functions/04_GLRLMtexturalFeatures.R")
source("./functions/04_GLRLMtexturalFeatures25Dmerged.R")
source("./functions/04_GLRLMtexturalFeaturesMerged.R")
source("./functions/05_GLSZMtexturalFeatures_dose.R")
source("./functions/05_GLSZMtexturalFeatures2-5D.R")

source("./functions/services.R")

dose.feature.extraction <- function(obj.dose ,feature.family = c("stat","morph","glcm","rlm","szm"), pixelSpacing = c(1,1,1), n_grey= 100){

  def <- c()
  
  if("stat" %in% feature.family) {
    cat("computing stat features...\n")
    stat.df <- statisticalFeatures(obj.dose$voxelCube.Dose)
    def <- c(def,stat.df)
  }
  # 
  if("glcm" %in% feature.family){
    cat("computing glcm features...\n")
    F_cm <-  glcmTexturalFeatures(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_cm <- do.call(data.frame,lapply(F_cm, function(x) replace(x, is.infinite(x),NA)))
    cm.df <- colMeans(F_cm,na.rm = T)
    def <- c(def,cm.df)
  }
  # 
  if("glcm" %in% feature.family){
    cat("computing glcm merged features...\n")
    F_cm_merged <-  glcmTexturalFeaturesMerged(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_cm_merged <- do.call(data.frame,lapply(F_cm_merged, function(x) replace(x, is.infinite(x),NA)))
    cm_merged.df <- colMeans(F_cm_merged,na.rm = T)
    def <- c(def,cm_merged.df)
  }
  
## da errore durante il calcolo del chiasma per pz 62330
  
  if("glcm" %in% feature.family){
    cat("computing glcm 2.5D features...\n")
    F_cm_25D <-  glcmTexturalFeatures25D(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_cm_25D <- do.call(data.frame,lapply(F_cm_25D, function(x) replace(x, is.infinite(x),NA)))
    cm_25d.df <- colMeans(F_cm_25D,na.rm = T)
    def <- c(def,cm_25d.df)
  }

  if("glcm" %in% feature.family){
    cat("computing glcm 2.5D merged features...\n")
    F_cm_25D_merged <-  glcmTexturalFeatures25Dmerged(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_cm_25D_merged <- do.call(data.frame,lapply(F_cm_25D_merged, function(x) replace(x, is.infinite(x),NA)))
    def <- c(def,F_cm_25D_merged)
  }

  if("rlm" %in% feature.family) {
    cat("computing glrm features...\n")
    F_rlm <-glrlmTexturalFeatures(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_rlm <- do.call(data.frame,lapply(F_rlm, function(x) replace(x, is.infinite(x),NA)))
    rlm.df <- colMeans(F_rlm,na.rm = T)
    def <- c(def,rlm.df)
  }
  #
  if("rlm" %in% feature.family) {
    cat("computing glrm merged features...\n")
    F_rlm_merged <-glrlmTexturalFeaturesMerged(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_rlm_merged <- do.call(data.frame,lapply(F_rlm_merged, function(x) replace(x, is.infinite(x),NA)))
    rlm_merged.df <- colMeans(F_rlm_merged,na.rm = T)
    def <- c(def,rlm_merged.df)
  }

  if("rlm" %in% feature.family) {
    cat("computing glrm 2.5D features...\n")
    F_rlm_25D <-glrlmTexturalFeatures25D(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_rlm_25D <- do.call(data.frame,lapply(F_rlm_25D, function(x) replace(x, is.infinite(x),NA)))
    rlm_25D.df <- colMeans(F_rlm_25D,na.rm = T)
    def <- c(def,rlm_25D.df)
  }
  #
  if("rlm" %in% feature.family) {
    cat("computing glrm 2.5D merged features...\n")
    F_rlm_25D_merged <-glrlmTexturalFeatures25Dmerged(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_rlm_25D_merged <- do.call(data.frame,lapply(F_rlm_25D_merged, function(x) replace(x, is.infinite(x),NA)))
    def <- c(def,F_rlm_25D_merged)
  }

  if("szm" %in% feature.family) {
    cat("computing szm features...\n")
    F_szm <- glszmTexturalFeatures(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_szm <- do.call(data.frame,lapply(F_szm, function(x) replace(x, is.infinite(x),NA)))
    szm.df <- colMeans(F_szm,na.rm = T)
    def <- c(def,szm.df)
  }
  #
  if("szm" %in% feature.family) {
    cat("computing szm 2.5D features...\n")
    F_szm_25D <- glszmTexturalFeatures25D(obj.dose$voxelCube.Dose, n_grey=n_grey)
    F_szm_25D <- do.call(data.frame,lapply(F_szm_25D, function(x) replace(x, is.infinite(x),NA)))
    szm_25D.df <- colMeans(F_szm_25D,na.rm = T)
   def <- c(def,szm_25D.df)
  }

  return(def)
}