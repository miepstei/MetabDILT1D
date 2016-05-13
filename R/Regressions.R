
#' performs univariate fold change regressions for each metabolite between two time points in DILT1D study
#' \deqn{\Delta M_i = X\Beta _i + \Epsilon _i}
#'
#' @param covariates The DILT1D covariate data frame
#' @param cSampleInfo The DILT1D dataframe for sample metadata
#' @param cMetDataLong The DILT1D dataframe for metabolite counts (or z-scores)
#' @param timepoints A vector of visits with which to calculate the fold changes
#' @param the predictor variable for the regression
#'
#' @importFrom plyr "."
#'
#' @return A list of coefficents, standard errors, t-values, p-values and the models for the regressions
#'
#' @export

RegressFoldchange <- function(covariates, cSampleInfo, cMetDataLong, timepoints = c("V0Pre","V0Post"), predictor = "DIFF"){

  if( !length(timepoints == 2 )){
    stop("Need two timepoints in regress_foldchange")
  }

  subject_lookup <- merge(covariates[,c("sex","strategyNEW","age_V0c","SUBJECT_ID")],cSampleInfo[,c("SUBJECT_ID","SAMPLE_NAME","VISIT"),drop=F], by="SUBJECT_ID")
  requiredSamples <- dplyr::filter(cSampleInfo,VISIT %in% timepoints)$SAMPLE_NAME

  #all the metabolite data for the two timepoints
  folds <- dplyr::filter(cMetDataLong, SAMPLE_NAME %in% requiredSamples)

  #select only relevent metabolites replresented in all selected samples
  missingMets <- plyr::ddply(folds, "METABOLITE_NAME", plyr::summarise, TIME_MISSING = MetabolonR::HasNA(METABOLITE_COUNT))
  notMissing <- dplyr::filter(missingMets, TIME_MISSING==0)$METABOLITE_NAME
  foldCompleteMets <- dplyr::filter(folds, METABOLITE_NAME %in% notMissing)
  foldCompleteMets <- plyr::ddply(foldCompleteMets, "METABOLITE_NAME", transform, METABOLITE_SCALED = MetabolonR::MeanCentred(log10(METABOLITE_COUNT)))
  foldCompleteMetsAndDose <- merge(subject_lookup, foldCompleteMets, by = "SAMPLE_NAME")
  foldCompleteMetsAndDose <- reshape2::dcast(foldCompleteMetsAndDose, METABOLITE_NAME + SUBJECT_ID + sex + strategyNEW + age_V0c ~ VISIT, value.var="METABOLITE_SCALED")
  foldCompleteMetsAndDose$DIFF <- foldCompleteMetsAndDose[, timepoints[[2]] ] - foldCompleteMetsAndDose[, timepoints[[1]] ]

  regress <- plyr::dlply(foldCompleteMetsAndDose, .(METABOLITE_NAME), lm, formula = as.formula (paste("DIFF", predictor,sep="~")))

  coeffs <- function (x) t(summary(x)$coefficients[,1,drop=F])
  se <-function (x) t(summary(x)$coefficients[,2,drop=F]) #drop is necessary in the case of 1_D
  tval <-function (x) t(summary(x)$coefficients[,3,drop=F])
  pval <-function (x)  t(summary(x)$coefficients[,4,drop=F])

  coefficents <- plyr::ldply(regress,coeffs)
  standardErrors <- plyr::ldply(regress,se)
  tStats <- plyr::ldply(regress,tval)
  pvals <- plyr::ldply(regress,pval)

  return(list(coefficents,standardErrors,tStats,pvals,regress))

}


#' performs univariate crosssectional regressions for each metabolite at a given timepoint in DILT1D study
#' \deqn{M_i = X\Beta _i + \Epsilon _i}
#'
#' @param covariates The DILT1D covariate data frame
#' @param cSampleInfo The DILT1D dataframe for sample metadata
#' @param cMetDataLong The DILT1D dataframe for metabolite counts (or z-scores)
#' @param timepoint A visit with which to calculate the crosssectional regression
#' @param the predictor variable for the regression
#'
#' @importFrom plyr "."
#'
#' @return  A list of coefficents, standard errors, tvalues and p values for the regressions
#'
#' @export


RegressCrosssection <- function(covariates, cSampleInfo, cMetDataLong, timepoint = "V0Pre", predictor ="factor(sex)"){
  #### Merge Datasets ####

  subject_lookup <- merge(covariates[,c("sex","strategyNEW","age_V0c","SUBJECT_ID")],cSampleInfo[,c("SUBJECT_ID","SAMPLE_NAME"), drop=F], by="SUBJECT_ID")

  #rownames(subject_lookup) <- rownames(sampleInfo[sampletype.DILT1D,])
  requiredSamples <- dplyr::filter(cSampleInfo,VISIT == timepoint)$SAMPLE_NAME
  crosssection <- cMetDataLong[cMetDataLong$SAMPLE_NAME %in% requiredSamples,]

  #select only relevent metabolites replresented in all selected samples
  missingMets <- plyr::ddply(crosssection, "METABOLITE_NAME", plyr::summarise, TIME_MISSING=MetabolonR::HasNA(METABOLITE_COUNT))
  notMissing <- missingMets[missingMets$TIME_MISSING == 0,]$METABOLITE_NAME
  crosssectionCompleteMets <- crosssection[ crosssection$METABOLITE_NAME %in% notMissing,]

  crosssectionCompleteMets <- plyr::ddply(crosssectionCompleteMets, "METABOLITE_NAME", transform, METABOLITE_SCALED=MetabolonR::MeanCentred(METABOLITE_COUNT))

  crosssectionCompleteMetsAndDose <- merge(subject_lookup, crosssectionCompleteMets, by = "SAMPLE_NAME")
  regress <- plyr::dlply(crosssectionCompleteMetsAndDose, .(METABOLITE_NAME), lm, formula = as.formula(paste("METABOLITE_SCALED", predictor, sep="~")))

  se <-function (x) summary(x)$coefficients[,2]
  tval <-function (x) summary(x)$coefficients[,3]
  pval <-function (x) summary(x)$coefficients[,4]

  coefficents <- plyr::ldply(regress,coef)
  standardErrors <- plyr::ldply(regress,se)
  tStats <- plyr::ldply(regress,tval)
  pvals <- plyr::ldply(regress,pval)

  return(list(coefficents,standardErrors,tStats,pvals,regress))

}
