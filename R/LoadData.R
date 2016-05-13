
#' Loads in a list of metabolites that are significantly differentially
#' expressed between sexes. This list is a result of a mWAS analysis
#' by Krumsiek et al
#'
#' @param sexMappingFile A character filename that documents a mapping between the
#' metabolites in the Krumsiek panel (3rd generation) and the Metabolites in DILT1D
#' (4th generation)
#' @param cMetInfo A Data Frame of metabolite metadata
#' @param cMetDataLong A Data Frame of metabolite counts
#'
#' @importFrom magrittr "%>%"
#'
#' @return A list of three Data Frames; the metabolite info of the sex metabolites;
#' the data for the sex metabolites in long format and the Krumsief lookup which has
#' the summary statistics from the analysis
#' @export



SexMetabolitesFromKrumsiek <- function(sexMappingFile = "Data/KrumsiekMapping.csv", field = "RAW", cMetInfo, cMetDataLong) {

  genderMetsKrum <- read.csv(sexMappingFile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  knownSexMets <- (dplyr::select_(genderMetsKrum, field) %>% na.omit %>% unlist %>% as.vector)

  cSMetInfo <- dplyr::filter(cMetInfo, BIOCHEMICAL %in%  knownSexMets)

  #DIL study biochemical data
  cSMetDataLong <- dplyr::filter(cMetDataLong, (METABOLITE_NAME %in% cSMetInfo$METABOLITE_NAME))

  #Krumsiek study merged with DIL biochemical data metabolite name for reference
  krumLookup <- merge(genderMetsKrum , dplyr::select(cMetInfo, c(METABOLITE_NAME,BIOCHEMICAL)),by.x=field, by.y="BIOCHEMICAL", all.y = TRUE)

  return(list(cSMetInfo,cSMetDataLong, krumLookup))

}


#' Loads the covariate data from the ipswich network share for DILT1D
#'
#' @param covariatesFileStem A filestem on the ipswich share that documents covariates for the
#' DILT1D trial participants
#' @param covariatesDate The date for the covariates file for which to choose the covariates
#' @param covariatesSourceFile An R file of functions to be sourced to retrieve the data, also on ipswich
#'
#' @return A Data Frame of covariate data for the DILT1D trial
#'
#' @examples
#' covariateData <- LoadDILT1DCovariates(covariatesFile = "/ipswich/path/to/file-stem-", covariatesDate = 'YYYY-MM-DD' covariatesFile = "/ipswich/path/to/file.R")
#' @export

LoadDILT1DCovariates <- function(covariatesFileStem, covariatesDate, covariatesSourceFile  ){

  #date.covariates="2014-04-07"
  #source("/ipswich/data/shared/DILT1D/RFunctions/readdoses.R")
  #dir.doses.file <- "/ipswich/data/shared/DILT1D/lookups/phenotype-"

  source(covariatesSourceFile)
  covariates <- prepare.doses.function(dir=covariatesFileStem, date=covariatesDate)
  covariates <- doses.formodelling(covariates,which.strategy="5")

  #rename trialid -> SUBJECT_NAME for consistency
  names(covariates)[names(covariates)=="trialid"] <- "SUBJECT_ID"
  covariates$strategyNEW <- MetabolonR::AsNumericFactor(covariates$strategyNEW)
  covariates$sex <- factor(covariates$sex,levels=c("M","F"))
  covariates <- dplyr::mutate(covariates,age_V0c = as.vector(scale(age_V0,center = T,scale = F)))
  return(covariates)
}


#' Loads the DILT1D dataset given the relevent sample, metabolite and count file
#' Cleans the files which have fields specific to the DILT1D output from Metabolon
#'
#' @param sampleFile A character filepath to the sample metadata
#' @param metaboliteFile A character filepath to the metabolite metadata
#' @param dataFile A character filepath to the count data
#'
#' @return A list of Data Frames with the sample and metabolite metadata and the
#' volume normalised count data (but not further normalised)
#'
#' @importFrom magrittr "%>%"
#'
#' @export

LoadDILT1DDataVolNormalised <- function(sampleFile, metaboliteFile, dataFile) {
  STANDARD_VOLUME <- 95

  data <- MetabolonR::LoadMetabolonData(sampleFile, metaboliteFile, dataFile)
  clean_data <- CleanMetabolonData (data[[1]],data[[2]],data[[3]])

  cSampleInfo <- clean_data[[1]]
  cMetInfo <- clean_data[[2]]
  cMetData <- clean_data[[3]]
  cMetDataLong <- clean_data[[4]]

  #generating dil dataset

  nonStandardMetaboliteNames <- cMetInfo %>% dplyr::select(METABOLITE_NAME) %>% unlist %>% as.vector
  dilt1dSamples <- dplyr::filter(cSampleInfo, SAMPLE_TYPE == "DILT1D") %>% dplyr::select(SAMPLE_NAME) %>% unlist %>% as.vector

  dilt1dData <- cMetDataLong %>% dplyr::filter(SAMPLE_NAME %in% dilt1dSamples) %>% dplyr::filter(METABOLITE_NAME %in% nonStandardMetaboliteNames)
  dilt1dSampleInfo <- cSampleInfo %>% dplyr::filter(SAMPLE_NAME %in% dilt1dSamples)

  #Volume normalise data here
  dilt1dDataNormalised <- MetabolonR::VolumeNormaliseDataset(dilt1dSampleInfo, cMetInfo, dilt1dData, STANDARD_VOLUME, "DILT1D")

  return(list(dilt1dSampleInfo, cMetInfo, dilt1dDataNormalised))
}

#' Loads the DILT1D dataset given the relevent sample, metabolite and count file
#' Cleans the files which have fields specific to the DILT1D output from Metabolon
#'
#' @param sampleFile A character filepath to the sample metadata
#' @param metaboliteFile A character filepath to the metabolite metadata
#' @param dataFile A character filepath to the count data
#'
#' @return A list of Data Frames with the sample and metabolite metadata and the
#' volume and median normalised count data
#' @importFrom magrittr "%>%"
#'
#' @export

LoadDILT1DData <- function(sampleFile, metaboliteFile, dataFile) {

  #these groups form the actual days over which the samples were analysed by Metabolon
  RUNDAY_GROUPS <- list(c(1, 2, 3), c(4, 5, 6, 7), c(8, 9, 10, 11))

  dilt1dDataVolNormalised <- LoadDILT1DDataVolNormalised(sampleFile, metaboliteFile, dataFile)

  dilt1dSampleInfo <- dilt1dDataVolNormalised[[1]]
  metaboliteInfo <- dilt1dDataVolNormalised[[2]]
  dilt1dDataVolNormalised <- dilt1dDataVolNormalised[[3]]

  #Median normalise data here
  dilt1dDataNormalised <- MetabolonR::MedianNormaliseDataset(dilt1dSampleInfo, dilt1dDataVolNormalised, "RUN_DAY", RUNDAY_GROUPS)

  return(list(dilt1dSampleInfo, metaboliteInfo, dilt1dDataNormalised))
}

#' Cleans the untidy dataset given the relevent sample, metabolite and count
#' Data Frames. These Frames have fields specific to the DILT1D xls output from Metabolon
#'
#' @param sampleInfo A Data Frame of sample metadata
#' @param metaboliteInfo A Data Frame of metabolite metadata
#' @param metaboliteData A Data Frame of count data
#'
#' @return A list of 4 Data Frames with the sample and metabolite metadata and the
#' Wide and Long versions of the Metabolite Count data
#' @export


CleanMetabolonData <- function(sampleInfo, metaboliteInfo, metaboliteData) {

  sampletype.POOLED<-grepl('MTRX',sampleInfo$SAMPLE_TYPE)
  sampletype.DILT1D <- as.logical(grepl('DIL',sampleInfo$PARAM_BOX) & grepl('EXPERIMENTAL',sampleInfo$SAMPLE_TYPE))
  sampletype.DGAP <- as.logical(grepl('DGAP',sampleInfo$PARAM_BOX)  & grepl('EXPERIMENTAL',sampleInfo$SAMPLE_TYPE))

  if (!all(sampletype.POOLED+sampletype.DILT1D+sampletype.DGAP)){
    #check that all the samples are partitioned
    stop("Discontiguous separation of IS, DILT1D and DGAP samples")
  }

  sampleInfo$SAMPLE_TYPE<-""
  sampleInfo$SAMPLE_TYPE[sampletype.POOLED]<-"POOLED"
  sampleInfo$SAMPLE_TYPE[sampletype.DILT1D]<-"DILT1D"
  sampleInfo$SAMPLE_TYPE[sampletype.DGAP]<-"DGAP"

  #clear up periods and empty strings in these columns
  metaboliteInfo[,c("SUPER_PATHWAY","SUB_PATHWAY","CAS","KEGG","HMDB")]<-apply(metaboliteInfo[c("SUPER_PATHWAY","SUB_PATHWAY","CAS","KEGG","HMDB")],2,function(x) gsub('\\.|^$',NA,x) )

  #reorder columns
  metaboliteInfoTidy<-cbind(metaboliteInfo[,c("METABOLITE_NAME","METABOLITE_TYPE")] ,metaboliteInfo[,!(colnames(metaboliteInfo) %in% c("METABOLITE_NAME","METABOLITE_TYPE"))])

  #need to tidy up the data about the sample information
  sampleInfo$PARAM_CLIENT_VOLUME_ML<-gsub('\\~| uL|\\.','',sampleInfo$PARAM_CLIENT_VOLUME_ML,perl = TRUE)

  #need to evaluate some expressions that have need entered as strings in the excel spreadsheet:
  sampleInfo$PARAM_CLIENT_VOLUME_ML[grepl('\\+',sampleInfo$PARAM_CLIENT_VOLUME_ML)]<-as.vector(sapply(sampleInfo$PARAM_CLIENT_VOLUME_ML[grepl('\\+',sampleInfo$PARAM_CLIENT_VOLUME_ML)],function(x) eval(parse(text=x))))
  sampleInfo$PARAM_CLIENT_VOLUME_ML<-as.numeric(sampleInfo$PARAM_CLIENT_VOLUME_ML)

  #need to tidy up and remove illegal strings from the PARAM fields
  sampleInfo[,colnames(dplyr::select(sampleInfo,contains("PARAM")))] <- apply(dplyr::select(sampleInfo,contains("PARAM")),2,function(x) gsub('^\\.$',NA,x) )

  #we are happy with warnings being surpressed in the following statement
  sampleInfo$PARAM_VOLUME_EXTRACTED_UL<-suppressWarnings(as.numeric(sampleInfo$PARAM_VOLUME_EXTRACTED_UL))

  #remove '.' and 'na' for missing values and replace with NA in the VISIT column

  sampleInfo[,c("PARAM_VISIT")] <- gsub('na',NA,sampleInfo[,c("PARAM_VISIT")])
  sampleInfo[,c("PARAM_VISIT")] <- gsub('\\.',NA,sampleInfo[,c("PARAM_VISIT")])
  sampleInfo$PARAM_VISIT <- factor(sampleInfo$PARAM_VISIT, levels = c('V0Pre','V0Post','V1','V2','V3','V4','V5','V6'))

  #order the levels in the RUN_DAY appropriately
  sampleInfo$PARAM_RUN_DAY <- factor( sampleInfo$PARAM_RUN_DAY ,levels = gtools::mixedsort(unique(sampleInfo$PARAM_RUN_DAY)))

  sampleInfoTidy<-sampleInfo[,c("SAMPLE_NAME","SAMPLE_TYPE","PARAM_BOX","PARAM_CLIENT_VOLUME_ML","PARAM_RUN_DAY","PARAM_SUBJECT_ID","PARAM_VISIT","PARAM_VOLUME_EXTRACTED_UL","SAMPLE_ID")]

  #add the following if they exist
  if ("PARAM_WELL" %in% colnames(sampleInfo) ){
    sampleInfoTidy <- cbind(sampleInfoTidy, dplyr::select(sampleInfo, PARAM_WELL))
  }

  if ("PARAM_LC_COLUMN" %in% colnames(sampleInfo)){
    sampleInfoTidy <- cbind(sampleInfoTidy, dplyr::select(sampleInfo, PARAM_LC_COLUMN))
  }

  if("CLIENT_IDENTIFIER" %in% colnames(sampleInfo)) {
    sampleInfoTidy <- cbind(sampleInfoTidy, dplyr::select(sampleInfo, CLIENT_IDENTIFIER))
  }

  #tidy up the column names
  colnames(sampleInfoTidy) <- gsub('^PARAM_','',colnames(sampleInfoTidy))

  #data in long format
  metDataLong <- reshape2::melt(metaboliteData, id.vars="SAMPLE_NAME")
  colnames(metDataLong) <- c("SAMPLE_NAME", "METABOLITE_NAME", "METABOLITE_COUNT")

  return(list(sampleInfoTidy,metaboliteInfoTidy,metaboliteData,metDataLong))
}


#' returns the CPeptide measurements for subjects in the DILT1D trial
#'
#' @param cPeptideFile A character file location for the CPeptide measurements
#'
#' @return A Data Frame of the C-peptide measurements for each subject per visit
#'
#' @export
#' @importFrom tidyr gather


readDILT1DCpeptideData <- function(cPeptideFile){

  cpepdata=read.table(cPeptideFile, header=TRUE, stringsAsFactors=FALSE)

  names(cpepdata)[names(cpepdata)=="trialid"] <-"SUBJECT_ID"
  names(cpepdata)[names(cpepdata)=="p1"] <- "C_Peptide,pmol/l"
  colnames(cpepdata) <- toupper(colnames(cpepdata))

  cpepdata$VISIT <- factor(cpepdata$VISIT)
  cpepdata$VISIT <- factor(cpepdata$VISIT, levels=c("Sc","V9","V10"))
  colnames(cpepdata) <- make.names(colnames(cpepdata))

  toChange <- reshape2::dcast(cpepdata, SUBJECT_ID~VISIT, value.var = "ELAPSED_DAYS") %>% dplyr::filter(is.na(V10)) %>% dplyr::select(SUBJECT_ID) %>% unlist %>% as.vector
  cpepdata[cpepdata$VISIT=="V9" & cpepdata$SUBJECT_ID %in% toChange,]$VISIT <- "V10"

  cpepdata
}

#' returns the HbA1c measurements for subjects in the DILT1D trial
#'
#' @param hbaFile A character file location for the HbA1c measurements
#'
#' @return A Data Frame of the HbA1c measurements for each subject per visit
#'
#' @export

readDILT1DHBAData <- function(hbaFile) {
  #date.file="2014-06-11"
  #dir.hba.file = "/ipswich/data/shared/DILT1D/datasets/clinical-HbA1c/HbA1c-results-by-visit-"
  hbadat = read.table(hbaFile, header = TRUE, stringsAsFactors=FALSE)
  names(hbadat)[names(hbadat)=="trialid"] <-"SUBJECT_ID"
  names(hbadat)[names(hbadat)=="a1"] <- "Clinical-HbA1c"
  colnames(hbadat) <- toupper(colnames(hbadat))

  hbadat$VISIT <- factor(hbadat$VISIT)
  hbadat$VISIT <- factor(hbadat$VISIT, levels=c("Sc","V9","V10"))
  colnames(hbadat) <- make.names(colnames(hbadat))

  toChange <- reshape2::dcast(hbadat, SUBJECT_ID~VISIT, value.var = "ELAPSED_DAYS") %>% dplyr::filter(is.na(V10)) %>% dplyr::select(SUBJECT_ID) %>% unlist %>% as.vector
  hbadat[hbadat$VISIT=="V9" & hbadat$SUBJECT_ID %in% toChange,]$VISIT <- "V10"

  hbadat
}

#' returns the Self Measured Blood Glucose (SMBG) measurements for subjects in the DILT1D trial
#'
#' @param smbgFile A character file location for the SMBG measurements
#'
#' @return A Data Frame of the SMBG measurements for each subject per visit
#'
#' @importFrom magrittr "%>%"
#'
#' @export

readDILT1DSMBGData <- function(smbgFile){
  #smbg - self-measured blood glucose

  smbgData = read.table(smbgFile, header = TRUE, stringsAsFactors=FALSE)
  names(smbgData)[names(smbgData)=="trialid"] <-"SUBJECT_ID"

  #standardise the levels and names of the VISIT variable
  colnames(smbgData) <- toupper(colnames(smbgData))
  smbgData$VISIT <- factor(smbgData$VISIT)
  smbgData$VISIT <- plyr::revalue(smbgData$VISIT, c("V0pre"="V0Pre", "V0post"="V0Post"))
  smbgData$VISIT <- factor(smbgData$VISIT, levels=c("Sc", "V0Pre","V0Post","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))

  #End V10 timepoint, agglomerate V9,V10 and remeasure analytes
  #need to do some reworking on the smbg,hba, cpeptide datasets to create the correct endpoints
  #due to effective abandoning of V5
  top <- reshape2::dcast(smbgData, SUBJECT_ID~VISIT, value.var = "ELAPSEDDAY") %>% dplyr::filter(!is.na(V10))
  bottom <- reshape2::dcast(smbgData, SUBJECT_ID~VISIT, value.var = "ELAPSEDDAY") %>% dplyr::filter(is.na(V10))
  smbg_shift <- rbind(top , MetabolonR::RenameColumns(bottom, c("V5","V6","V7","V8","V9"), c("V6","V7","V8","V9","V10")))
  smbg_shift <- reshape2::melt(smbg_shift,id.vars="SUBJECT_ID", variable.name="VISIT", value.name="ELAPSED_DAYS")

  top <- reshape2::dcast(smbgData,SUBJECT_ID~VISIT,value.var = "SMBG") %>% dplyr::filter(!is.na(V10))
  bottom <- reshape2::dcast(smbgData,SUBJECT_ID~VISIT,value.var = "SMBG") %>% dplyr::filter(is.na(V10))
  smbg_shift2 <- rbind(top , MetabolonR::RenameColumns(bottom, c("V5","V6","V7","V8","V9"), c("V6","V7","V8","V9","V10")))
  smbg_shift2 <- reshape2::melt(smbg_shift2,id.vars="SUBJECT_ID", variable.name="VISIT", value.name="SMBG")
  smbgDatFixed <- dplyr::arrange(plyr::join(smbg_shift, smbg_shift2,by=c("SUBJECT_ID","VISIT")), SUBJECT_ID)

  smbgDatFixed
}

#' returns the Clinical Chemistry measurements (ClinChem) for subjects in the DILT1D trial
#'
#' @param smbgFile A character file location for the ClinChem measurements
#'
#' @return A Data Frame of the ClinChem measurements for each subject per visit
#'
#' @importFrom magrittr "%>%"
#'
#' @export

readDILT1DClinchemData <- function(clinChemFile){

  ccdat = read.table(clinChemFile,header = TRUE,stringsAsFactors=FALSE)

  refsFile <- read.csv("/ipswich/data/shared/DILT1D/datasets/clinical-Biochemistry/biochemistry-references.csv",header = TRUE,stringsAsFactors=FALSE)
  colnames(ccdat)[which(grepl("^c\\d+",colnames(ccdat)))] <- refsFile$reportline

  names(ccdat)[names(ccdat)=="trialid"] <-"SUBJECT_ID"
  colnames(ccdat) <- toupper(colnames(ccdat))

  ccdat$VISIT <- factor(ccdat$VISIT)
  ccdat$VISIT <- plyr::revalue(ccdat$VISIT, c("V0pre"="V0Pre"))
  ccdat$VISIT <- factor(ccdat$VISIT, levels=c("Sc", "V0Pre","V1","V2","V3","V5","V6","V7","V8","V9","V10"))
  colnames(ccdat) <- make.names(colnames(ccdat))
  ccdat <- tidyr::gather(ccdat,CLINCHEM_MEASUREMENT,CLINCHEM_VALUE,ALT..U.L:UREA..MMOL.L)
  ccdat
}

#' returns a particular Clinical Chemistry analyte from the ClinChem dataframe for subjects in the DILT1D trial
#'
#' @param clinChem A Data Frame of the ClinChem measurements for each subject per visit
#'
#' @return A character vector of the analyto to be returned
#'
#' @importFrom magrittr "%>%"
#'
#' @export

fetchClinchemAnalyte <- function(clinChem, analyte){

  analyteData <- dplyr::filter(clinChem, CLINCHEM_MEASUREMENT == analyte) #"GLUCOSE..MMOL.L")
  top <- reshape2::dcast(analyteData,SUBJECT_ID~VISIT, value.var = "ELAPSED_DAYS") %>% dplyr::filter(!is.na(V10))
  bottom <- reshape2::dcast(analyteData,SUBJECT_ID~VISIT, value.var = "ELAPSED_DAYS") %>% dplyr::filter(is.na(V10))
  analyteShift <- rbind(top, MetabolonR::RenameColumns(bottom, c("V5","V6","V7","V8","V9"), c("V6","V7","V8","V9","V10")))
  analyteShift <- reshape2::melt(analyteShift, id.vars="SUBJECT_ID", variable.name="VISIT", value.name="ELAPSED_DAYS")

  top <- reshape2::dcast(analyteData,SUBJECT_ID~VISIT, value.var = "CLINCHEM_VALUE") %>% dplyr::filter(!is.na(V10))
  bottom <- reshape2::dcast(analyteData,SUBJECT_ID~VISIT, value.var = "CLINCHEM_VALUE") %>% dplyr::filter(is.na(V10))
  analyteShift2 <- rbind(top, MetabolonR::RenameColumns(bottom, c("V5","V6","V7","V8","V9"), c("V6","V7","V8","V9","V10")))
  analyteShift2 <- reshape2::melt(analyteShift2, id.vars="SUBJECT_ID", variable.name="VISIT", value.name="CLINCHEM_VALUE")

  analyteDatFixed <- dplyr::arrange(plyr::join(analyteShift, analyteShift2, by=c("SUBJECT_ID","VISIT")), SUBJECT_ID)
  analyteDatFixed
}

