#' Performs A comparison of QC methods on Metabolon data
#'
#' The method reads in the data, cleans it, and estimates
#' the impact of different QC methods (NOMIS, CCMN and median)
#' on a pca and positive control (glucose and sex metabolites)
#'
#' @param fileLocationsList A list of named file locations
#' @param outputDir A character filepath store the output from the analysis
#'
#' @export

RunQCAnalysis <- function(fileLocationsList, outputDir) {
  PROP_MISSING <- 0.2

  message(paste0("Creating / cleaning Output directory at ", outputDir))
  MetabolonR::CreateOrClean(outputDir)

  #we want to load the data and clean it, and get the volume normalised DILT1D samples
  message("Loading Volume Normalised DILT1D data")
  dilt1dData <- LoadDILT1DDataVolNormalised(fileLocationsList$sampleFile, fileLocationsList$metaboliteFile, fileLocationsList$metaboliteData)

  sampleInfo <- dilt1dData[[1]]
  metaboliteInfo <- dilt1dData[[2]]
  volNormMetabCounts <- dilt1dData[[3]]

  message(paste0("Removing metabolites missing in more than ", 100*PROP_MISSING , "% of the samples"))
  DILMetDataLongProp <- MetabolonR::FilterMissingMets(volNormMetabCounts, PROP_MISSING)
  message("[TODO]: Volume normalised data")

  message("NOMIS normalising data")
  nomisNormalised <- MetabolonR::NOMISNormaliseDataset(metaboliteInfo, DILMetDataLongProp)
}


#' Performs the variance components analysis across DILT1D dataset
#'
#' The method reads in the data, cleans it, and estimates
#' a variance components model for the DILT1D dataset using
#' visits V0Pre, V4, V5 assuming drug washout
#' The output is placed in a directory specified by the user
#'
#' @param fileLocationsList A list of named file locations
#' @param outputDir A character filepath store the output from the analysis
#'
#' @export

RunVarianceComponentsAnalysis <- function (fileLocationsList, outputDir) {
  PROP_MISSING <- 0.2

  message(paste0("Creating / cleaning Output directory at ", outputDir))
  MetabolonR::CreateOrClean(outputDir)

  message("Loading and Normalising dataset")
  dilt1dData <-  LoadDILT1DData(fileLocationsList$sampleFile, fileLocationsList$metaboliteFile, fileLocationsList$metaboliteData)

  sampleInfo <- dilt1dData[[1]]
  metaboliteInfo <- dilt1dData[[2]]
  normMetabCounts <- dilt1dData[[3]]

  message(paste0("Removing metabolites missing in more than ", 100*PROP_MISSING , "% of the samples"))
  DILMetDataLongProp <- MetabolonR::FilterMissingMets(normMetabCounts, PROP_MISSING)
  #DILMetDataLongProp <- MetabolonR::TransformLog10Mets(DILMetDataLongProp)

  message(paste0("Loading DILT1D covariate data"))
  covariates <- LoadDILT1DCovariates(fileLocationsList$covariateStem, fileLocationsList$covariateDate, fileLocationsList$covariateSourceFile)

  message(paste0("Transforming count data to z-scores"))
  metDataLongZScore <- MetabolonR::TransformZscoreMets(DILMetDataLongProp)

  metaboliteNames <- unique(metDataLongZScore$METABOLITE_NAME)

  fixedIntercept <- vector(mode = "list", length = length(metaboliteNames))
  ICC <- vector(mode = "list", length = length(metaboliteNames))
  phi_B_BW <- vector(mode = "list", length = length(metaboliteNames))
  phi_B_T <- vector(mode = "list", length = length(metaboliteNames))
  residualVar <- vector(mode = "list", length = length(metaboliteNames))
  #residualConfidence <- matrix(nrow=length(metaboliteNames),ncol=2)
  models <- vector(mode = "list", length = length(metaboliteNames))

  b <- 0
  for (i in metaboliteNames) {
    b <- b + 1
    print(paste("processing" , i))
    data <- dplyr::filter(metDataLongZScore, METABOLITE_NAME==i)

    DILSubSampleInfo <- dplyr::filter(sampleInfo, VISIT %in% c("V0Pre","V4","V5")) %>% dplyr::select(SAMPLE_NAME,VISIT,SUBJECT_ID)
    DILSubSampleInfo <- dplyr::mutate(DILSubSampleInfo, TIME = as.numeric(ifelse(VISIT=="V0Pre",0,ifelse(VISIT=="V0Post",0.16,ifelse(VISIT=="V1",1,ifelse(VISIT=="V2",2,ifelse(VISIT=="V3",3,ifelse(VISIT=="V4",4,ifelse(VISIT=="V5",5,6)))))))))
    testdata <- plyr::join(DILSubSampleInfo, data,by = "SAMPLE_NAME")

    p1 <- ggplot2::ggplot(testdata,ggplot2::aes(x=factor(TIME), y=METABOLITE_COUNT)) + ggplot2::geom_boxplot()
    p2 <- ggplot2::ggplot(testdata,ggplot2::aes(x=SUBJECT_ID, y=METABOLITE_COUNT)) + ggplot2::geom_boxplot()

    #model.fit <- lmer(METABOLITE_COUNT ~ 1 + (TIME | SUBJECT_ID), data=testdata)
    model.fit <- lme4::lmer(METABOLITE_COUNT ~ 1 +  (1 | TIME) + (1|SUBJECT_ID) , data = testdata, REML = T)
    varComponents <- as.data.frame(lme4::VarCorr(model.fit))

    betweenSubjectVar <- varComponents %>% dplyr::filter(grp == "SUBJECT_ID", var1=="(Intercept)", is.na(var2)) %>% dplyr::select(vcov) %>% unlist %>% as.vector
    #withinSubjectVar <- varComponents %>% filter(grp == "SUBJECT_ID",var1=="TIME",is.na(var2)) %>% select(vcov) %>% unlist %>% as.vector
    withinSubjectVar <- varComponents %>% dplyr::filter(grp == "TIME",var1=="(Intercept)", is.na(var2)) %>% dplyr::select(vcov) %>% unlist %>% as.vector
    technicalVariability <- varComponents %>% dplyr::filter(grp == "Residual") %>% dplyr::select(vcov) %>% unlist %>% as.vector

    totalVar <- betweenSubjectVar + withinSubjectVar + technicalVariability

    ICC[[b]] <- (betweenSubjectVar + withinSubjectVar) / totalVar
    phi_B_BW[[b]] <- betweenSubjectVar / (betweenSubjectVar + withinSubjectVar)
    phi_B_T[[b]] <- betweenSubjectVar / totalVar
    residualVar[[b]] <- technicalVariability
    models[[b]] <- model.fit

    #save the fixed effects here
    fixedIntercept[[b]] <- unname(lme4::fixef(model.fit))

    #bootstrap confidence intervals for residual variance parameters
    #a<-confint(model.fit,"sigma",method="boot",oldNames=FALSE,level=0.95,nsim = 100)
    #residualConfidence[b,1] <- 1#a[1,1]
    #residualConfidence[b,2] <- 2#a[1,2]
  }

  ICCs <- data.frame(METABOLITE_NAME=metaboliteNames,ICC=unlist(ICC), phi_B_BW = unlist(phi_B_BW), phi_B_T = unlist(phi_B_T), FIXED_INTERCEPT=unlist(fixedIntercept), MODEL_VARIANCE = unlist(residualVar))
  ICCs <- plyr::join(ICCs, dplyr::select(metaboliteInfo, METABOLITE_NAME, SUPER_PATHWAY, SUB_PATHWAY, BIOCHEMICAL, PATHWAY_SORTORDER), by="METABOLITE_NAME")
  ICCs <- dplyr::mutate(ICCs, MODEL_SD = sqrt(MODEL_VARIANCE))

  plotList <- PlotMetaboliteVarianceComponents(ICCs)

  ggsave(file.path(outputDir, "DILT1D_Individual_of_Total.pdf"), plotList$MedianIndividual, width = 7, height = 7)
  ggsave(file.path(outputDir, "DILT1D_Between_Individual_Variation_Of_Total.pdf"), plotList$MedianBetweenIndividualOfTotalIndividual, width = 7, height = 7)
  ggsave(file.path(outputDir, "DILT1D_Between_Individual_Variation.pdf"), plotList$MedianBetweenIndividualOfTotal, width = 7, height = 7)
  ggsave(file.path(outputDir, "DILT1D_Technical_Variation.pdf"), plotList$MedianTechnicalVariation, width = 7, height = 7)

}


#' Runs the DILT1D univariate analysis across visits
#'
#' The method reads in the data, cleans it, and estimates
#' a univariate model for each metabolite \deqn{y = X\Beta + \Epsilon}
#' between the V0Pre visit and subsequent timepoints. The output
#' is placed in a directory specified by the user
#'
#' @param fileLocationsList A list of named file locations
#' @param outputDir A character filepath store the output from the analysis
#'
#' @export

RunUniVariateAnalysis <- function (fileLocationsList, outputDir) {

  PROP_MISSING <- 0.2

  message(paste0("Creating / cleaning Output directory at ", outputDir))
  MetabolonR::CreateOrClean(outputDir)

  message("Loading and Normalising dataset")
  dilt1dData <-  LoadDILT1DData(fileLocationsList$sampleFile, fileLocationsList$metaboliteFile, fileLocationsList$metaboliteData)

  sampleInfo <- dilt1dData[[1]]
  metaboliteInfo <- dilt1dData[[2]]
  normMetabCounts <- dilt1dData[[3]]

  message(paste0("Removing metabolites missing in more than ", 100*PROP_MISSING , "% of the samples"))
  DILMetDataLongProp <- MetabolonR::FilterMissingMets(normMetabCounts, PROP_MISSING)
  DILMetDataLongProp <- MetabolonR::TransformLog10Mets(DILMetDataLongProp)

  message(paste0("Loading DILT1D covariate data from ", fileLocationsList$covariateStem , fileLocationsList$covariateDate))
  covariates <- LoadDILT1DCovariates(fileLocationsList$covariateStem, fileLocationsList$covariateDate, fileLocationsList$covariateSourceFile)

  #plot the covariate data here
  message("Plotting covariate data")
  p <- PlotCovariates(covariates)
  ggplot2::ggsave(filename = file.path(outputDir, "Covariates.pdf"), p, width = 7, height = 7)

  #load sex specific metabolites for positive control
  message("Deriving sex specific metabolites from Krumsiek et al.")
  sexMetabolites <- SexMetabolitesFromKrumsiek(fileLocationsList$krumsiekFile, field = "RAW", cMetInfo = metaboliteInfo, cMetDataLong = DILMetDataLongProp)
  sexRegression <- RegressCrosssection(covariates, sampleInfo, sexMetabolites[[2]], timepoint = "V0Pre", predictor ="factor(sex)")

  p <- PlotSexMetaboliteComparison(sexRegression[[1]], sexRegression[[4]], sexMetabolites[[3]], metaboliteInfo)[[1]]
  ggplot2::ggsave(filename = file.path(outputDir, "SexMetComparison.pdf"), p, width = 7, height = 7)


  #run fold change analysis
  message("Running fold change analysis over patient visits")

  #remove sex metabolites from fold change analysis
  sexMetaboliteNames <- sexMetabolites[[1]]$METABOLITE_NAME
  neutralMetInfo <- dplyr::filter(metaboliteInfo, !(METABOLITE_NAME %in% sexMetaboliteNames))
  neutralMetData <- dplyr::filter(DILMetDataLongProp, !(METABOLITE_NAME %in% sexMetaboliteNames))

  # Now perform the fold change regressions and volcano plots
  times = list(c("V0Pre", "V0Post"), c("V0Post", "V1"), c("V1", "V2"), c("V2", "V3"), c("V3", "V4"), c("V4", "V5"))
  timeVolcanos <- NULL
  for (time in times) {
    #fold change regressions
    foldChangeComplexModel <- RegressFoldchange(covariates, sampleInfo, neutralMetData, timepoints = time, predictor = "age_V0c + strategyNEW + factor(sex)")

    #plot_volcano
    volcano <- PlotVolcano(foldChangeComplexModel, "(Intercept)", neutralMetInfo)
    ggplot2::ggsave(filename = file.path(outputDir, paste0(paste0(time, collapse = ""), "VolcanoPlot.pdf")), volcano[[1]], width = 7, height = 7)
    foldChg <- volcano[[3]]
    foldChg <- dplyr::mutate(foldChg, TIMEPOINT = paste0(time, collapse = "-"))
    timeVolcanos <- rbind(timeVolcanos, foldChg)
  }

  #make a faceted plot of the fold changes
  timeVolcanos$TIMEPOINT <- factor(timeVolcanos$TIMEPOINT, levels=c("V0Pre-V0Post","V0Post-V1","V1-V2","V2-V3","V3-V4","V4-V5"))
  p <- ggplot2::ggplot(timeVolcanos, ggplot2::aes(x = LOG_2_FOLD,y = M_LOG10P_VAL)) + ggplot2::geom_point() + ggplot2::facet_wrap(~TIMEPOINT , nrow = 2)
  p <- p + ggplot2::theme(text = ggplot2::element_text(size=20), plot.title = ggplot2::element_text(size = 46), legend.position = "right")
  p <- p + ggplot2::ylab("-log10(p-value)") + ggplot2::xlab("log2 Fold Change")
  p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "longdash")
  p <- p + theme(axis.text=element_text(size=16),
             axis.title=element_text(size=24), legend.text=element_text(size=24), legend.title=element_text(size=24))


  ggplot2::ggsave(filename = file.path(outputDir, "VolcanoEvolution.pdf"),  p, width = 12, height = 8)

  #pathway analysis
  message("Running metabolic pathway analysis for fold changes in metabolites for V0Pre to V0Post")
  foldChgV0PreV0Post <- dplyr::filter(timeVolcanos, TIMEPOINT=="V0Pre-V0Post")
  foldChgV0PreV0Post <- dplyr::left_join(foldChgV0PreV0Post, dplyr::select(metaboliteInfo, METABOLITE_NAME, SUPER_PATHWAY, SUB_PATHWAY, KEGG), by = "METABOLITE_NAME")
  foldChgV0PreV0Post <- dplyr::mutate(foldChgV0PreV0Post, SUPER_PATHWAY = ifelse(is.na(SUPER_PATHWAY), "Unknown", SUPER_PATHWAY))

  #overall by super_pathway
  p <- ggplot2::ggplot(foldChgV0PreV0Post, ggplot2::aes(x = LOG_2_FOLD, y = M_LOG10P_VAL))
  p <- p + ggplot2::geom_point(ggplot2::aes(colour = SUPER_PATHWAY, alpha = SIG)) + ggplot2::facet_wrap(~SUPER_PATHWAY , nrow = 2)
  p <- p + ggplot2::scale_alpha_discrete(c(0.1, 1), guide = FALSE) + ggplot2::theme(legend.position = "bottom")
  p <- p + ggplot2::scale_colour_discrete(guide = FALSE)
  p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dotted")
  p <- p + ggplot2::theme(text = ggplot2::element_text(size=24), plot.title = ggplot2::element_text(size = 32), legend.position = "bottom")
  p <- p + ggplot2::ylab("-log10(p-value)") + ggplot2::xlab("log2 Fold Change") + scale_x_continuous(breaks = c(-0.1, 0, 0.1))
  p <- p + ggplot2::ggtitle("Significant metabolites by KEGG Super Pathway")

  #pick out at most top 3 significant metabolites in each class
  top3 <- foldChgV0PreV0Post %>% dplyr::filter(SIG==TRUE) %>% dplyr::group_by(SUPER_PATHWAY) %>% dplyr::top_n(3, M_LOG10P_VAL) %>% dplyr::ungroup()

  p <- p + ggrepel::geom_text_repel(data = top3, ggplot2::aes(label=as.character(BIOCHEMICAL)), size=5,
                                    box.padding = ggplot2::unit(0.2, "lines"), point.padding = ggplot2::unit(1, "lines"), force=3, max.iter = 2e4, nudge_y = 0.05, nudge_x = -0.05)
  p <- p + theme(axis.text=element_text(size=16),
                 axis.title=element_text(size=24), legend.text=element_text(size=24), legend.title=element_text(size=24))


  ggplot2::ggsave(filename = file.path(outputDir, "VolcanoV0PreV0Post_ColourByPathway.pdf"),  p, width = 16, height = 8)

  #lots of lipids, lets check by sub_pathway
  p <- ggplot2::ggplot(dplyr::filter(foldChgV0PreV0Post, SUPER_PATHWAY=="Lipid"), ggplot2::aes(x = LOG_2_FOLD,y = M_LOG10P_VAL)) + ggplot2::geom_point(ggplot2::aes(colour = SIG)) + ggplot2::facet_wrap(~SUB_PATHWAY , nrow = 5)
  ggplot2::ggsave(filename = file.path(outputDir, "VolcanoV0PreV0Post_LipidSubPathway.pdf"),  p, width = 12, height = 8)

  #save the fold change V0Pre-V0Post
  saveRDS(foldChgV0PreV0Post, file = file.path(outputDir, "V0Post-V0PreFC.RDS"))
}

#' Runs the DILT1D analyte analysis for comparing Ag1,5
#'
#' The method reads in the data, cleans it, and estimates
#' the correlations between CPeptide, Ag1,5, clinical glucose
#' and Self-Measured glucose in the DILT1D trial
#'
#' @param fileLocationsList A list of names file locations
#' @param outputDir A character filepath store the output from the analysis
#'
#' @export

RunAnalyteAnalysis <- function (fileLocationsList, outputDir){

  PROP_MISSING <- 0.2

  #define some file locations

  message(paste0("Creating / cleaning Output directory at ", outputDir))
  MetabolonR::CreateOrClean(outputDir)

  message("Loading and Normalising dataset")
  dilt1dData <-  LoadDILT1DData(fileLocationsList$sampleFile, fileLocationsList$metaboliteFile, fileLocationsList$metaboliteData)

  sampleInfo <- dilt1dData[[1]]
  metaboliteInfo <- dilt1dData[[2]]
  normMetabCounts <- dilt1dData[[3]]

  message(paste0("Removing metabolites missing in more than ", 100*PROP_MISSING , "% of the samples"))
  DILMetDataLongProp <- MetabolonR::FilterMissingMets(normMetabCounts, PROP_MISSING)
  DILMetDataLongProp <- MetabolonR::TransformLog10Mets(DILMetDataLongProp)

  message(paste0("Loading DILT1D covariate data from ", fileLocationsList$covariateStem, fileLocationsList$covariateDate))
  covariates <- LoadDILT1DCovariates(fileLocationsList$covariateStem, fileLocationsList$covariateDate, fileLocationsList$covariateSourceFile)

  #pull out Ag1,5 metabolite info
  metabolite15AG <- dplyr::select(dplyr::filter(metaboliteInfo, grepl("1,5-anhydroglucitol",BIOCHEMICAL)), METABOLITE_NAME)
  metaboliteGlucose <- dplyr::select(dplyr::filter(metaboliteInfo, grepl("glucose",BIOCHEMICAL)), METABOLITE_NAME)

  #metabolon Ag1,5 data for all patient visits merged with patient data
  data15AG <- dplyr::filter(DILMetDataLongProp, METABOLITE_NAME == metabolite15AG$METABOLITE_NAME)

  #metabolon glucose data for all VISITs merged with patient data
  dataGlucose <- dplyr::filter(DILMetDataLongProp, METABOLITE_NAME == metaboliteGlucose$METABOLITE_NAME)

  message("Loading clinical data")
  cpeptideData <- readDILT1DCpeptideData(fileLocationsList$cpeptideFile)
  hbaData <- readDILT1DHBAData(fileLocationsList$hbaFile)
  smbgData <- readDILT1DSMBGData(fileLocationsList$smbgFile)
  clinChemData <- readDILT1DClinchemData(fileLocationsList$clinChemFile)
  gluData <- fetchClinchemAnalyte(clinChemData, analyte = "GLUCOSE..MMOL.L")

  message("Plotting data availability...")

  #summarise what data we have about which visits
  #sampleInfo, covariates, cpeptideData, hbaData, agData, gluData, smbgData

  p <- PlotAnalyteAvailability(sampleInfo, covariates, cpeptideData, hbaData, data15AG, gluData, smbgData)
  ggsave(file.path(outputDir, "SummaryDataAvailability.pdf"), p, width = 10, height = 8)

  #we need to rename the cpeptide and hba screening visit to V0Pre to merge pre-IL-2 administration, both random samples.
  cpeptideData$VISIT <- plyr::revalue(cpeptideData$VISIT, c("Sc" = "V0Pre"))
  hbaData$VISIT <- plyr::revalue(hbaData$VISIT, c("Sc"="V0Pre"))

  message("Plotting baseline figures")
  message("...Ag1,5 panel")
  ag15plots <- PlotAg15Panel(sampleInfo, covariates, data15AG, gluData, hbaData, cpeptideData)
  ggsave(file.path(outputDir, "Ag15_Glucose_Correlation.pdf"), ag15plots$Glucose, width = 7, height = 7)
  ggsave(file.path(outputDir, "Ag15_HbA1c_Correlation.pdf"), ag15plots$HbA1c, width = 7, height = 7)
  ggsave(file.path(outputDir, "Ag15_CPeptide_Correlation.pdf"), ag15plots$cpep, width = 7, height = 7)

  message("...C-peptide panel")
  cpepPlots <- PlotCPepPanel(cpeptideData, hbaData, gluData)
  ggsave(file.path(outputDir, "CPeptide_Glucose_Correlation.pdf"), cpepPlots$Glucose, width = 7, height = 7)
  ggsave(file.path(outputDir, "CPeptide_HbA1c_Correlation.pdf"), cpepPlots[[1]], width = 7, height = 7)

  message("...Glucose panel")
  gluPlots <- PlotGlucosePanel(sampleInfo, covariates, dataGlucose, gluData, smbgData)

  ggsave(file.path(outputDir, "Metabolon_SMBG_Glucose_Correlation.pdf"), gluPlots$MetabolonSMBG, width = 7, height = 7)
  ggsave(file.path(outputDir, "Clinical_SMBG_Glucose_Correlation.pdf"), gluPlots$ClinicalSMBG, width = 7, height = 7)
  ggsave(file.path(outputDir, "Clinical_SMBG_Glucose_Correlation_V10.pdf"), gluPlots$ClinicalSMBGV10, width = 7, height = 7)
  ggsave(file.path(outputDir, "Metabolon_Clinical_Glucose_Correlation_V10.pdf"), gluPlots$MetabolonClinical, width = 7, height = 7)

}
