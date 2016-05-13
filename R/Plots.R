
#' Plots the covariates of sex, age and dose in DILT1D study
#'
#' @param dilt1dCovariates a Data Frame of the DILT1D covariates from the doses file
#' @return A ggplot of the covariates faceted by sex
#'
#'
#' @importFrom ggplot2 ggplot aes theme scale_x_continuous scale_y_continuous geom_point facet_grid scale_size_discrete geom_jitter element_text
#' @export

PlotCovariates <- function(dilt1dCovariates) {
  #### Trial information ####

  partition <- dilt1dCovariates[,c("sex","age_V0","strategyNEW")]

  p <- ggplot(partition,aes(x=age_V0,y=strategyNEW,colour=sex)) + geom_point() + geom_jitter()
  p <- ggplot(partition,aes(x=age_V0,y=sex,colour=strategyNEW)) + geom_point()
  partition.freq <- plyr::count(partition, c("age_V0","strategyNEW","sex"))
  partition.freq$freq <- factor(partition.freq$freq)
  colnames(partition.freq)[colnames(partition.freq) == "freq"] <- "num_of_individuals"

  #partition.freq takes into account participants with the same sex, age and dose

  p <- ggplot(partition.freq,aes(x=age_V0,y=strategyNEW)) + geom_point(aes(size = num_of_individuals,colour=sex),alpha=0.5) + scale_size_discrete(range = c(3, 6)) + facet_grid(.~sex)
  p <- p + scale_x_continuous(name="Age at baseline") + scale_y_continuous(name = "dose of IL-2, IU/m^2")
  p <- p + theme(axis.title.x = element_text(face="bold", size=16),
                 axis.text.x  = element_text(size=10))
  p <- p + theme(axis.title.y = element_text(face="bold", size=16),
                 axis.text.y  = element_text(size=10))
  p <- p + theme(legend.text = element_text( size = 10, face = "bold"), legend.title = element_text(size=16, face="bold"), legend.position="bottom")
  return(p)
}



#' Plots a presentation ready plot comparing the effect sizes of the Krumsiek and DILT1D
#' sex lined metabolites
#'
#' @param sexRegressionResults A Data Frame containing cleaned information about the samples
#' @return A ggplot volcano plot containing the log2 fold change vs -log10P value of the fold changes across Visits
#' in the DILT1D data trial
#' @examples
#' p <- PlotSexMetaboliteComparison(sexRegressionResults)
#'
#' @export

PlotVolcano <- function(DFregression, regressor, cMetInfo, plotLabels = TRUE, plotColour = TRUE, correction = "fdr"){

  dfFoldChange <- merge(dplyr::select(DFregression[[1]],METABOLITE_NAME, starts_with(regressor)),dplyr::select(DFregression[[4]],METABOLITE_NAME, starts_with(regressor)),by="METABOLITE_NAME")
  dfFoldChange <- merge(dfFoldChange , dplyr::select(cMetInfo, METABOLITE_NAME, METABOLITE_TYPE), by = "METABOLITE_NAME")
  colnames(dfFoldChange) <- c("METABOLITE_NAME","COEFFICIENT","P_VAL","METABOLITE_TYPE")
  dfFoldChange$M_LOG10P_VAL <- -log10(dfFoldChange$P_VAL)
  dfFoldChange$P_VAL_CORRECTED <- p.adjust(dfFoldChange$P_VAL,correction)
  dfFoldChange$SIG <- dfFoldChange$P_VAL_CORRECTED < 0.05
  dfFoldChange$LOG_2_FOLD <- log2(10^dfFoldChange$COEFFICIENT)
  dfFoldChange$PATHWAY_SORTORDER <- cMetInfo$PATHWAY_SORTORDER[match(dfFoldChange$METABOLITE_NAME, cMetInfo$METABOLITE_NAME)]
  dfFoldChange$BIOCHEMICAL <- cMetInfo$BIOCHEMICAL[match(dfFoldChange$METABOLITE_NAME, cMetInfo$METABOLITE_NAME)]

  #for the labelling...
  dfFoldChange <- dfFoldChange %>% dplyr::mutate(METABOLITE_TYPE=ifelse(METABOLITE_TYPE == "UNKNOWN" | METABOLITE_TYPE == "KNOWN","BIOLOGICAL","IS"))

  sigMetNames <- dplyr::filter(dfFoldChange,SIG==T)

  p <- ggplot2::ggplot(dfFoldChange,ggplot2::aes(x=LOG_2_FOLD,y=M_LOG10P_VAL))

  if(plotLabels) {
    p  <- p +  ggplot2::geom_text(ggplot2::aes(label=ifelse(SIG==T, as.character(BIOCHEMICAL),'')),hjust=0,vjust=0,size=3)
  }

  p <- p + ggplot2::labs(x="Log2 Fold Change in Metabolite Expression",y="-log10 p-value (adjusted)")

  p1 <- p + ggplot2::geom_point(ggplot2::aes(colour = METABOLITE_TYPE, order=as.factor(METABOLITE_TYPE)))
  p2 <- p
  if(plotColour) {
    p2 <- p + ggplot2::geom_point(ggplot2::aes(colour = SIG),size=0.5)
  }
  else {
    p2 <- p + ggplot2::geom_point(size=0.5)
  }
  return (list(p1,p2,dfFoldChange))
}


#' Plots a presentation ready plot comparing the effect sizes of the Krumsiek and DILT1D
#' sex lined metabolites
#'
#' @param dilt1dCoeffs A Data Frame containing coefficents for the sex metabolite regressions
#' @param dilt1dCoeffs A Data Frame containing p-values for the sex metabolite regressions
#' @param dilt1dCoeffs A Data Frame literature results from Krumsiek sex metabolite regressions
#' @param metaboliteInfo A Data Frame containing the metabolite metadata
#'
#'
#' @return A 2 element List. A ggplot containing the effect sizes of the sex metabolites in Krumsiek and DILT1D
#' and the data frame used for the plotting
#'
#' @examples
#' p <- PlotSexMetaboliteComparison(dilt1dCoeffs, dilt1dPVals, krumResults)
#'
#' @export

PlotSexMetaboliteComparison <- function(dilt1dCoeffs, dilt1dPVals, krumResults, metaboliteInfo) {

  krumsiek <- dplyr::select(krumResults, METABOLITE_NAME, KRUM_FOLD=log2.fold.change, KRUM_PVAL=p.value..gender.)
  dilt1d <- dplyr::select(dilt1dCoeffs, DIL_FOLD=contains("factor(sex)F"), METABOLITE_NAME)

  oursVsKrum <- merge(krumsiek, dilt1d, by = "METABOLITE_NAME")
  oursVsKrum <- merge(oursVsKrum,dplyr::select(metaboliteInfo,c(METABOLITE_NAME,PATHWAY_SORTORDER,BIOCHEMICAL)))
  oursVsKrum <- dplyr::filter(oursVsKrum,!grepl('^X',PATHWAY_SORTORDER))
  oursVsKrum$ANDRO <- grepl('andro', oursVsKrum$BIOCHEMICAL)
  oursVsKrum$DIL_FOLD <- log2(10^(oursVsKrum$DIL_FOLD))
  oursVsKrum$M_LOG10_P <- -log10(oursVsKrum$KRUM_PVAL)

  p <- ggplot2::ggplot(oursVsKrum, ggplot2::aes(x = KRUM_FOLD, y = DIL_FOLD, colour=ANDRO)) + ggplot2::geom_point()
  p <- p + ggplot2::xlim(-2,1) + ggplot2::ylim(-2,1)
  p <- p + ggplot2::labs(x="Log2 Fold Change in Krumsiek Sex Metabolites", y="Log2 Fold Change in DILT1D Sex Metabolites")
  p <- p + ggplot2::labs(title="Replication in fold changes", colour="androgens")
  return (list(p, oursVsKrum))

}

#' Plots a presentation ready plot representing the number of sampled analytes available at each time
#' point in DILT1D
#'
#' @param sampleInfo A Data Frame containing sample data for DILT1D Metabolon data
#' @param covariateData A Data Frame containing covariate data for DILT1D
#' @param cpeptideData A Data Frame containing C-Peptide expression data
#' @param hbaData A Data Frame containing HbA1c expression data
#' @param agData A Data Frame containing Metabolon Ag1,5 expression data
#' @param gluData A Data Frame containing clinical glucose measurements
#' @param smbgData A Data Frame containing clinical smbg measurements
#'
#' @return A ggplot containing a tiled plot of the data availability
#' @examples
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient geom_text theme labs ggsave
#' @importFrom magrittr %>%


PlotAnalyteAvailability <- function(sampleInfo, covariates, cpeptideData, hbaData, agData, gluData, smbgData) {

  covariatesAnd15AG <- dplyr::left_join(agData, dplyr::select(sampleInfo, c(SAMPLE_NAME, SUBJECT_ID, VISIT)), by="SAMPLE_NAME")
  covariatesAnd15AG <- dplyr::left_join(covariatesAnd15AG, dplyr::select(covariates, c(SUBJECT_ID,sex, dose_level, age_V0, strategyNEW, strategyNEW.colours, strategyNEW.labels)), by="SUBJECT_ID")

  a <- plyr::ddply(cpeptideData,c("VISIT"), plyr::summarise, SUBJECT_COUNT_PER_VISIT = length(SUBJECT_ID),DATA_TYPE = "CPeptide")
  b <- plyr::ddply(hbaData,c("VISIT"), plyr::summarise, SUBJECT_COUNT_PER_VISIT = length(SUBJECT_ID),DATA_TYPE = "HBA")
  c <- plyr::ddply(covariatesAnd15AG,c("VISIT"), plyr::summarise, SUBJECT_COUNT_PER_VISIT = length(SUBJECT_ID),DATA_TYPE = "AG15")
  d <- plyr::ddply(gluData,c("VISIT") ,plyr::summarise, SUBJECT_COUNT_PER_VISIT = length(SUBJECT_ID),DATA_TYPE = "BLOOD_GLUCOSE")
  e <- plyr::ddply(smbgData,c("VISIT"), plyr::summarise, SUBJECT_COUNT_PER_VISIT = length(SUBJECT_ID),DATA_TYPE = "SELF_MONITORED_BLOOD_GLUCOSE")

  summaryVisits <- rbind(a,b,c,d,e)
  summaryVisits$VISIT <- factor(summaryVisits$VISIT, levels = c('Sc', 'V0Pre', 'V0Post', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10'))
  summaryVisits <- reshape2::dcast(summaryVisits,VISIT ~ DATA_TYPE,value.var = "SUBJECT_COUNT_PER_VISIT")
  summaryVisits[is.na(summaryVisits)] <- 0
  summaryVisits <- tidyr::gather(summaryVisits, DATA_TYPE, SUBJECT_COUNT_PER_VISIT, AG15:SELF_MONITORED_BLOOD_GLUCOSE)

  #tidy up biochemical names for plotting
  summaryVisits$DATA_TYPE <- gsub(x = summaryVisits$DATA_TYPE,"SELF_MONITORED_BLOOD_GLUCOSE","SMBG")
  summaryVisits$DATA_TYPE <- gsub(x = summaryVisits$DATA_TYPE,"CPeptide","C-Peptide")
  summaryVisits$DATA_TYPE <- gsub(x = summaryVisits$DATA_TYPE,"AG15","Ag1,5")
  summaryVisits$DATA_TYPE <- gsub(x = summaryVisits$DATA_TYPE,"BLOOD_GLUCOSE","Clinical Glucose")
  summaryVisits$DATA_TYPE <- gsub(x = summaryVisits$DATA_TYPE,"HBA","HbA1c")

  p <- ggplot(summaryVisits, aes(VISIT, DATA_TYPE)) + geom_tile(aes(fill = SUBJECT_COUNT_PER_VISIT))
  p <- p + scale_fill_gradient(name = "Number of Observations", low = "white", high = "steelblue") + geom_text(aes(label=SUBJECT_COUNT_PER_VISIT))
  p <- p + labs(x="VISIT", y="BIOCHEMICAL",title="Datapoints for DILT1D subjects") + theme(text = ggplot2::element_text(size=20))
  p <- p + theme(legend.position = "bottom", legend.title = element_text(size=16, face="bold"))

  return(p)
}

#' Plots a presentation ready plot representing the correlations between Ag1,5 and
#' point in DILT1D
#'
#' @param sampleInfo A Data Frame containing sample data for DILT1D Metabolon data
#' @param covariateData A Data Frame containing covariate data for DILT1D
#' @param data15AG A Data Frame containing Metabolon Ag1,5 expression data
#' @param gluData A Data Frame containing clinical glucose measurements for DILT1D patients
#' @param hbaData A Data Frame containing clinical HbA1c measurements for DILT1D patients
#' @param cpeptideData A Data Frame containing clinical C-peptide measurements for DILT1D patients
#'
#' @return A list of ggplots containing correlations between Ag1,5 and the analytes
#' @examples
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point theme labs ggsave
#' @importFrom magrittr %>%

PlotAg15Panel <- function(sampleInfo, covariateData, data15AG, gluData, hbaData, cpeptideData) {
  plottingSet <- dplyr::left_join(data15AG, dplyr::select(sampleInfo, c(SAMPLE_NAME, SUBJECT_ID, VISIT)), by="SAMPLE_NAME")
  plottingSet <- dplyr::left_join(plottingSet, dplyr::select(covariateData, c(SUBJECT_ID, sex, dose_level, age_V0, strategyNEW, strategyNEW.colours, strategyNEW.labels)), by="SUBJECT_ID")
  gluPlotData <- suppressWarnings(dplyr::left_join(plottingSet, gluData, by=c("SUBJECT_ID","VISIT")))

  p <- (ggplot(dplyr::filter(gluPlotData, VISIT=="V0Pre"),aes(x = CLINCHEM_VALUE, y = METABOLITE_COUNT)) + geom_point() +
        labs(x="Clinical Glucose, mmol/L",y="Log10(1,5 Ag), Raw Count",title="glucose vs 1,5Ag"))

  hbaPlotData <- suppressWarnings(dplyr::left_join(plottingSet, hbaData, by=c("SUBJECT_ID","VISIT")))
  q <- (ggplot(dplyr::filter(hbaPlotData, VISIT=="V0Pre"),aes(x = CLINICAL.HBA1C, y = METABOLITE_COUNT)) + geom_point() +
    labs(x="Clinical HbA1c mmol/mol, Screening Visit",y="Log10(1,5 Ag), Raw ion-count",title="Clinical HbA1c vs 1,5Ag"))

  cpepPlotData <- suppressWarnings(dplyr::left_join(plottingSet, cpeptideData, by=c("SUBJECT_ID","VISIT")))
  r <- (ggplot(dplyr::filter(cpepPlotData, VISIT=="V0Pre"),aes(x = C_PEPTIDE.PMOL.L, y = METABOLITE_COUNT)) + geom_point() +
         labs(x="Clinical C-peptide pmol/L, Screening Visit",y="Log10(1,5 Ag), Raw ion-count",title="Clinical C-peptide vs 1,5Ag"))

  analytePlots <- list(p, q, r)
  names(analytePlots) <- c("Glucose", "HbA1c", "cpep")

  return(analytePlots)

}

#' Plots a presentation ready plot representing the correlations between C-peptide
#' and HbA1c and Clinical Glucose measurements
#'
#' @param gluData A Data Frame containing clinical glucose measurements for DILT1D patients
#' @param hbaData A Data Frame containing clinical HbA1c measurements for DILT1D patients
#' @param cpeptideData A Data Frame containing clinical C-peptide measurements for DILT1D patients
#'
#' @return A list of ggplots containing correlations between C-peptide and the analytes
#' @examples
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap theme labs ggsave ggtitle xlab ylab

PlotCPepPanel <- function(cpeptideData, hbaData, gluData) {

  dat <- dplyr::left_join(dplyr::filter(cpeptideData, VISIT %in% c("V0Pre","V10")), dplyr::filter(hbaData, VISIT %in% c("V0Pre","V10")), by=c("SUBJECT_ID","VISIT"))
  p <- ggplot(dat, aes(x = C_PEPTIDE.PMOL.L, y = CLINICAL.HBA1C)) + geom_point() + xlab("Cpeptide pM/M") + ylab("hba1c")
  p <- p + ggtitle("C-peptide versus HbA1c: Visits V0Pre and V10") + facet_wrap(~VISIT)

  dat <- dplyr::left_join(dplyr::filter(cpeptideData, VISIT %in% c("V0Pre","V10")), dplyr::filter(gluData, VISIT %in% c("V0Pre","V10")),by=c("SUBJECT_ID","VISIT"))
  q <- ggplot(dat,aes(x = C_PEPTIDE.PMOL.L,y = CLINCHEM_VALUE)) + geom_point() + xlab("Cpeptide pM/M") + ylab("Clinical Glucose, mmol/L")
  q <- q + ggtitle("C-peptide versus Clinical Glucose: Visits V0Pre and V10") + facet_wrap(~VISIT)

  cpepPlots <- list(p, q)
  names(cpepPlots) <- c("Hba1c", "Glucose")
  return(cpepPlots)
}

#' Plots a presentation ready plot representing the correlations between Ag1,5 and
#' point in DILT1D
#'
#' @param sampleInfo A Data Frame containing sample data for DILT1D Metabolon data
#' @param covariateData A Data Frame containing covariate data for DILT1D
#' @param metabolonGlu A Data Frame containing Metabolon glucose expression data
#' @param gluData A Data Frame containing clinical glucose measurements for DILT1D patients
#' @param smbgData A Data Frame containing SMBG measurements for DILT1D patients
#'
#' @return A list of ggplots containing correlations between the various glucose measurements
#' @examples
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point theme labs ggsave ggtitle xlab ylab
#' @importFrom magrittr %>%

PlotGlucosePanel <- function(sampleInfo, covariateData, metabolonGlu, gluData, smbgData) {

  plottingSet <- dplyr::left_join(metabolonGlu, dplyr::select(sampleInfo, c(SAMPLE_NAME, SUBJECT_ID, VISIT)), by="SAMPLE_NAME")
  plottingSet <- dplyr::left_join(plottingSet, dplyr::select(covariateData, c(SUBJECT_ID, sex, dose_level, age_V0, strategyNEW, strategyNEW.colours, strategyNEW.labels)), by="SUBJECT_ID")

  # Baseline Metabolon glucose vs SMBG at V0Pre
  dat <- dplyr::left_join(dplyr::filter(plottingSet, VISIT %in% c("V0Pre")), dplyr::filter(smbgData, VISIT %in% c("V0Pre")), by=c("SUBJECT_ID","VISIT"))
  p <- ggplot(dat, aes(x=METABOLITE_COUNT, y=SMBG)) + geom_point() + xlab("Metabolon log10 ion-count") + ylab("SMBG, mmol/L")
  p <- p + ggtitle("Metabolon glucose versus SMBG (V0Pre)")

  #clinical glucose vs SMBG V0Pre and V10
  dat <- dplyr::left_join(dplyr::filter(gluData,VISIT %in% c("V0Pre")), dplyr::filter(smbgData,VISIT %in% c("V0Pre","V10")), by=c("SUBJECT_ID","VISIT"))
  q <- ggplot(dat, aes(x=CLINCHEM_VALUE, y=SMBG)) + geom_point() + xlab("Clinical Glucose, mmol/L") + ylab("SMBG, mmol/L")
  q <- q + ggtitle("Clinical glucose versus SMBG (V0Pre)")

  dat <- dplyr::left_join(dplyr::filter(gluData, VISIT %in% c("V10")), dplyr::filter(smbgData, VISIT %in% c("V0Pre","V10")), by=c("SUBJECT_ID","VISIT"))
  r <- ggplot(dat, aes(x=CLINCHEM_VALUE, y=SMBG)) + geom_point() + xlab("Clinical Glucose, mmol/L") + ylab("SMBG, mmol/L")
  r <- r + ggtitle("Clinical glucose versus SMBG (V10)")

  #metabolon glucose vs clinical glucose at V0Pre
  dat <- dplyr::left_join(dplyr::filter(plottingSet, VISIT %in% c("V0Pre")), dplyr::filter(gluData, VISIT %in% c("V0Pre")),by=c("SUBJECT_ID","VISIT"))
  s <- ggplot(dat, aes(x=METABOLITE_COUNT, y=CLINCHEM_VALUE)) + geom_point() + xlab("Metabolon log10 ion-count") + ylab("Clinical Glucose, mmol/L")
  s <- s + ggtitle("Metabolon glucose versus Clinical (V0Pre)")

  glucosePlots <- list(p, q, r, s)
  names(glucosePlots) <- c("MetabolonSMBG", "ClinicalSMBG", "ClinicalSMBGV10", "MetabolonClinical")
  return(glucosePlots)
}

#' Plots a presentation ready plot representing the variance components of the
#' metabolite data in DILT1D
#'
#' @param ICCs A Data Frame containing variance components data from the Linear Mixed Model
#'
#' @return A list of ggplots containing plots of the technical, within subject and between subject variation
#' @examples
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point theme labs ggsave ggtitle xlab ylab element_blank
#' @importFrom magrittr %>%

PlotMetaboliteVarianceComponents <- function(ICCs) {

  ICCs <- dplyr::arrange(ICCs,ICC)
  ICCs$metnamesOrdered <- factor(ICCs$METABOLITE_NAME, levels = ICCs$METABOLITE_NAME)
  p1 <- ggplot(ICCs,aes(x=metnamesOrdered,y=ICC)) + geom_point(aes(colour=SUPER_PATHWAY))+ theme(axis.text.x=element_blank(),
                                                                                                axis.ticks=element_blank())
  lMedianICC <- floor(median(seq(1,nrow(ICCs))))
  uMedianICC <- ceiling(median(seq(1,nrow(ICCs))))
  medianICC <- (ICCs[lMedianICC,]$ICC + ICCs[uMedianICC,]$ICC) / 2
  p1 <- p1 + ggtitle(paste("Median individual (within and total) variation as % of total individual variation =",format(medianICC*100, digits = 3 ),"%"))


  ICCs <- dplyr::arrange(ICCs,phi_B_BW)
  ICCs$metnamesOrdered <- factor(ICCs$METABOLITE_NAME, levels = ICCs$METABOLITE_NAME)
  p2 <- ggplot(ICCs,aes(x=metnamesOrdered,y=phi_B_BW)) + geom_point(aes(colour=SUPER_PATHWAY))+ theme(axis.text.x=element_blank(),
                                                                                                     axis.ticks=element_blank())
  medianB_BW <- (ICCs[lMedianICC,]$phi_B_BW + ICCs[uMedianICC,]$phi_B_BW) /2
  p2 <- p2 + ggtitle(paste("Median between individual variation as % of total individual variation =",format(medianB_BW*100, digits = 3 ),"%"))


  ICCs <- dplyr::arrange(ICCs,phi_B_T)
  ICCs$metnamesOrdered <- factor(ICCs$METABOLITE_NAME, levels = ICCs$METABOLITE_NAME)
  p3 <- ggplot(ICCs, aes(x=metnamesOrdered,y=phi_B_T)) + geom_point(aes(colour=SUPER_PATHWAY))+ theme(axis.text.x=element_blank(),
                                                                                                    axis.ticks=element_blank())
  medianB_T <- (ICCs[lMedianICC,]$phi_B_T + ICCs[uMedianICC,]$phi_B_T) /2
  p3 <- p3 + ggtitle(paste("Median between individual variation as % of total =", format(medianB_T*100, digits = 3 ),"%"))


  ICCs$residual <- 1-ICCs$ICC
  ICCs <- dplyr::arrange(ICCs, residual)
  ICCs$metnamesOrdered <- factor(ICCs$METABOLITE_NAME, levels = ICCs$METABOLITE_NAME)
  p4 <- ggplot(ICCs,aes(x=metnamesOrdered,y=residual)) + geom_point(aes(colour=SUPER_PATHWAY))

  medianResidual <- (ICCs[lMedianICC,]$residual + ICCs[uMedianICC,]$residual) /2
  p4 <- p4 + ggtitle(paste("Median technical variation =",format(medianResidual * 100, digits = 3 ),"%"))

  plotList <- list(p1, p2, p3, p4)
  names(plotList) <- c("MedianIndividual", "MedianBetweenIndividualOfTotalIndividual", "MedianBetweenIndividualOfTotal", "MedianTechnicalVariation")
  return(plotList)

}
