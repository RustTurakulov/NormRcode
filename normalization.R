library(data.table)
library(parallel)
library(tidyr)
library(openxlsx)
library(meffil)
library(dplyr)
# Bunch for visualization/dimension reduction 
library(tibble)
library(uwot)
library(ggplot2)
library(Rtsne)
library(gridExtra)
library(tidyverse)
library(densvis)

# parsing folder locations  
args = commandArgs(trailingOnly=TRUE)
samplesheet <- args[1];
batchdir    <- sub("(_20[0-9][0-9]-[01][0-9]-[0-9][0-9].xlsx)", "", samplesheet, perl=T);
batchdir    <- paste0("TRANSFER/SAMPLESHEETS/", batchdir);
samplesheet <- paste0("TRANSFER/SAMPLESHEETS/", samplesheet);
dir.create(batchdir);

#create sample sheet
#anno_base <- openxlsx::read.xlsx("NormRcode/Sample_sheet_test.xlsx")
anno_base <- openxlsx::read.xlsx(samplesheet)

#filter cases based on string 
#(in this case, we are interested in all non-duplicate cases)
#use 'unique(anno_base$panCNS_study)' to see choices
anno <- anno_base %>% filter(panCNS_study=="Case" |
                               panCNS_study=="Recurrent_nonmatched" |
                               panCNS_study=="Recurrent_matched" |
                               panCNS_study=="Metastasis_matched" |
                               panCNS_study=="Metastasis_nonmatched" |
                               panCNS_study=="Non_neoplastic_bulk" |
                               panCNS_study=="Case_nonCNS")

#other filtering options
#anno <- anno_base[!is.na(anno_base$Purity_known_class) | grepl("CELLS",anno_base$Purity_cell),]
#anno <- anno_base %>% filter(!is.na(Purity_study))

#format sample sheet for meffil
anno <- anno[,c( "idat_filename", "Sex", "Sentrix_ID", "Sentrix_position", "Basename_Ubuntu","Platform_methy","material_prediction")]
anno <- separate(anno, Sentrix_position, into = c("sentrix_row", "sentrix_col"), sep = 3, remove = TRUE)
anno$sentrix_row = gsub("\\R","",anno$sentrix_row)
anno$sentrix_col = gsub("\\C","",anno$sentrix_col)
names(anno) <- c("Sample_Name", "Sex", "Slide", "sentrix_row", "sentrix_col", "Basename", "Array", "Material")
anno$Slide <- as.numeric(as.character(anno$Slide))
anno <- anno %>% filter(!is.na(Basename))
anno <- anno %>% filter(!duplicated(Basename))
#anno <- anno %>% filter(Array=="HumanMethylationEPIC")

# #for copy number only
# anno_450k <- anno %>% filter(Array == "HumanMethylation450")
# anno_EPIC <- anno %>% filter(Array == "HumanMethylationEPIC")

#define number of cores for parallelization (linux-based machines)
#code for parallel implcd ..ementation on windows-based machines is different
#options(mc.cores=26)
options(mc.cores=1)

cores <- options()$mc.cores


# #list of available chips and feature sets for meffil
# meffil.list.chips()
# meffil.list.featuresets()
# 
# #list of available cell type references
# meffil.list.cell.type.references()

#qc object: perform background correction, dye bias correction, sex prediction and cell count estimates
qc.objects <- meffil.qc(anno, 
                        featureset = "common", 
                        cell.type.reference = NA, 
                        verbose = T)

#in any bioinformatic process, it is good practice to append saved files with
#details of the sample numbers, cohort name, functions applied, etc.
saveRDS(qc.objects, file = file.path(paste0(batchdir,"/qc.objects.rds")))
qc.objects <- readRDS(paste0(batchdir,"/qc.objects.rds"))
gc()

# #use this code if error results from 'meffil.qc'
# #in my experience, this occurs in two situations: missing idat file(s) in folder
# #and incompatible idat file
# qc.objects_error <- qc.objects[which(sapply(qc.objects, class) == 'try-error')]
# qc.objects_error <- as.data.frame(do.call(c, qc.objects_error))
# write.table(qc.objects_error, "meffil_output/error.txt", sep = "\t")

# #Obtain the matrix of genotypes for comparison with those measured on the microarray
# annotation <- qc.objects[[1]]$featureset
# writeLines(meffil.snp.names(annotation), con="snp-names.txt")
# #command shell > ./plink2 --bfile dataset --extract snp-names.txt --recode A --out genotypes
# filenames <- "genotypes.raw"
# genotypes <- meffil.extract.genotypes(filenames)

#QC analysis of the raw data
qc.summary <- meffil.qc.summary(qc.objects)
saveRDS(qc.summary, file = file.path(paste0(batchdir,"/qc.summary.rds")))
qc.summary <- readRDS(paste0(batchdir,"/qc.summary.rds"))

#extract bad cpgs from qc.summary
bad.cpgs <- qc.summary$bad.cpgs$name
saveRDS(bad.cpgs, file = file.path(paste0(batchdir,"/bad.cpgs.rds")))
bad.cpgs <- readRDS(paste0(batchdir,"/bad.cpgs.rds"))

#create html-based report of qc results
meffil.qc.report(qc.summary, output.file=paste0(batchdir,"/qc_report.html"))


#this section outlines what I consider to be the important qc checks for sample quality
#outlier samples whose predicted median methylated signal is more than 3 standard deviations from the expected
bad_samples_meth_unmeth <- as.data.frame(qc.summary$meth.unmeth.summary$tab)
bad_samples_meth_unmeth <- bad_samples_meth_unmeth[bad_samples_meth_unmeth$outliers==TRUE,]
anno_bad_meth_unmeth <- anno_base[anno_base$idat_filename %in% bad_samples_meth_unmeth$sample.name,]

# #remove specific samples
# bad_samples_meth_unmeth <- bad_samples_meth_unmeth[!bad_samples_meth_unmeth$sample.name=="GSM2403462_8622007027_R02C01",]

#deviations from mean values for control probes
bad_samples_control <- as.data.frame(qc.summary$controlmeans.summary$tab)
bad_samples_control <- bad_samples_control[bad_samples_control$outliers==TRUE,]
bad_samples_control <- bad_samples_control[bad_samples_control$variable=="bisulfite1" | bad_samples_control$variable=="bisulfite2",]
anno_bad_control <- anno_base[anno_base$idat_filename %in% bad_samples_control$sample.name,]

#samples with proportion of probes with detection p-value > 0.01 is > 0.2
bad_samples_detectionp <- as.data.frame(qc.summary$sample.detectionp.summary$tab)
bad_samples_detectionp <- bad_samples_detectionp[bad_samples_detectionp$outliers==TRUE,]
anno_bad_detectionp <- anno_base[anno_base$idat_filename %in% bad_samples_detectionp$sample.name,]

#smaples with proportion of probes with bead number < 3 is > 0.2
bad_samples_lowbead <- as.data.frame(qc.summary$sample.beadnum.summary$tab)
bad_samples_lowbead <- bad_samples_lowbead[bad_samples_lowbead$outliers==TRUE,]
anno_bad_lowbead <- anno_base[anno_base$idat_filename %in% bad_samples_lowbead$sample.name,]

#combine all bad samples into vector
bad_samples <- c(
  bad_samples_meth_unmeth$sample.name,
  bad_samples_control$sample.name,
  bad_samples_detectionp$sample.name,
  bad_samples_lowbead$sample.name)

#remove bad samples prior to performing quantile normalization
qc.objects <- meffil.remove.samples(qc.objects, bad_samples)

#determine the number of principal components of the control matrix to include in the quantile normalization
#plots the quantile residuals remaining after fitting different numbers of control matrix principal components
#I use the 'elbow' to determine the optimal number of pcs
pc.fit <- meffil.plot.pc.fit(qc.objects)
setwd(batchdir)
print(pc.fit$plot)
setwd("/data/MDATA")
pc <- 8
gc()

#remove control probe variance from the sample quantiles
#additional fixed and random effects can be included
norm.objects <- meffil.normalize.quantiles(qc.objects,
                                           fixed.effects = "Material",
                                           random.effects= "Slide", 
                                           number.pcs=pc,
                                           verbose = TRUE)

saveRDS(norm.objects, file = file.path(paste0(batchdir,"/norm.objects_pc8_fixedeffectsMaterial.rds")))
norm.objects <- readRDS(paste0(batchdir,"/norm.objects_pc8_fixedeffectsMaterial.rds"))
gc()

#remove objects not needed for further analysis (i.e. free up memory)
rm(list=setdiff(ls(), c("norm.objects", "bad.cpgs", "cores", "batchdir", "samplesheet")))
gc()

#normalize samples using their normalized quality control objects and remove bad CpGs (from the QC analysis)
norm.beta <- meffil.normalize.samples(norm.objects,
                                         just.beta = T,
                                         cpglist.remove = bad.cpgs,
                                         verbose = TRUE)

saveRDS(norm.beta, file = file.path(paste0(batchdir,"/norm.betas.rds"))); ## Original saveRDS(norm.signals, file = file.path("meffil_output/norm.signals.rds")) 
norm.beta <- readRDS(paste0(batchdir,"/norm.betas.rds"))

#convert normalized signals to beta values (if 'just.beta' = F)
#beta <- meffil.get.beta(norm.signals$M, norm.signals$U, pseudo = 100)
beta <- as.data.frame(norm.beta)
rm(norm.beta)
gc()

fwrite(beta, paste0(batchdir,"/betas_pc8_fixedeffectsMaterial.txt"), 
       row.names=TRUE, sep = "\t")

## summary report of the normalization performance
#pcs <- meffil.methylation.pcs(as.matrix(beta))
#norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs)


#####################################################
######## Make dimension reduction with UWOT  ########
anno <- openxlsx::read.xlsx(samplesheet)

PC <- prcomp(t(beta), center = TRUE, scale = FALSE) 
saveRDS(PC, paste0(batchdir,"/pc_prcomp_puritycorrected.rds"))
PC <- readRDS(paste0(batchdir,"/pc_prcomp_puritycorrected.rds"))

#unsupervised dimensionality reduction with UMAP
umap <- uwot::umap(PC$x,
                   #n_components = 3,
                   #pca = 100,
                   n_neighbors = 5,
                   #y = anno$Combined_class_match,
                   spread = 2,
                   min_dist = 0.1,
                   local_connectivity = 1,
                   bandwidth = 1);

umap <- as.data.frame(umap)
rownames(umap) <- rownames(PC$x)
umap <- umap %>% rownames_to_column("idat_filename")
anno <- anno[anno$idat_filename %in% umap$idat_filename,]

outf <- merge(umap, anno[,c(
    "idat_filename",
	"Sample",
	"Primary_study",
	"material_prediction",
	"Platform_methy",
	"Center_methy",
	"Primary_study",
	"material_prediction",
	"Platform_methy",
	"Center_methy",
	"Age",
	"Sex",
	"Sex_prediction",
	"OS_months",
	"OS_status",
	"Primary_category",
	"Histology",
	"Molecular",
	"Combined_class_match",
	"Secondary_class_match",
	"C19MC_segments",
	"fusion_matches",
	"breakpoint_genes",
	"amplifications",
	"deletions",
	"Likely_integrated_diagnosis")],
    by = "idat_filename");

fwrite(outf, paste0(batchdir,"/umap_1.txt"), 
    row.names=TRUE, sep = "\t");

file.rename(samplesheet, gsub(".xlsx",  ".done.xlsx", samplesheet));
message(Sys.time());
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
message(" End of normalization and UMAP script ");
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
