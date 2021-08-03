## Drew's normalization pipeline with small tweaks by Rust
## August 2021 :: pancancer set

args = commandArgs(trailingOnly = TRUE)
# args = "Placenta:NewCNSdir:TRANSFER/newbatch.txt"
if(length(args) != 1){
      stop("!!! Crash !!!\nParameter settings missed in sbatch: [cancer:outputfolder:newsamplesheet] \nWhere:\n   cancer \t  -- [Primary_category] from annotation table column [CC]\n   outputfolder   -- directory to save at /data/MDATA/TRANSFER\n   newsamplesheet -- full path to the new batch metadata\nTry to rerun sbatch with appropriate settings for like this:\n\nsbatch /data/MDATA/NormRcode/norm.sh  Thoracic|Bone and soft tissue:NewCNSdir:TRANSFER/newbatch.txt\n\n")
}else{
     parameter = as.character(args[1])
     parameter = unlist(strsplit(parameter, ":"))
     message("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", Sys.time(),"  Started with parameters:\nCancer:     \t", parameter[1], "\nOutputdir:\t", parameter[2], "\nSamples:\t", parameter[3], "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
}

options(scipen = 999)
library(data.table)
library(parallel)
library(tidyr)
library("openxlsx") #, lib.loc="/home/prattdw/R/4.0/library/")
library(meffil)
library(dplyr)
library(Rfast)
library(tibble)
library(readxl)
library(densvis)
setwd("/data/MDATA")
options(mc.cores=70)
cores <- options()$mc.cores
gc()

### Set paths around
idatrootfolder = "pancancer/iScan_raw";
batchdir       = "pancancer/NORMALIZATION"; 
batchdirout    = file.path("TRANSFER", parameter[2]);    # working folder
dir.create(batchdirout, recursive = TRUE)
samplesheet    = "bams_RNAseq/Sample_sheet_master.xlsm"; # the database
newsamplesheet = parameter[3];


#############  Inject static samplesheet
### Part 0 ##  anno_base -- rich file for the ref. set
#############  anno -- minimal for meffil normalization
anno_base <- openxlsx::read.xlsx(samplesheet)
anno <- anno_base[grep(parameter[1], anno_base$Primary_category), ] ## 1st filter for cancer type set at command line 
anno <- anno[!is.na(anno$pan_study),]
anno <- anno[,c("idat_filename", "Sex", "Sentrix_ID", "Sentrix_position", "Basename_Biowulf","Platform_methy","material_prediction")]
anno <- separate(anno, Sentrix_position, into = c("sentrix_row", "sentrix_col"), sep = 3, remove = TRUE)
anno$sentrix_row = gsub('R',"",anno$sentrix_row)
anno$sentrix_col = gsub('C',"",anno$sentrix_col)
names(anno) <- c("Sample_Name", "Sex", "Slide", "sentrix_row", "sentrix_col", "Basename", "Array", "Material")
anno$Slide <- as.numeric(as.character(anno$Slide))
anno <- anno %>% filter(!is.na(Basename))
anno <- anno %>% filter(!duplicated(Basename))

######### Read major reference set 
qc.objects <- readRDS(file.path(batchdir,"ref.qc.objects.rds"))
qc.objects <- qc.objects[names(qc.objects) %in% anno$Sample_Name]
gc()

if(is.na(newsamplesheet)){
	message("\nNo samplesheet for new samples. Skip normalization.\nSlowly (ETA up to ~30min) devour pancancer beta value matrix...");   ### 
	beta <- fread("pancancer/NORMALIZATION/Betas_pan_pc11_fixedeffectsMaterialArray_badcgpsremoved_badsamplesremoved.txt", nThread = cores)
	beta <- beta %>% column_to_rownames("V1")
    beta <- beta[, names(beta) %>% anno$idat_filename]
}else{
    newidats <- read.csv(newsamplesheet, header=F)
	message("\n", nrow(newidats), " New samples to normalize")

	anno_base_new <- anno_base[anno_base$idat_filename %in% newidats$V1,]
	anno_new <- anno_base_new[,c("idat_filename", "Sex", "Sentrix_ID", "Sentrix_position", "Basename_Biowulf","Platform_methy","material_prediction")]
	anno_new <- separate(anno_new, Sentrix_position, into = c("sentrix_row", "sentrix_col"), sep = 3, remove = TRUE)
	anno_new$sentrix_row = gsub('R',"",anno_new$sentrix_row)
	anno_new$sentrix_col = gsub('C',"",anno_new$sentrix_col)
	names(anno_new) <- c("Sample_Name", "Sex", "Slide", "sentrix_row", "sentrix_col", "Basename", "Array", "Material")
	anno_new$Slide <- as.numeric(as.character(anno_new$Slide))
	anno_new <- anno_new %>% filter(!is.na(Basename))
	anno_new <- anno_new %>% filter(!duplicated(Basename))
    
	files <- paste0(anno_new$Basename, "_Grn.idat.gz")
	files <- c(files, paste0(anno_new$Basename, "_Red.idat.gz"))
    if(sum(file.access(files))!=0){
		message("\nCan not read some files: ",  files[file.access(files)!=0]);
		stop("Crashed dummy error", call. = TRUE)
	}
    new.qc.objects <- meffil.qc(anno_new,
                             featureset = "common",
                             cell.type.reference = NA,
                             verbose = T)

    saveRDS(new.qc.objects, file = file.path(batchdirout,"new.qc.objects.rds"))
    new.qc.objects <- readRDS(file.path(batchdirout,"new.qc.objects.rds"))

# error detection: identify samples with errors and other samples not processed
# because they were on the same core with missed or corrupted file

	 new.qc.objects_error <- new.qc.objects[which(sapply(new.qc.objects, class) == 'try-error')]
	 new.qc.objects_error <- as.data.frame(do.call(c, new.qc.objects_error))
	 write.table(new.qc.objects_error, file.path(batchdirout, "newsamples.error.txt"), sep = "\t")
     if (nrow(new.qc.objects_error>0)) {
		 message("\nDamaged samples detected in the batch:\n", new.qc.objects_error)
     	 new.qc.objects <- new.qc.objects[-new.qc.objects_error]
	 }else{
         message("\nGood: no damaged samples in new batch.\n")
	 }
# isolate samples with errors in reference set
	qc.objects_error <- qc.objects[which(sapply(qc.objects, class) == 'try-error')]
	qc.objects_error <- as.data.frame(do.call(c, qc.objects_error))
# remove damaged samples from the reference qc.objects
	if (nrow(qc.objects_error>0)) {
		 message("\nDamaged samples detected in reference set:\n", qc.objects_error)
		 qc.objects <- qc.objects[-qc.objects_error]
	}else{
		 message("\nGood: No damaged samples detected in reference set.\n")
	}
	 
#identify bystander samples
    missing <- as.data.frame(setdiff(anno$Sample_Name, names(qc.objects)))
    names(missing) <- "missing"
    anno_missing <- anno_base[anno_base$idat_filename %in% missing$missing,]
    anno_missing <- anno_missing[,c("Sample", "Histology")]
    anno_missing ;
####combine qc objects and save
    qc.objects <- do.call(c, list(qc.objects, new.qc.objects))
    qc.objects <- qc.objects[!duplicated(names(qc.objects))]
    saveRDS(qc.objects,    file = file.path(batchdirout,"qc.objects.rds"))
    rm(new.qc.objects)
    qc.objects <- readRDS(file.path(batchdirout, "qc.objects.rds"))
	gc()
#########qc summary######## ~slow on large set
	qc.summary <- meffil.qc.summary(qc.objects)
	
	##extract bad CpGs for later
	bad.cpgs <- qc.summary$bad.cpgs$name
	saveRDS(bad.cpgs, file = file.path(batchdirout, "bad.cpgs.rds"))
	
	##identify and filter bad samples
	bad_samples_meth_unmeth <- as.data.frame(qc.summary$meth.unmeth.summary$tab)
	write.table(bad_samples_meth_unmeth, file = file.path(batchdirout, "qc_meth_unmeth.txt"), sep = "\t", row.names = FALSE)
	##bad_samples_meth_unmeth <- fread(file.path(batchdir, "qc_meth_unmeth.txt"))
	bad_samples_meth_unmeth <- bad_samples_meth_unmeth[bad_samples_meth_unmeth$outliers==TRUE,]
	bad_samples_control <- as.data.frame(qc.summary$controlmeans.summary$tab)
	write.table(bad_samples_control, file = file.path(batchdirout, "qc_control.txt"), sep = "\t", row.names = FALSE)
	##bad_samples_control <- fread(file.path(batchdir, "qc_control.txt"))
	bad_samples_control <- bad_samples_control[bad_samples_control$outliers==TRUE,]
	##bad_samples_control <- bad_samples_control[bad_samples_control$variable=="bisulfite1" | bad_samples_control$variable=="bisulfite2",]
	bad_samples_detectionp <- as.data.frame(qc.summary$sample.detectionp.summary$tab)
	write.table(bad_samples_detectionp, file = file.path(batchdirout, "qc_detectionp.txt"), sep = "\t", row.names = FALSE)
	##bad_samples_detectionp <- fread(file.path(batchdir, "qc_detectionp.txt"))
	bad_samples_detectionp <- bad_samples_detectionp[bad_samples_detectionp$outliers==TRUE,]
	bad_samples_lowbead <- as.data.frame(qc.summary$sample.beadnum.summary$tab)
	write.table(bad_samples_lowbead, file = file.path(batchdirout, "qc_lowbead.txt"), sep = "\t", row.names = FALSE)
	##bad_samples_lowbead <- fread(file.path(batchdir, "qc_lowbead.txt"))
	bad_samples_lowbead <- bad_samples_lowbead[bad_samples_lowbead$outliers==TRUE,]

    bad_samples <- c(
       bad_samples_meth_unmeth$sample.name,
       bad_samples_control$sample.name,
       bad_samples_detectionp$sample.name,
       bad_samples_lowbead$sample.name)
	message("\nExcluded samples with failed QCs:\n", bad_samples)

##  remove bad samples prior to performing quantile normalization
	qc.objects <- meffil.remove.samples(qc.objects, bad_samples)
    rm(qc.summary)
    gc()

## Check PC distribution for the optimum number of PC removal 
## The default settings is 11 in meffil.normalize.quantiles function
	pc.fit <- meffil.plot.pc.fit(qc.objects)
	pdf(file.path(batchdirout, "PC.rplot.pdf")) 
	print(pc.fit$plot)
	dev.off()

##normalize sample quantiles
    message("\n\nStart quantiles normalization:   ", Sys.time(),"...");
	norm.objects <- meffil.normalize.quantiles( qc.objects,
                      fixed.effects  = c("Material", "Array"), verbose = TRUE,
####                  random.effects = "Slide",      ## not working properly sometimes and adds >1hr
                      number.pcs = 11 );             ## between 4 and 10 depends how much variability you need to trim
	saveRDS(norm.objects, file = file.path(batchdirout,"norm.objects_pc11.rds"));
	norm.objects <- readRDS(file.path(batchdirout,"norm.objects_pc11.rds"))
	message("\n\nDone quantiles normalization:   ", Sys.time(), "\n now extracting beta...");
	options(mc.cores=70)
    cores <- options()$mc.cores

#### Generating normalized beta to the number of PCA (Part 0 from step above) for later reuse with new batches. 
    bad.cpgs <- readRDS(file.path(batchdirout, "bad.cpgs.rds"))
	norm.beta <- meffil.normalize.samples(norm.objects,
                                          just.beta = T,
                                          cpglist.remove = bad.cpgs,
                                          verbose = FALSE)
	saveRDS(norm.beta, file = file.path(batchdirout, "norm.signals_pc11.rds"))
    norm.signals   <- readRDS(file.path(batchdirout, "norm.signals_pc11.rds"))
    beta <- as.data.frame(norm.beta)
    fwrite(beta,  file.path(batchdirout,"Betas_pc11_fixedeffectsMaterialArray_badcgpsremoved_badsamplesremoved.txt"),
       sep = "\t", row.names = TRUE, nThread = 10 )
    beta <- fread(file.path(batchdirout,"Betas_pc11_fixedeffectsMaterialArray_badcgpsremoved_badsamplesremoved.txt"))
	beta <- beta %>% column_to_rownames("V1")
}


##probe filtering ~124K of probes (from both types of arrays)
HM450_hg19_Zhou <- fread(file.path("/data/MDATA/NormRcode/CNS13Krefset","HM450.hg19.manifest.tsv"),
                          data.table = FALSE,
                          stringsAsFactors = FALSE,
                          check.names = FALSE)

HM450_hg19_Zhou_MASK <- HM450_hg19_Zhou$probeID[HM450_hg19_Zhou$MASK_general == "TRUE"]
HM450_hg19_Zhou_chrmX <- HM450_hg19_Zhou$probeID[HM450_hg19_Zhou$CpG_chrm == "chrX"]
HM450_hg19_Zhou_chrmY <- HM450_hg19_Zhou$probeID[HM450_hg19_Zhou$CpG_chrm == "chrY"]
HM450_hg19_Zhou_rs <- HM450_hg19_Zhou$probeID[grepl("rs",HM450_hg19_Zhou$probeID)]
HM450_hg19_Zhou_ch <- HM450_hg19_Zhou$probeID[grepl("ch",HM450_hg19_Zhou$probeID)]

EPIC_hg19_Zhou <- fread(file.path("/data/MDATA/NormRcode/CNS13Krefset","EPIC.hg19.manifest.tsv"),
                         data.table = FALSE,
                         stringsAsFactors = FALSE,
                         check.names = FALSE)
EPIC_hg19_Zhou_MASK  <- EPIC_hg19_Zhou$probeID[EPIC_hg19_Zhou$MASK_general == "TRUE"]
EPIC_hg19_Zhou_chrmX <- EPIC_hg19_Zhou$probeID[EPIC_hg19_Zhou$CpG_chrm == "chrX"]
EPIC_hg19_Zhou_chrmY <- EPIC_hg19_Zhou$probeID[EPIC_hg19_Zhou$CpG_chrm == "chrY"]
EPIC_hg19_Zhou_rs <- EPIC_hg19_Zhou$probeID[grepl("rs",EPIC_hg19_Zhou$probeID)]
EPIC_hg19_Zhou_ch <- EPIC_hg19_Zhou$probeID[grepl("ch",EPIC_hg19_Zhou$probeID)]

EPIC_Illumina <- fread(file.path("/data/MDATA/NormRcode/CNS13Krefset","infinium-methylationepic-v-1-0-b5-manifest-file.csv"))
EPIC_Illumina_MFG <- EPIC_Illumina$IlmnID[EPIC_Illumina$MFG_Change_Flagged == "TRUE"]
EPIC_Illumina_chrmX <- EPIC_Illumina$IlmnID[EPIC_Illumina$CHR == "X"]
EPIC_Illumina_chrmY <- EPIC_Illumina$IlmnID[EPIC_Illumina$CHR == "Y"]
EPIC_Illumina_rs <- EPIC_Illumina$IlmnID[grepl("rs",EPIC_Illumina$IlmnID)]
EPIC_Illumina_ch <- EPIC_Illumina$IlmnID[grepl("ch",EPIC_Illumina$IlmnID)]

bad.cpgs_filter <- unique(c(HM450_hg19_Zhou_MASK,
                            HM450_hg19_Zhou_chrmX,
                            HM450_hg19_Zhou_chrmY,
                            HM450_hg19_Zhou_rs,
                            HM450_hg19_Zhou_ch,
                            EPIC_hg19_Zhou_MASK,
                            EPIC_hg19_Zhou_chrmX,
                            EPIC_hg19_Zhou_chrmY,
                            EPIC_hg19_Zhou_rs,
                            EPIC_hg19_Zhou_ch,
                            EPIC_Illumina_MFG,
                            EPIC_Illumina_chrmX,
                            EPIC_Illumina_chrmY,
                            EPIC_Illumina_rs,
                            EPIC_Illumina_ch))

bad.cpgs <- unique(c(bad.cpgs, bad.cpgs_filter))

beta <- beta[!rownames(beta) %in% bad.cpgs_filter,]
beta_t <- t(beta)
rm(beta)
gc()
beta_reduced <- beta_t[,Rfast::colVars(as.matrix(beta_t), std = TRUE, parallel = TRUE)>0.227]
beta_reduced <- as.data.frame(beta_reduced)
fwrite(beta_reduced, file.path(batchdirout,"Betas_pc11_fixedeffectsMaterialArray_badcgpsremoved_badsamplesremoved_0.227SDvariableprobes.txt"),
       sep = "\t", row.names = TRUE )

#big clenup
rm(list=setdiff(ls(), c("beta_reduced", "anno_base", "parameter", "batchdir", "newidats", "batchdirout", "cores")))
gc()

#####  ?? Optional reduction with principal components
# PC <- prcomp(beta_reduced, 
#             center = TRUE, 
#			 scale = FALSE,
#			 rank. = 200); 
#saveRDS(PC, file.path(batchdirout,"combo_200pc_prcomp5.rds"));
#PC <- readRDS(file.path(batchdirout,"combo_200pc_prcomp5.rds"));
#message("\nDone PCA:   ", Sys.time(), "\nGenerating Umap ...\n");

message("\nUmap straight on beta values   ", Sys.time(), "\ngenerating X and Y with uwot  ...");
umap <- uwot::umap(beta_reduced,    #PC$x -- if pca done 
                   #n_components = 2,
                   #pca = 25,
                   n_neighbors = 10,
                   metric = "cosine",
                   #y = anno_pan$Combined_class_match,
                   spread = 1,
                   min_dist = 0,
                   local_connectivity = 1,
                   bandwidth = 1)

message("\nDone unsupervised umap:   ", Sys.time(), "\n now densmap...");
dmap <- densmap(
        beta_reduced,
        n_components = 2L,
        dens_frac = 0.3,
        dens_lambda = 0.1,
        #var_shift = 0.1,
        n_neighbors = 10L,
        metric = "cosine",
        n_epochs = 750L,
        learning_rate = 1,
        init = "spectral",
        Y_init = NULL,
        min_dist = 0.0001,
        spread = 1,
        set_op_mix_ratio = 1,
        local_connectivity = 1L,
        repulsion_strength = 1,
        negative_sample_rate = 5L,
        transform_queue_size = 4,
        random_state = NULL,
        angular_rp_forest = FALSE,
        target_n_neighbors = -1,
        target_weight = 0.5
)
umap <-cbind(umap, dmap)
umap <- as.data.frame(umap)
rownames(umap) <- rownames(beta_reduced)
umap <- umap %>% rownames_to_column("idat_filename")
names(umap) <- c("idat_filename", "x1","y1","x2","y2");

newanno <- anno_base[anno_base$idat_filename %in% umap$idat_filename,];
newanno <- newanno[!duplicated(newanno$idat_filename),];

## Pick you own but do not forget to fix Rshiny app file later
metadatatouse <- c(
   "idat_filename",
   "Sample",
   "Sex",
   "Age",
   "material_prediction",
   "ABSOLUTE",
   "ESTIMATE",
   "LUMP",
   "Location_general",
   "Location_specific",
   "OS_months",
   "OS_status",
   "PFS_months",
   "PFS_status",
   "Histology",
   "Variants",
   "Fusions/translocations",
   "Assay",
   "Primary_study",
   "MCF_v11b6",
   "MCF_v11b6_score",
   "MC_v11b6",
   "MC_v11b6_score")

#Some text conditioning for the r-shiny code compartibility
outf <- merge(umap, newanno[,metadatatouse],   by = "idat_filename");
outf$Molecular <- paste(outf$Variants,  
                 '|',  outf$"Fusions/translocations",
                 '|',  outf$Assay)
outf[grep("NA | NA | NA", outf$Molecular), "Molecular"] <- ""; 
outf <- outf[, -c(20,21,22)]
outf$Study <- "RefPool";
outf[outf$idat_filename %in% newidats$V1, "Study"] <- "compass"
names(outf) <- c(
   "idat_filename",
   "x1","y1","x2","y2",
   "Sample", "Sex", "Age",
   "material_prediction",
   "RFpurity.ABSOLUTE", "RFpurity.ESTIMATE",  "LUMP",
   "Location_general",   "Location_specific", 
   "OS_months",  "OS_status", "PFS_months",  "PFS_status",
   "Histology",  "Primary_study",
   "CNS.MCF",      "CNS.MCF.score",
   "CNS.Subclass",  "CNS.Subclass.score",
   "Molecular",  "Study")
outf <- arrange(outf, Study);

fwrite(outf, file.path(batchdirout,"/umap_2shiny.txt"), 
    row.names=TRUE, sep = "\t", nThread = cores);

message(Sys.time());
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
message(" End of normalization and UMAP script ");
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");