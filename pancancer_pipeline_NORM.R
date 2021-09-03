## Drew's normalization pipeline with small tweaks by Rust
## August 2021 :: pancancer set

args = commandArgs(trailingOnly = TRUE)
# args = "Placenta:NewCNSdir:TRANSFER/newbatch.txt"
# args = "Neuropathology:CNS_Aug13:TRANSFER/newbatch_aug13.txt" 
# args = "Sarcoma:Sarcoma_Aug13:TRANSFER/newbatch_aug13.txt" 
if(length(args) != 1){
      stop("!!! Crash !!!\nParameter settings missed in sbatch: [cancer:outputfolder:newsamplesheet] \nWhere:\n   cancer \t  -- [Primary_category] from annotation table column [CC]\n   outputfolder   -- directory to save at /data/MDATA/TRANSFER\n   newsamplesheet -- full path to the new batch metadata\nTry to rerun sbatch with appropriate settings for like this:\n\nsbatch /data/MDATA/NormRcode/norm.sh  Thoracic|Bone and soft tissue:NewCNSdir:TRANSFER/newbatch.txt\n\n")
}else{
     parameter = as.character(args[1])
     parameter = unlist(strsplit(parameter, ":"))
     message("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", Sys.time(),"  Started with parameters:\nCancer:     \t", parameter[1], "\nOutputdir:\t", parameter[2], "\nSamples:\t", parameter[3], "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
}
options(scipen = 999)
library(data.table)
library(parallel)
library(readxl) #read_excel
suppressMessages(library("openxlsx")) #, lib.loc="/home/prattdw/R/4.0/library/")
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(meffil))
suppressMessages(library(Rfast))
suppressMessages(library(tibble))
library(RSpectra)
library(densvis)
setwd("/data/MDATA")
options(mc.cores=140)
cores <- options()$mc.cores
gc()

### Set paths around
idatrootfolder = "pancancer/iScan_raw";
batchdir       = "pancancer/NORMALIZATION"; 
batchdirout    = file.path("TRANSFER", parameter[2]);    # working folder
newsamplesheet = as.character(parameter[3])
dir.create(batchdirout, recursive = TRUE)
samplesheet    = "bams_RNAseq/Sample_sheet_master.xlsm"; # the database


#############  Inject static samplesheet
### Part 0 ##  anno_base -- rich file for the ref. set
#############  anno -- minimal for meffil normalization
message("\nLoading excel masterfile...")
anno_base <- openxlsx::read.xlsx(samplesheet)
if (parameter[1] == "Sarcoma") {
   anno <- anno_base[grep("Case", anno_base$panSARCOMA), ]             ## Sarcoma flag is in different column    
}else{
   anno <- anno_base[grep(parameter[1], anno_base$Primary_category), ] 
}
anno <- anno[!is.na(anno$pan_study),]
anno <- anno[,c("idat_filename", "Sex", "Sentrix_ID", "Sentrix_position", "Basename_Biowulf","Platform_methy","material_prediction")]
anno <- separate(anno, Sentrix_position, into = c("sentrix_row", "sentrix_col"), sep = 3, remove = TRUE)
anno$sentrix_row = gsub('R',"",anno$sentrix_row)
anno$sentrix_col = gsub('C',"",anno$sentrix_col)
names(anno) <- c("Sample_Name", "Sex", "Slide", "sentrix_row", "sentrix_col", "Basename", "Array", "Material")
anno$Slide <- as.numeric(as.character(anno$Slide))
anno <- anno %>% filter(!is.na(Basename))
anno <- anno %>% filter(!duplicated(Basename))

######### Read major reference set (from prebuilt file)
message("\nLoading ref.qc.objects.rds...")
qc.objects <- readRDS(file.path(batchdir,"ref.qc.objects.rds"))
qc.objects <- qc.objects[names(qc.objects) %in% anno$Sample_Name]
gc()


if(is.na(newsamplesheet)){
	message("\nNo samplesheet for new samples. Skip normalization.\nSlowly (ETA up to ~30min) devour pancancer beta value matrix...");   ### 
	beta <- fread("pancancer/NORMALIZATION/Betas_pan_pc11_fixedeffectsMaterialArray_badcgpsremoved_badsamplesremoved.txt", nThread = 32)
	beta <- beta %>% column_to_rownames("V1")
    beta <- beta[, names(beta) %>% anno$idat_filename]
}else{
    anno_new <- read.csv(newsamplesheet, header=F)  ## getting barcodes for new samples
    allnewsamples <- as.character(anno_new$V1)
	message("\n", nrow(anno_new), " New samples to normalize")
    newdirs = unique(matrix(unlist(strsplit(anno_new$V1, "_")), ncol=2, byrow=TRUE)[,1])

	nx = 0; ## All chips counter
	kk = 0; ## With samplesheet counter 

	newsamplesheet = c();
	for(centrix in newdirs){
	   nx     = nx+1 
	   Sample = c(); 
	   idat_filename =c();
	   idat   = c();
	   Sentrix_ID = c();
	   Sentrix_position = c();
	   material_prediction  = c();
	   Platform_methy = c();
	   Age = c();
	   Sex = c();
	   matched_cases = c();
	   Location_general = c();
	   Location_specific = c();
	   Neoplastic = c();
	   Primary_category = c();
	   CNS_study = c();
	   OS_months = c();
	   OS_status = c();
	   PFS_months = c();
	   PFS_status = c();
	   Histology = c();
	   Molecular = c();
	   Study = c();
	   predFFPE = c();
	   CNS.MCF = c();
	   CNS.MCF.score = c();
	   CNS.Subclass = c();
	   CNS.Subclass.score = c();
	   RFpurity.ABSOLUTE = c();
	   RFpurity.ESTIMATE = c();
	   LUMP = c();
	   knnsheet <- paste0("/data/MDATA/compass/iScan_raw/", centrix, "/", centrix,"_KNN.combined.csv");
	   if((file.exists(knnsheet))&(file.info(knnsheet)$size)>1){
		   kk = kk +1;
			 rawsheetfull  <- read.csv(knnsheet, row.names=1);
			 sampletype = rawsheetfull[,grep('Sample_Plate', names(rawsheetfull), ignore.case = T, perl = T)];
#			 rawsheet <- rawsheetfull[grepl("Clinical|Brain|Sarcoma|CBTN", sampletype, ignore.case = TRUE, perl = TRUE), ]; ## Only some rows/types go in
			 rawsheet <- rawsheetfull; ## All in 
			 nnn                  = rep(kk, nrow(rawsheet));
			 Sample               = rawsheet[,grep('SAMPLE_NAME', names(rawsheet), ignore.case = T, perl = T)];
			 idat_filename        = rawsheet[,which( names(rawsheet) == 'ID')];
			 idat                 = idat_filename;
			 Sentrix_ID           = rawsheet[,grep('SENTRIX_ID', names(rawsheet), ignore.case = T, perl = T)];
			 Sentrix_position     = rawsheet[,grep('SENTRIX_POSITION', names(rawsheet), ignore.case = T, perl = T)];
			 material_prediction  = rawsheet[,grep('PREDFFPE', names(rawsheet), ignore.case = T, perl = T)];
			 Platform_methy       = rawsheet[,grep('SAMPLE_GROUP', names(rawsheet), ignore.case = T, perl = T)];
			   Platform_methy     = replace(Platform_methy, grep('EPIC', Platform_methy, ignore.case = T), "HumanMethylationEPIC");
			   Platform_methy     = replace(Platform_methy, grep('450', Platform_methy, ignore.case = T, perl = T), "HumanMethylation450");
			 Age                  = rawsheet[,grep('AGE', names(rawsheet), ignore.case = T, perl = T)];
			   Age                = replace(Age, grep("\\D", Age), NA);
			 Sex                  = rawsheet[,grep('GENDER|SEX', names(rawsheet), ignore.case = T, perl = T)];
			   Sex                = replace(Sex, grep("FEMALE",  Sex, ignore.case = T, perl = T), "F" );
			   Sex                = replace(Sex, grep("MALE",    Sex, ignore.case = T, perl = T), "M" );
			   Sex                = replace(Sex, grep("unknown", Sex, ignore.case = T, perl = T),  NA );
			   Sex                = replace(Sex, grep("TBD",     Sex, ignore.case = T, perl = T),  NA );
			 matched_cases        = rawsheet[,grep('NEARESTNEIGHBOR', names(rawsheet), ignore.case = T, perl = T)];
			 Location_general     = rawsheet[,grep('TUMOR_SITE',      names(rawsheet), ignore.case = T, perl = T)];
			 Location_specific    = rawsheet[,grep('SURGICAL_CASE',   names(rawsheet), ignore.case = T, perl = T)];
			 if(parameter[1] == "Sarcoma") {
                 Neoplastic           = rep("Neuropathology", nrow(rawsheet));
			 }else{
                 Neoplastic           = rep("Sarcoma", nrow(rawsheet));
			 }
			 Primary_category     = rep("Case",           nrow(rawsheet));
			 CNS_study            = rep("Case",           nrow(rawsheet));
			 OS_months            = rep(NA,               nrow(rawsheet));
			 OS_status            = rep(NA,               nrow(rawsheet));
			 PFS_months           = rep(NA,               nrow(rawsheet));
			 PFS_status           = rep(NA,               nrow(rawsheet));
			 Histology            = rawsheet[,grep('DIAGNOSIS', names(rawsheet), ignore.case = T, perl = T)];
			 Molecular            = rawsheet[,grep('NOTES',     names(rawsheet), ignore.case = T, perl = T)];
			 Study                = rep("compass",        nrow(rawsheet));
			 NCI_METRIC           = rep("",        nrow(rawsheet));
			 predFFPE             = rawsheet[,grep('MATERIAL_TYPE',     names(rawsheet), ignore.case = T, perl = T)];
			 CNS.MCF              = rawsheet[,which( names(rawsheet) == 'MCF1')]; 
			 CNS.MCF.score        = rawsheet[,grep('MCF1.SCORE',        names(rawsheet), ignore.case = T)];
			 CNS.Subclass         = rawsheet[,which( names(rawsheet) == 'Class1')]; 
			 CNS.Subclass.score   = rawsheet[,grep('CLASS1.SCORE',      names(rawsheet), ignore.case = T)];
			 RFpurity.ABSOLUTE    = rawsheet[,grep('RFPURITY.ABSOLUTE', names(rawsheet), ignore.case = T, perl = T)];
			 RFpurity.ESTIMATE    = rawsheet[,grep('RFPURITY.ESTIMATE', names(rawsheet), ignore.case = T, perl = T)];
			 LUMP                 = rawsheet[,grep('LUMP',              names(rawsheet), ignore.case = T, perl = T)];
			 Basename             = rawsheet$Basename;
			clms <- list(nnn, Sample,idat_filename,idat,Sentrix_ID,Sentrix_position,material_prediction,Platform_methy,Age,Sex,matched_cases, Location_general, Location_specific, Neoplastic, Primary_category, CNS_study, OS_months, OS_status,PFS_months, PFS_status, Histology,Molecular, Study,NCI_METRIC,predFFPE,CNS.MCF, CNS.MCF.score,CNS.Subclass, CNS.Subclass.score, RFpurity.ABSOLUTE, RFpurity.ESTIMATE, LUMP,Basename)
			names(clms) <- c("nn", "Sample", "idat_filename", "idat", "Sentrix_ID", "Sentrix_position", "material_prediction", "Platform_methy", "Age", "Sex", "matched_cases", "Location_general", "Location_specific", "Neoplastic", "Primary_category", "CNS_study", "OS_months","OS_status", "PFS_months", "PFS_status", "Histology", "Molecular", "Study", "NCI_METRIC", "predFFPE","CNS.MCF", "CNS.MCF.score", "CNS.Subclass", "CNS.Subclass.score", "RFpurity.ABSOLUTE", "RFpurity.ESTIMATE", "LUMP", "Basename");
			colszs  <- sapply(clms,"length");  
			for(i in 1:length(colszs)){
				if(colszs[i] == 0){
				  clms[[names(colszs[i])]] <- c(rep(NA, length.out = max(colszs)));
				}
			 };
			 newsamplesheet = rbind( newsamplesheet , as.data.frame(clms) );
			 message(centrix, " : has ", nrow(rawsheet), " samples in samplesheet");
	   }else{
			 message(knnsheet, " : is not here");
	   }
	};
    
	## Intersect samples only in original list and  no missed metadata in _KNN.combined.csv 
	newsamplesheet <- newsamplesheet[newsamplesheet$idat_filename %in% anno_new$V1,]
    anno_new <-  anno_new[anno_new$V1 %in% newsamplesheet$idat_filename,]
    ms <- length(setdiff(allnewsamples, newsamplesheet$idat_filename ))
	message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\t", 
	 nx, "\t-- chips processed\n\t", 
	 kk, "\t-- with KNN.combined.csv\n\t", 
     ms, "\t-- sample without annotations\n\t",
     nrow(newsamplesheet),"\t-- new samples collected\n",
	"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
   	write.csv(newsamplesheet, paste0(batchdirout,"/newsamples.csv", row.names=F))
    # Some chips may not have all samples included *_KNN.combined.csv those samples will be filtered out from newsamplesheet
    # run roboreporter.pl script prior, to make *_KNN.combined.csv spreadsheet 
    # sample should be in Clinical category in samplesheet.csv)

##### Blending new samples with Drew database {anno_base}
	anno_base_new                  <- anno_base[FALSE,]
    anno_base_new[1:nrow(newsamplesheet),] <- NA
    anno_base_new$order  <- (tail(anno_base$order,1)+1):(nrow(newsamplesheet)+tail(anno_base$order,1))
	anno_base_new$Sample           <- newsamplesheet$Sample
    anno_base_new$idat_filename    <- newsamplesheet$idat_filename
    anno_base_new$Basename_Biowulf <- newsamplesheet$Basename
    anno_base_new$ESTIMATE         <- newsamplesheet$RFpurity.ESTIMATE
    anno_base_new$ABSOLUTE         <- newsamplesheet$RFpurity.ABSOLUTE
    anno_base_new$LUMP             <- newsamplesheet$LUMP
    anno_base_new$GSM_accession    <- "compass"
    anno_base_new$Primary_study    <- "compass"             #used in umap to label new samples
    anno_base_new$Primary_category <- parameter[1]          #Neuropathology:Endocrine:Placenta etc 
    anno_base_new$pan_study        <- newsamplesheet$Study  #compass too asigned above
    anno_base_new$idat             <- newsamplesheet$idat
    anno_base_new$Sentrix_ID       <- newsamplesheet$Sentrix_ID
    anno_base_new$Sentrix_position <- newsamplesheet$Sentrix_position
    anno_base_new$material_prediction <- newsamplesheet$material_prediction
    anno_base_new$Platform_methy   <- newsamplesheet$Platform_methy
    anno_base_new$Center_methy     <- "NIH"
    anno_base_new$Filename_methy   <- newsamplesheet$idat_filename
    anno_base_new$Accession_methy  <- "NIH"
    anno_base_new$Age              <- newsamplesheet$Age
    anno_base_new$Sex              <- newsamplesheet$Sex
	anno_base_new$matched_cases    <- newsamplesheet$matched_cases
	anno_base_new$Location_general <- newsamplesheet$Location_general
    anno_base_new$Location_specific<- newsamplesheet$Location_specific
    anno_base_new$Neoplastic       <- newsamplesheet$Neoplastic
    anno_base_new$Primary_category <- newsamplesheet$Primary_category
    anno_base_new$NCI_METRIC       <- newsamplesheet$NCI_METRIC
    anno_base_new$panCNS           <- newsamplesheet$CNS_study
    anno_base_new$OS_months        <- newsamplesheet$OS_months
    anno_base_new$OS_status        <- newsamplesheet$OS_status
    anno_base_new$PFS_months       <- newsamplesheet$PFS_months
    anno_base_new$PFS_status       <- newsamplesheet$PFS_status
    anno_base_new$Histology        <- newsamplesheet$Histology
    anno_base_new$Variants         <- newsamplesheet$Molecular
    anno_base_new$Primary_study    <- newsamplesheet$Study
    anno_base_new$material_prediction <- newsamplesheet$predFFPE
    anno_base_new$MCF1_v11b6       <- newsamplesheet$CNS.MCF
    anno_base_new$MCF1_v11b6_score <- newsamplesheet$CNS.MCF.score
    anno_base_new$v11b6            <- newsamplesheet$CNS.Subclass
    anno_base_new$MCF2_v11b6_score <- newsamplesheet$CNS.Subclass.score  ##This score is not in drew metadata

	anno_base <- rbind(anno_base, anno_base_new)  ## for later  joining with UMAP 
    
### prepare minimal annotation for the meffil
	anno_new <- anno_base_new[,c("idat_filename", "Sex", "Sentrix_ID", "Sentrix_position", "Basename_Biowulf","Platform_methy","material_prediction")]
	anno_new <- separate(anno_new, Sentrix_position, into = c("sentrix_row", "sentrix_col"), sep = 3, remove = TRUE)
	anno_new$sentrix_row = gsub('R',"",anno_new$sentrix_row)
	anno_new$sentrix_col = gsub('C',"",anno_new$sentrix_col)
	names(anno_new) <- c("Sample_Name", "Sex", "Slide", "sentrix_row", "sentrix_col", "Basename", "Array", "Material")
	anno_new$Slide <- as.numeric(as.character(anno_new$Slide))
	anno_new <- anno_new %>% filter(!is.na(Basename))
	anno_new <- anno_new %>% filter(!duplicated(Basename))

	files <- paste0(anno_new$Basename, "_Grn.idat")
	files <- c(files, paste0(anno_new$Basename, "_Red.idat"))
    if(sum(file.access(files))!=0){
		message("\nCan not read some files: ",  files[file.access(files)!=0]);
		stop("Crashed dummy error", call. = TRUE)
	}else{message("\nFound all input idat files for new batch.\n")}

    ## Generate QC objects for new batch
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
#### !! Merging new batch and reference set here

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

    bad_samples_all <- c(
       bad_samples_meth_unmeth$sample.name,
       bad_samples_control$sample.name,
       bad_samples_detectionp$sample.name,
       bad_samples_lowbead$sample.name)
    newsmplsfailedqc <- intersect(anno_base_new$idat_filename, bad_samples_all)
    bad_samples <- setdiff(bad_samples_all,  newsmplsfailedqc)
	qc.objects <- meffil.remove.samples(qc.objects, bad_samples)
    rm(qc.summary)
    gc()

## Check PC distribution for the optimum number of PC removal 
## The default settings is 11 in meffil.normalize.quantiles function
    message("\nGenerating PDF with principal components plot:\n")
	pc.fit <- meffil.plot.pc.fit(qc.objects)
	pdf(file.path(batchdirout, "PC.rplot.pdf")) 
	print(pc.fit$plot)
	dev.off()

##normalize sample quantiles
    message("\n\nStart quantiles normalization:   ", Sys.time(),"...");
	norm.objects <- meffil.normalize.quantiles( qc.objects,
                      fixed.effects  = c("Material", "Array"), verbose = TRUE,
 ##                   random.effects = "Slide",      ## not working properly sometimes and adds >1hr
                      number.pcs = 11 );             ## between 4 and 10 depends how much variability you need to trim
	saveRDS(norm.objects, file = file.path(batchdirout,"norm.objects_pc11.rds"));
	norm.objects <- readRDS(file.path(batchdirout,"norm.objects_pc11.rds"))
	message("\n\nDone quantiles normalization:   ", Sys.time(), "\n now extracting beta...");
	options(mc.cores=90)
    cores <- options()$mc.cores

#### Generating normalized beta to the number of PCA (Part 0 from step above) for later reuse with new batches. 
    bad.cpgs <- readRDS(file.path(batchdirout, "bad.cpgs.rds"))
	norm.beta <- meffil.normalize.samples(norm.objects,
                                          just.beta = T,
                                          cpglist.remove = bad.cpgs,
                                          verbose = FALSE)
	saveRDS(norm.beta, file = file.path(batchdirout, "norm.signals_pc11.rds"))
    norm.beta   <- readRDS(file.path(batchdirout, "norm.signals_pc11.rds"))
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
message("\nRemoving ", length(bad.cpgs), " of bad cpgs: mt, X, Y, SNPs");
beta <- beta[!rownames(beta) %in% bad.cpgs_filter,]
beta_t <- t(beta)
rm(beta)
beta_reduced <- beta_t[,Rfast::colVars(as.matrix(beta_t), std = TRUE, parallel = TRUE)>0.227]
beta_reduced <- as.data.frame(beta_reduced)
fwrite(beta_reduced, file.path(batchdirout,"Betas_pc11_fixedeffectsMaterialArray_badcgpsremoved_badsamplesremoved_0.227SDvariableprobes.txt"),
       sep = "\t", row.names = TRUE )
message("\nSaved transposed reduced file with beta for ", dim(beta_reduced)[2], " top variable cpgs with SD>0.227\n");
#beta_reduced <- read.csv(file.path(batchdirout,"Betas_pc11_fixedeffectsMaterialArray_badcgpsremoved_badsamplesremoved_0.227SDvariableprobes.txt"),sep = "\t", row.names =1 )

message("\nExcluded samples with failed QCs:\n")
print(bad_samples);
if (length(newsmplsfailedqc)>0) {
	message("\n^^^^^^^^^^^^^^^^^\n!!! ATTENTION !!!  You have QC outliers in your batch.\n^^^^^^^^^^^^^^^^^  Check QCs for those samples:")
	BADAS <- anno_base_new[anno_base_new$idat_filename %in%  newsmplsfailedqc, c("Sample","idat_filename", "Primary_category", "Primary_study")]
	print(BADAS)
}else{message("\nGood: No QC outliers in new batch\n")}; 

#big cleanup
rm(list=setdiff(ls(), c("beta_reduced", "anno_base", "newsmplsfailedqc", "newsamplesheet", "parameter", "batchdir", "batchdirout", "cores")))
gc()

##### fast PCA via RSpectra SVD by Martin Sill m.sill@dkfz.de (2018)
##### if(!require(RSpectra)) install.packages("RSpectra")

prcomp_svds <-
  function(x, retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL, k=nrow(x), ...)
  {
    chkDots(...)
    x <- as.matrix(x)
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0))
      stop("cannot rescale a constant/zero column to unit variance")
    s <- svds(x, k)
    s$d <- s$d / sqrt(max(1, nrow(x) - 1))
    if (!is.null(tol)) {
      ## we get rank at least one even for a 0 matrix.
      rank <- sum(s$d > (s$d[1L]*tol))
      if (rank < ncol(x)) {
        s$v <- s$v[, 1L:rank, drop = FALSE]
        s$d <- s$d[1L:rank]
      }
    }
    dimnames(s$v) <-
      list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
    r <- list(sdev = s$d, rotation = s$v,
              center = if(is.null(cen)) FALSE else cen,
              scale = if(is.null(sc)) FALSE else sc)
    if (retx) r$x <- x %*% s$v
    class(r) <- "prcomp"
    r
  }
#####  This is optional reduction variable reduction. Makes picture sharper for big sets
message("\nRun principal component reduction on reduced beta", Sys.time(), "...");
 PC <- prcomp_svds(beta_reduced, 
             center = TRUE, 
			 scale = FALSE,
			 rank. = 200); 
saveRDS(PC, file.path(batchdirout,"combo_200pc_prcomp.rds"));
PC <- readRDS(file.path(batchdirout,"combo_200pc_prcomp.rds"));
beta_reduced <- PC$x;

message("\nUmap straight on beta values   ", Sys.time(), "\ngenerating X and Y with uwot  ...");
umap <- uwot::umap(beta_reduced,  
                   #n_components = 2,
                   #pca = 25,
                   n_neighbors = 10,     # 10
                   metric = "cosine",   
                   #y = anno_pan$Combined_class_match,
                   spread = 1,          # 1
                   min_dist = 0.1,      # 0
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
if(parameter[1] == "Sarcoma"){
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
	   "MCF1_v11b6",
	   "MCF1_v11b6_score",
	   "NCI_METRIC",        ## ? those are actually close with "v11b6" ? 
	   "NCI_METRIC")
}else{
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
	   "MCF1_v11b6", 
	   "MCF1_v11b6_score",
	   "v11b6",
	   "MCF2_v11b6_score")
}

#Some text conditioning for the r-shiny code compartibility
outf <- merge(umap, newanno[,metadatatouse],   by = "idat_filename");
outf$Molecular <- paste(outf$Variants,  
                 '|',  outf$"Fusions/translocations",
                 '|',  outf$Assay)
outf[grep("NA . NA . NA", outf$Molecular), "Molecular"] <- ""; 
# Add attention label for molecular column for QC failed samples
outf[outf$idat_filename %in% newsmplsfailedqc, "Molecular"] <- paste(outf[outf$idat_filename %in% newsmplsfailedqc, "Molecular"], "Bad QC")
outf <- outf[, -c(20,21,22)]
## Override compass labels could be some sample info loss due to duplications with the reference data
outf$Study <- "RefPool";
outf[outf$idat_filename %in% newsamplesheet$idat_filename, "Study"] <- "compass"
names(outf) <- c(
   "idat_filename",
   "x1","y1","x2","y2",  
   "Sample", 
   "Sex", 
   "Age",
   "material_prediction",
   "RFpurity.ABSOLUTE", 
   "RFpurity.ESTIMATE",  
   "LUMP",
   "Location_general",  
   "Location_specific", 
   "OS_months",
   "OS_status",
   "PFS_months",
   "PFS_status",
   "Histology",       
   "Primary_study",
   "CNS.MCF",
   "CNS.MCF.score",
   "CNS.Subclass",
   "CNS.Subclass.score",
   "Molecular",
   "Study")
outf <- arrange(outf, Study);
fwrite(outf, file.path(batchdirout,"umap_2shiny.txt"), 
    row.names=TRUE, sep = "\t", nThread = 32);

message(Sys.time());
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
message(" End of normalization and UMAP script ");
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");