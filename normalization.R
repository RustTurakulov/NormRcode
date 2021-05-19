

## 5/18/2021 :: With two UMAP plots on output 

options(scipen = 999)
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
library(data.table)

options(mc.cores=70)
cores <- options()$mc.cores

## compass dir
idatrootfolder = "compass/iScan_raw/";
batchdir      <- "NormRcode/CNS13Krefset/"; 
batchdirout   <- "TRANSFER/SAMPLESHEETS/compass"; 
samplesheet   <- "NormRcode/CNS13Krefset/Compass_13K.xlsx";

#############  Inject static samplesheet
### Part 0 ##  anno_base -- rich file for the ref. set
#############  anno -- minimal for meffil normalization

anno_base <- openxlsx::read.xlsx(samplesheet)
anno <- anno_base %>% filter( CNS_study == "Case" );
anno <- anno[,c( "idat_filename", "Sex", "Sentrix_ID", "Sentrix_position", "Platform_methy","material_prediction")];
anno$Basename <- paste0("idat/", anno$idat_filename);
anno$material_prediction[anno$material_prediction!="FFPE"]<-"Frozen";
anno$Sex = gsub("FEMALE", "F", toupper(anno$Sex));
anno$Sex = gsub("MALE",   "M", toupper(anno$Sex));
anno$Sex[!(anno$Sex %in% c("M","F"))]<- NA;
anno <- separate(anno, Sentrix_position, into = c("sentrix_row", "sentrix_col"), sep = 3, remove = TRUE)
anno$sentrix_row = gsub("\\R","",anno$sentrix_row)
anno$sentrix_col = gsub("\\C","",anno$sentrix_col)
names(anno) <- c("Sample_Name", "Sex", "Slide", "sentrix_row", "sentrix_col", "Array", "Material", "Basename");
anno$Slide <- as.numeric(as.character(anno$Slide))
anno <- anno %>% filter(!is.na(Basename))
anno <- anno %>% filter(!duplicated(Basename))


#############  Collect new samples from compass folder 
### Part 1 ##  and check those against reference dataset samplesheet    
###        ##  note samples should be preprocessed with Zhichao's pipeline to generate DKFZ scores
#############  use 'roboreporter.pl' this one runs on compass folder in MDATA 

compdirs  <- list.files(idatrootfolder, pattern = "^[0-9]+");
newdirs   <- setdiff(compdirs, anno$Slide)

nn = 0; 
kk = 0;

newsamplesheet = c();
for(centrix in newdirs){
   nn     = nn+1 
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
   knnsheet <- paste0(idatrootfolder, centrix, "/", centrix,"_KNN.combined.csv");
   if(file.exists(knnsheet)){
	   kk = kk +1;
       rawsheet  <- read.csv(knnsheet, row.names=1);
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
		 matched_cases        = rawsheet[,grep('NEARESTNEIGHBOR', names(rawsheet), ignore.case = T, perl = T)];
		 Location_general     = rawsheet[,grep('TUMOR_SITE',      names(rawsheet), ignore.case = T, perl = T)];
		 Location_specific    = rawsheet[,grep('SURGICAL_CASE',   names(rawsheet), ignore.case = T, perl = T)];
		 Neoplastic           = rep("Neuropathology", nrow(rawsheet));
		 Primary_category     = rep("Case",           nrow(rawsheet));
		 CNS_study            = rep("Case",           nrow(rawsheet));
		 OS_months            = rep(NA,               nrow(rawsheet));
		 OS_status            = rep(NA,               nrow(rawsheet));
		 PFS_months           = rep(NA,               nrow(rawsheet));
		 PFS_status           = rep(NA,               nrow(rawsheet));
		 Histology            = rawsheet[,grep('DIAGNOSIS', names(rawsheet), ignore.case = T, perl = T)];
		 Molecular            = rawsheet[,grep('NOTES',     names(rawsheet), ignore.case = T, perl = T)];
		 Study                = rep("compass",        nrow(rawsheet));
		 predFFPE             = rawsheet[,grep('MATERIAL_TYPE',     names(rawsheet), ignore.case = T, perl = T)];
		 CNS.MCF              = rawsheet[,which( names(rawsheet) == 'MCF1')]; 
		 CNS.MCF.score        = rawsheet[,grep('MCF1.SCORE',        names(rawsheet), ignore.case = T)];
		 CNS.Subclass         = rawsheet[,which( names(rawsheet) == 'Class1')]; 
		 CNS.Subclass.score   = rawsheet[,grep('CLASS1.SCORE',      names(rawsheet), ignore.case = T)];
		 RFpurity.ABSOLUTE    = rawsheet[,grep('RFPURITY.ABSOLUTE', names(rawsheet), ignore.case = T, perl = T)];
		 RFpurity.ESTIMATE    = rawsheet[,grep('RFPURITY.ESTIMATE', names(rawsheet), ignore.case = T, perl = T)];
		 LUMP                 = rawsheet[,grep('LUMP',              names(rawsheet), ignore.case = T, perl = T)];           
		 Basename             = rawsheet$Basename;
		 reshaped = data.frame(nn = nn, Sample = Sample, idat_filename = idat_filename, idat = idat, Sentrix_ID = Sentrix_ID, 
                        Sentrix_position = Sentrix_position, material_prediction = material_prediction, Platform_methy = Platform_methy,
                        Age = Age, Sex = Sex, matched_cases = matched_cases,  Location_general = Location_general, Location_specific = Location_specific,
			            Neoplastic = Neoplastic, Primary_category = Primary_category, CNS_study = CNS_study, OS_months = OS_months, OS_status = OS_status,
                        PFS_months = PFS_months, PFS_status, Histology = Histology, Molecular = Molecular, Study = Study, predFFPE = predFFPE, CNS.MCF = CNS.MCF,
                        CNS.MCF.score = CNS.MCF.score, CNS.Subclass = CNS.Subclass, CNS.Subclass.score = CNS.Subclass.score, RFpurity.ABSOLUTE = RFpurity.ABSOLUTE,
                        RFpurity.ESTIMATE = RFpurity.ESTIMATE, LUMP = LUMP, Basename = Basename );
        
		 newsamplesheet = rbind( newsamplesheet , reshaped );

		 message(centrix, " : has ", nrow(reshaped), " samples in samplesheet");
   }else{
	   message(knnsheet, " : is not here");
   }
};

## Putting both sampleshhets together without Basename column 
combo_samplesheet <- rbind(newsamplesheet[,-32], anno_base);
fwrite(combo_samplesheet, file.path(batchdirout, "combo_samplesheet.tsv"),
       sep = "\t", row.names = TRUE,  nThread = cores);

message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n  ", 
        nn, "\t-- chips processed\n  ", 
		kk, "\t-- with CNS sample sheets\n ", 
        nrow(newsamplesheet),"\t-- new samples collected\n",
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

## Reshaping of reshaped new annotation for meffil
new_anno <- newsamplesheet %>% filter( CNS_study == "Case" );
new_anno <- new_anno[,c( "idat_filename", "Sex", "Sentrix_ID", "Sentrix_position", "Platform_methy","material_prediction", "Basename")];
new_anno$material_prediction[new_anno$material_prediction!="FFPE"]<-"Frozen";
new_anno$Sex[!(new_anno$Sex %in% c("M","F"))]<- NA;
new_anno <- separate(new_anno, Sentrix_position, into = c("sentrix_row", "sentrix_col"), sep = 3, remove = TRUE)
new_anno$sentrix_row = gsub("\\R","",new_anno$sentrix_row)
new_anno$sentrix_col = gsub("\\C","",new_anno$sentrix_col)
names(new_anno) <- c("Sample_Name", "Sex", "Slide", "sentrix_row", "sentrix_col", "Array", "Material", "Basename");
new_anno$Slide <- as.numeric(as.character(new_anno$Slide))
new_anno <- new_anno %>% filter(!is.na(Basename))
new_anno <- new_anno %>% filter(!duplicated(Basename))


#############  Quantile normalization 
### Part 2 ##  https://github.com/perishky/meffil/wiki/Functional-normalizing-separate-datasets   
#############  One of the slowest part ~ 3hrs for 13K for new set 100 samples ~20-40min
message("\nPart 2. ", Sys.time(), "\nBackground correction: meffil.qc ...\n");

new.qc.objects <- meffil.qc(new_anno, 
                          featureset = "common", 
                          cell.type.reference = NA, 
                          verbose = T)
saveRDS(new.qc.objects, file = file.path(batchdirout,"/new.qc.objects.rds"))
new.qc.objects   <-    readRDS(file.path(batchdirout,"/new.qc.objects.rds"))
gc()

## Read major reference set 13K samples
qc.objects <- readRDS(paste0(batchdir,"/qc.objects.13K.rds"))

## Pick some random 1K samples from reference set and combine two sets for normalization 
##combo.qc.objects   <- c(new.qc.objects, qc.objects[sample(1:length(names(qc.objects)), 1000, replace=T)]);
combo.qc.objects <- c(new.qc.objects, qc.objects);

combo.norm.objects <- meffil.normalize.quantiles( combo.qc.objects,
                       fixed.effects  = c("Material", "Array"), verbose = TRUE, 
                       random.effects = "Slide",    ## very slow part: on 13K set
                       number.pcs = 4 );            ## 23 may be way too much

black.cpgs    <- readRDS(paste0(batchdir,"/black.cpgs.rds"))
# Generating normalized reference beta to the number of PCA (#5 from step above) for later reuse with new batches. 
# Subsetting it back to reference data only (for number.pcs) checks (PC on ref has to be match with new batch)
ref.norm.objects <- combo.norm.objects[names(combo.norm.objects) %in% names(qc.objects)]
ref.norm.beta <- meffil.normalize.samples(ref.norm.objects ,
                                       just.beta = T,
                                       cpglist.remove = black.cpgs,
                                       verbose = TRUE)
saveRDS(ref.norm.beta, file = file.path(batchdir,"/norm.betas4se.rds"))
ref.norm.betas <- readRDS(paste0(batchdir,"/norm.betas4se.rds"))

# Subsetting it back to new data only 
new.norm.objects <- combo.norm.objects[names(combo.norm.objects) %in% names(new.qc.objects)]
new.norm.beta <- meffil.normalize.samples(new.norm.objects ,
                                       just.beta = T,
                                       cpglist.remove = black.cpgs,
                                       verbose = TRUE)


saveRDS(new.norm.beta, file = file.path(batchdirout,"/new.norm.beta4se.rds"))
new.norm.beta  <- readRDS(paste0(batchdirout,"/new.norm.beta4se.rds"))


#remove objects not needed for further analysis (i.e. free up memory)
rm(list=setdiff(ls(), c("new.norm.beta", "ref.norm.betas", "combo_samplesheet", "batchdir", "batchdirout", "cores")));
gc()


beta <- as.data.frame(cbind(new.norm.beta, ref.norm.betas[row.names(new.norm.beta),]))
rm(new.norm.beta, ref.norm.betas)
gc()


#############  Weeding extra XY/SNPs/etc probes based on 
### Part 3 ##  https://zwdzwd.github.io/InfiniumAnnotation   
#############  Another ~50K probes should be cleared.
message("\nPart 3.", Sys.time(), "\nProbes cleansing extra ...\n");

annotations_hg19 <- fread(file.path(batchdir,"HM450.hg19.manifest.tsv"),
                          sep="\t", 
                          verbose = TRUE,
                          data.table = FALSE,
                          stringsAsFactors = FALSE,
                          check.names = FALSE,
						  nThread = cores);

annotations_hg19 <- annotations_hg19 %>%
                    filter(MASK_general == FALSE) %>%
                    #filter(MASK_sub30_copy == FALSE) %>%
                    #filter(MASK_snp5_GMAF1p == FALSE) %>%
                    filter(!CpG_chrm == "chrX") %>%
                    filter(!CpG_chrm == "chrY") %>%
                    filter(!grepl("rs",probeID)) %>%
                    filter(!grepl("ch",probeID))

beta_filtered <- beta[rownames(beta) %in% annotations_hg19$probeID,]
rm(beta)
gc()

#reduce beta values by variable probes
beta_filtered <- t(beta_filtered)
beta_reduced <- beta_filtered[,order(-apply(beta_filtered,2,sd))[1:10000]]
beta_reduced <- as.data.frame(beta_reduced)
fwrite(beta_reduced, file.path(batchdirout,"combo_beta_top10Kprobes4se.tsv"),  sep = "\t", row.names = TRUE,  nThread = cores);
beta_reduced <- fread(file.path(batchdirout,"combo_beta_top10Kprobes4se.tsv"), sep = "\t", nThread = cores);
beta_reduced <- beta_reduced %>% column_to_rownames("V1");
rm(beta_filtered)
gc()


#############  Make dimension reduction with UWOT  
### Part 4 ##  PCA is slow single threaded need better package
#############  >5hrs
message("\nPart 4. Initiated:   ", Sys.time(), "\nGenerating PCA with top 200 componenets ...");
PC <- prcomp(beta_reduced, 
             center = TRUE, 
			 scale = FALSE,
			 rank. = 200); 

saveRDS(PC, paste0(batchdirout,"/combo_200pc_prcomp4se.rds"))
PC <- readRDS(paste0(batchdirout,"/combo_200pc_prcomp4se.rds"))
message("\nDone PCA:   ", Sys.time(), "\nGenerating Umap ...");

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
newanno <- combo_samplesheet[combo_samplesheet$idat_filename %in% umap$idat_filename,];
newanno$CNS.Subclass <- as.factor(newanno$CNS.Subclass);

message("\nDone unsupervised:   ", Sys.time(), "\n now supervised ...");
umap1 <- uwot::umap(PC$x,
                   #n_components = 3,
                   #pca = 100,
                   n_neighbors = 10,
                   y = anno$Combined_class_match,
                   spread = 2,
                   min_dist = 0.1,
                   local_connectivity = 1,
                   bandwidth = 1)
umap1 <- as.data.frame(umap1)
rownames(umap1) <- rownames(PC$x)
umap1 <- umap1 %>% rownames_to_column("idat_filename")

umap <- cbind(umap, umap1$V1, umap1$V2 );
names(umap) <- c("idat_filename", "x1","y1","x2","y2");

outf <- merge(umap, newanno[,c(
   "nn",
   "idat_filename",
   "Sample",
   "Sex",
   "Age",
   "material_prediction",
   "RFpurity.ABSOLUTE",
   "RFpurity.ESTIMATE",
   "LUMP",
   "Location_specific",
   "OS_months",
   "OS_status",
   "PFS_months",
   "PFS_status",
   "Histology",
   "Molecular",
   "Study",
   "CNS.MCF",
   "CNS.MCF.score",
   "CNS.Subclass",
   "CNS.Subclass.score")],
   by = "idat_filename");

outf <- arrange(outf, Study);

fwrite(outf, file.path(batchdirout,"/umap_2-PCA4se.txt"), 
    row.names=TRUE, sep = "\t", nThread = cores);

message(Sys.time());
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
message(" End of normalization and UMAP script ");
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");