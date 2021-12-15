## Drews normilization pipeline with small tweaks by Rust
## 11/26/2021:: With two UMAP plots on output 
options(scipen = 999, "width"=180)
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

options(mc.cores=140) ## 140/54
cores <- options()$mc.cores

## compass dir
idatrootfolder = "compass/iScan_raw";
batchdir      <- "NormRcode/CNS13Krefset"; 
batchdirout   <- "TRANSFER/SAMPLESHEETS/compass"; 
samplesheet   <- "NormRcode/CNS13Krefset/Compass_15K.xlsx";

#############  Inject static samplesheet
### Part 0 ##  anno_base -- rich file for the ref. set
#############  anno -- minimal for meffil normalization

anno_base <- openxlsx::read.xlsx(samplesheet)
anno <- anno_base[,c( "idat_filename", "Sex", "Sentrix_ID", "Sentrix_position", "Platform_methy","material_prediction")];
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

### New samples sheet can be prepared with samplegenerator.R 
newsamplesheet <- read.csv("/data/MDATA/NormRcode/CNS13Krefset/newsamples.csv");
newsamplesheet$CNS.MCF      <- gsub(" |,|/|__", "_", newsamplesheet$CNS.MCF); 
newsamplesheet$CNS.MCF      <- gsub("__", "_",       newsamplesheet$CNS.MCF); 
newsamplesheet$CNS.MCF      <- gsub("__", "_",       newsamplesheet$CNS.MCF); 
newsamplesheet$CNS.Subclass <- gsub(" |,|/|__", "_", newsamplesheet$CNS.Subclass); 
newsamplesheet$CNS.Subclass <- gsub("__", "_",       newsamplesheet$CNS.Subclass);
newsamplesheet$CNS.Subclass <- gsub("__", "_",       newsamplesheet$CNS.Subclass);
newsamplesheet$CNS.Subclass <- gsub("_$", "",        newsamplesheet$CNS.Subclass, perl=TRUE);
newsamplesheet$NCI_METRIC   <- gsub(" |,|/|__", "_", newsamplesheet$NCI_METRIC);
newsamplesheet$NCI_METRIC   <- gsub("__", "_",     newsamplesheet$NCI_METRIC);
newsamplesheet$NCI_METRIC   <- gsub("__", "_",     newsamplesheet$NCI_METRIC);
newsamplesheet$NCI_METRIC   <- gsub("_$", "",      newsamplesheet$NCI_METRIC, perl=TRUE);

## Putting both samplesheets together without Basename column 
combo_samplesheet <- rbind(newsamplesheet[,-33], anno_base);

fwrite(combo_samplesheet, file.path(batchdirout, "combo_samplesheet.tsv"),
       sep = "\t", row.names = TRUE,  nThread = cores);

## Reshaping of truncated rawsamplesheet with new annotation for meffil format
#new_anno <- newsamplesheet %>% filter( CNS_study == "Case" );
new_anno <- newsamplesheet[,c( "idat_filename", "Sex", "Sentrix_ID", "Sentrix_position", "Platform_methy","material_prediction", "Basename")];
new_anno$material_prediction[new_anno$material_prediction!="FFPE"]<-"Frozen";
new_anno$Sex[!(new_anno$Sex %in% c("M","F"))]<- NA;
new_anno <- separate(new_anno, Sentrix_position, into = c("sentrix_row", "sentrix_col"), sep = 3, remove = TRUE)
new_anno$sentrix_row = gsub("\\R", "", new_anno$sentrix_row)
new_anno$sentrix_col = gsub("\\C", "", new_anno$sentrix_col)
names(new_anno) <- c("Sample_Name", "Sex", "Slide", "sentrix_row", "sentrix_col", "Array", "Material", "Basename");
new_anno$Slide <- as.numeric(as.character(new_anno$Slide))
new_anno <- new_anno %>% filter(!is.na(Basename))
new_anno <- new_anno %>% filter(!duplicated(Basename))


#############  Quantile normalization 
### Part 2 ##  https://github.com/perishky/meffil/wiki/Functional-normalizing-separate-datasets   
#############  One of the slowest part ~ 3hrs for 13K for new set 100 samples ~20-40min
#message("\nPart 2. ", Sys.time(), "\nBackground correction on new set: meffil.qc ...\n");
new.qc.objects <- meffil.qc(new_anno, 
                          featureset = "common", 
                          cell.type.reference = NA, 
                          verbose = T);
saveRDS(new.qc.objects, file = file.path(batchdirout,"new.qc.objects.rds"))
new.qc.objects   <-  readRDS(file.path(batchdirout,"new.qc.objects.rds"))
gc()


## Make HTML report for new samples
new.qc.summary <- meffil.qc.summary(new.qc.objects);
meffil.qc.report(new.qc.summary, output.file=file.path(batchdirout,"qc_report_new.html"))

## Generate  major reference set 13K samples
##  Has to be done once a month when toped up with new samples other times just comment 
#   message("\nPart 2 continue. ", Sys.time(), "\nBackground correction on the reference set: meffil.qc ...\n");
#	qc.objects <- meffil.qc(anno, 
#							featureset = "common", 
#							cell.type.reference = NA, 
#							verbose = T);
#	saveRDS(qc.objects, file = file.path(batchdir,"qc.objects.15K.rds"))
#	gc()
### Read major reference set 13K samples
qc.objects <- readRDS(file.path(batchdir,"qc.objects.15K.rds"))
qc.objects <- qc.objects[!names(qc.objects)=="Error"];  # We may loose quite a number of innocent samples here.
gc();

#############################################################################
#  Repeat injection for dropped /missed /new samples
lostsample <- setdiff(anno$Sample_Name, names(qc.objects))
print(lostsample); 
bigpilelst <- list.files("idat/", pattern = "_Grn.idat.gz", ignore.case = T )
bigpilelst <- gsub("_Grn.idat.gz", "", bigpilelst)
if(length(setdiff(lostsample, bigpilelst))>1){
	  print(setdiff(lostsample, bigpilelst));
	  stop("Sampes above no files here: /data/MDATA/idat\n Died!");
}else{print("OK no idat files missed, topping up the reference set.");
    if(length(lostsample > 0)){
		lost.qc.objects <- meffil.qc(anno[anno$Sample_Name %in% lostsample, ], 
									featureset = "common", 
									cell.type.reference = NA, 
									verbose = T);
		lost.qc.objects <- lost.qc.objects[!names(lost.qc.objects)=="Error"]
		tryCatch(
		  {
			 message("\ncombine reference qc objects ...\n")
			 qc.objects <- do.call(c, list(qc.objects, lost.qc.objects))
		  },
			 error = function(e){  NULL  }  
		)
	}else{
			 message("\nwill continue filtering any failed samples with error flag.\n")
	}
};

## Any more ?
qc.objects <- qc.objects[!names(qc.objects)=="Error"]; 
lostsample <- setdiff(anno$Sample_Name, names(qc.objects))
print(c("Idat are still needed for:", lostsample));

saveRDS(qc.objects, file = file.path(batchdir,"qc.objects.15K.rds"))

#############  Now clearing up reference set (only)
### Part 3 ##  This part is not too bad in time but it is better run 
#############  on each batch rather than keep track all bad things within reference set
#DREW: this section outlines what Drew consider to be the important qc checks for sample quality
#outlier samples whose predicted median methylated signal is more than 3 standard deviations from the expected
#this section outlines what I consider to be the important qc checks for sample quality
#outlier samples whose predicted median methylated signal is more than 3 standard deviations from the expected
 message("\n\nPart 3: started quantiles normalization:   ", Sys.time(), "\n");
     qc.summary <- meffil.qc.summary(qc.objects); # ~20-30min 
     saveRDS(qc.summary, file = file.path(batchdir,"qc.summary.rds"))
     meffil.qc.report(qc.summary, output.file=file.path(batchdir,"qc_report.html"))

qc.summary <- readRDS(file.path(batchdir,"qc.summary.rds"))
meffil.qc.report(qc.summary, output.file=file.path(batchdir,"qc_report.html"));
bad.cpgs   <- qc.summary$bad.cpgs$name;
saveRDS(bad.cpgs, file = file.path(batchdir,"bad.cpgs.rds"))

### Never turn off this line:            !!!!!!!!!!!!!!
bad.cpgs   <- readRDS(file.path(batchdir,"bad.cpgs.rds"))

bad_samples_meth_unmeth <- as.data.frame(qc.summary$meth.unmeth.summary$tab)
bad_samples_meth_unmeth <- bad_samples_meth_unmeth[bad_samples_meth_unmeth$outliers==TRUE,]
anno_bad_meth_unmeth    <- newsamplesheet[newsamplesheet$idat_filename %in% bad_samples_meth_unmeth$sample.name,]
bad_samples_meth_unmeth <- as.data.frame(qc.summary$meth.unmeth.summary$tab)
bad_samples_meth_unmeth <- bad_samples_meth_unmeth[bad_samples_meth_unmeth$outliers==TRUE,]
anno_bad_meth_unmeth    <- newsamplesheet[newsamplesheet$idat_filename %in% bad_samples_meth_unmeth$sample.name,]

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
anno_bad_detectionp <- newsamplesheet[newsamplesheet$idat_filename %in% bad_samples_detectionp$sample.name,]

#samples with proportion of probes with bead number < 3 is > 0.2
bad_samples_lowbead <- as.data.frame(qc.summary$sample.beadnum.summary$tab)
bad_samples_lowbead <- bad_samples_lowbead[bad_samples_lowbead$outliers==TRUE,]
anno_bad_lowbead <- newsamplesheet[newsamplesheet$idat_filename %in% bad_samples_lowbead$sample.name,]

####combine all bad samples into vector ~144
bad_samples <- c(
  bad_samples_meth_unmeth$sample.name,
  bad_samples_control$sample.name,
  bad_samples_detectionp$sample.name,
  bad_samples_lowbead$sample.name)
write.csv(bad_samples, file.path(batchdir,"bad_samples_qc_failed.csv"));

####remove ~150 bad samples prior to performing quantile normalization
qc.objects <- meffil.remove.samples(qc.objects, bad_samples)
saveRDS(qc.objects, file = file.path(batchdir,"qc.objects.withoutbads.15K.rds"))
qc.objects <- readRDS(file.path(batchdir,"qc.objects.withoutbads.15K.rds"))

### Optional part 
 message("\nGenerating PDF with principal components plot:\n")
	pc.fit <- meffil.plot.pc.fit(qc.objects)
	pdf(file.path(batchdir, "PC.rplot.pdf")) 
	print(pc.fit$plot)
	dev.off()


## Combine new batch and reference set 
combo.qc.objects <- c(new.qc.objects, qc.objects);

options(mc.cores=32)  ## 32 for large memory 
cores <- options()$mc.cores
gc();
message("\nCores number changed to: ",  cores, "\n\n")
combo.norm.objects <- meffil.normalize.quantiles( combo.qc.objects,
                       fixed.effects  = "Array", verbose = TRUE, 
                       number.pcs = 14 ); 
#                      random.effects = "Slide",      ## may not work properly: on 13K set

message("\n\nDone quantiles normalization:   ", Sys.time(), "\n now extracting beta...");

## Generating normalized reference beta to the number of PCA (#4 from step above) for later reuse with new batches. 
## Subsetting it back to reference data only (for number.pcs) checks (PC on ref has to be match with new batch)
	ref.norm.objects <- combo.norm.objects[names(combo.norm.objects) %in% names(qc.objects)]
	ref.norm.beta    <- meffil.normalize.samples(ref.norm.objects ,
										   just.beta = T,
										   cpglist.remove = bad.cpgs,
										   verbose = TRUE)
 saveRDS(ref.norm.beta, file = file.path(batchdir,"/norm.betas_pc14.rds"))
ref.norm.betas <- readRDS(paste0(batchdir,"/norm.betas_pc14.rds"))

# Subsetting it back to new data only 
message("\n\nDone quantiles normalization for reference set:   ", Sys.time(), "\n now for new batch...");

new.norm.objects <- combo.norm.objects[names(combo.norm.objects) %in% names(new.qc.objects)]
new.norm.objects <- new.norm.objects[!duplicated(names(new.norm.objects))]
new.norm.beta    <- meffil.normalize.samples(new.norm.objects,
                                       just.beta = T,
                                       cpglist.remove = bad.cpgs,
                                       verbose = FALSE)

saveRDS(new.norm.beta, file = file.path(batchdirout,"new.norm.beta_pc14.rds"))
new.norm.beta  <- readRDS(file.path(batchdirout,"new.norm.beta_pc14.rds"))

#remove objects not needed for further analysis (i.e. free up memory)
rm(list=setdiff(ls(), c("new.norm.beta", "ref.norm.betas", "combo_samplesheet", "batchdir", "batchdirout", "cores")));
gc()

beta <- as.data.frame(cbind(new.norm.beta, ref.norm.betas[row.names(new.norm.beta),]))
rm(new.norm.beta, ref.norm.betas)
gc()

#############  Weeding extra XY/SNPs/etc (~50k) probes based on 
### Part 4 ##  https://zwdzwd.github.io/InfiniumAnnotation   
#############  then reducing matrix to most to top most variable probes 
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
#beta_reduced  <- beta_filtered[,order(-apply(beta_filtered,2,sd))[1:10000]]
beta_reduced <- beta_filtered[,Rfast::colVars(as.matrix(beta_filtered), std = TRUE, parallel = TRUE)>0.227]; 
write.table(beta_reduced, file.path(batchdirout, "combo_beta_top10Kprobes_pc14.tsv"),  sep = "\t", row.names = TRUE);
beta_reduced  <- fread(file.path(batchdirout,"combo_beta_top10Kprobes_pc14.tsv"), sep = "\t", nThread = cores);
beta_reduced  <- beta_reduced[!duplicated(beta_reduced$V1),];
beta_reduced  <- beta_reduced %>% column_to_rownames("V1");
rm(beta_filtered)
gc()


#############  Make dimension reduction with UWOT  
### Part 5 ##  PCA is slow / single threaded need better package
#############  >2hrs
message("\nPart 4. Initiated:   ", Sys.time(), "\nGenerating PCA with top 150 componenets ...");
#PC <- prcomp(as.matrix(beta_reduced), 
#            center = TRUE, 
#			 scale = FALSE,
#			 rank. = 200); 

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

PC <- prcomp_svds(beta_reduced, k=150) 
saveRDS(PC, file.path(batchdirout,"combo_150pc_prcomp14.rds"));
PC <- readRDS(file.path(batchdirout,"combo_150pc_prcomp14.rds"));

newanno <- combo_samplesheet[combo_samplesheet$idat_filename %in% rownames(PC$x),];

message("\nDone PCA:   ", Sys.time(), "\n now embedding ...");

dmap <- densmap(
        PC$x,
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
message("\nDone densmap:   ", Sys.time(), "\nGenerating Umap ...\n");

umap <- uwot::umap(PC$x,
                   #n_components = 2,
                   #pca = 25,
                   n_neighbors = 10,
                   metric = "cosine",
                   #y =  newanno$CNS.Subclass,
                   spread = 1,
                   min_dist = 0,
                   local_connectivity = 1,
                   bandwidth = 1,
				   fast_sgd = FALSE,
				   pcg_rand = FALSE,
				   n_sgd_threads =1,
				   approx_pow = FALSE
				   );
# It's safer to leave fast_sgd = FALSE. If fast_sgd = TRUE, then user-supplied values of 
# pcg_rand, n_sgd_threads, and approx_pow are ignored.


##  Supervised UMAP (combo_samplesheet
X  <- data.frame(index=1:length(rownames(PC$x)), idat_filename=rownames(PC$x))
XX <- merge(X, newanno, by="idat_filename")  ## Scuko !!!
newanno <- XX[order(X$index),-1]
mlabels <- as.numeric(as.factor(newanno$CNS.Subclass));
mlabels[is.na(mlabels)] <- 0;
umap1 <- uwot::umap(PC$x,
                   #n_components = 3,
                   #pca = 100,
                   n_neighbors = 10,
                   y = mlabels,
                   spread = 1,
                   min_dist = 0,
                   local_connectivity = 1,
                   bandwidth = 1,
				   fast_sgd = FALSE,
				   pcg_rand = FALSE,
				   n_sgd_threads =1,
				   approx_pow = FALSE)

umap <- as.data.frame(umap)
rownames(umap) <- rownames(PC$x)
rownames(umap) <- make.names(rownames(PC$x), unique=TRUE)
umap <- cbind(umap, dmap, umap1)
umap <- umap %>% rownames_to_column("idat_filename")
names(umap) <- c("idat_filename",  "x1","y1",  "x2","y2",  "x3","y3");
umap$idat_filename <- gsub("^X", "", umap$idat_filename, perl=TRUE)


newanno <- combo_samplesheet[combo_samplesheet$idat_filename %in% rownames(PC$x),];
outf <- merge(umap, newanno[,c(
   "idat_filename",
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
   "Molecular",
   "Study",
   "CNS.MCF",
   "CNS.MCF.score",
   "CNS.Subclass",
   "CNS.Subclass.score",
   "NCI_METRIC")],
   by = "idat_filename");

outf <- arrange(outf, Study);

fwrite(outf, file.path(batchdirout,"umap_2-PCA5.txt"), 
    row.names=TRUE, sep = "\t", nThread = cores);

tstamp <- format(Sys.time(), "%a%d%b%Y");
rmarkdown::render(file.path(batchdir, "plot_CNS_umap_X.Rmd"), 
    output_file = paste0("/data/MDATA/", batchdirout, "/plot_CNS_umap_",tstamp,".html"));

message(Sys.time());
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
message(" End of normalization and UMAP script ");
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");