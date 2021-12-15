## 12/14/2021
## Generates individual umaps in separate HTML files. 
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
      stop("!!! Crash !!!\nParameters missed in sbatch: [$1:outputfolder $2:newsampleslist] \n\n")
}else{
     batchdirout = as.character(args[1]);
     newsamplesheetname = as.character(args[2]);
     message("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
		 Sys.time(),"  Started with parameters:\nDump dir:\t", 
		 batchdirout, "\nNew samples:\t", newsamplesheetname,
		 "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
}
options(scipen = 999, "width"=180, mc.cores=54);
cores <- options()$mc.cores
set.seed(1234)

suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(tidyr))
suppressMessages(library(openxlsx))
suppressMessages(library(meffil))
suppressMessages(library(dplyr))
# Bunch for visualization/dimension reduction 
suppressMessages(library(tibble))
suppressMessages(library(uwot))
suppressMessages(library(ggplot2))
suppressMessages(library(Rtsne))
suppressMessages(library(gridExtra))
suppressMessages(library(tidyverse))
suppressMessages(library(densvis))

## 
idatrootfolder = "/data/MDATA/compass/iScan_raw";
batchdir      <- "/data/MDATA/NormRcode/CNS13Krefset"; 
batchdir1     <- "/data/MDATA/NormRcode/FixedCNSset";
samplesheet   <- "/data/MDATA/NormRcode/CNS13Krefset/Compass_15K.xlsx";
# batchdirout    <- "/data/MDATA/TRANSFER/SAMPLESHEETS/WeeklyCNS"; 
# newsamplesheetname <- "/data/MDATA/TRANSFER/HTMLS/dna_newsamples_Dec10.txt"; 

if(newsamplesheetname == "/data/MDATA/NormRcode/CNS13Krefset/newsamples.csv"){
	newsamplesheet <- read.csv(newsamplesheetname);
	message("\n !!! Running on premade samplesheet !!!\n");
}else{
	newsamplesheet <- read.csv(newsamplesheetname, header=F, sep="\t");
	message("\n Loaded DNA and sentrix IDs from the text file.\n");
	names(newsamplesheet) <- c("Sample", "idat");
    newsamplesheet$Sentrix_ID       = NA;
	newsamplesheet$Sentrix_position = NA;
    newsamplesheet$Platform_methy   = NA;
	for(i in 1:length(newsamplesheet$idat)){
	  newsamplesheet$Sentrix_ID[i]       = head(tail(unlist(str_split(newsamplesheet$idat[i], "_")), n=2), n=1);
	  newsamplesheet$Sentrix_position[i] = tail(unlist(str_split(newsamplesheet$idat[i], "_")), n=1);
	  greenfile = list.files(path = paste0(idatrootfolder,"/",newsamplesheet$Sentrix_ID[i]), pattern = "_Grn.idat.gz$|_Grn.idat$")[1];
	  greensize = file.info(paste0(idatrootfolder,"/",newsamplesheet$Sentrix_ID[i],"/",greenfile))$size;
	  if(greensize < 5000000){
        newsamplesheet$Platform_methy[i] = "HumanMethylation450";
	  }else{
        newsamplesheet$Platform_methy[i] = "MethylationEPIC";
	  }
	}
    newsamplesheet$Material_Type[1:nrow(newsamplesheet)] <- "FFPE"    
    newsamplesheet$Basename = paste0(idatrootfolder,"/",newsamplesheet$Sentrix_ID,"/", newsamplesheet$idat)

## Some bunch for classifier
	suppressMessages(library(mnp.v11b6))
	suppressMessages(library(RFpurify))
	suppressMessages(library(LUMP))
	source("/app/R/cns_classifier.vNCI.R")	         #"classifier" "get_predicted"
	source("/app/R/predict_umap.R")                  #"predict_umap"
	source("/app/R/process_prediction_result.R")	 #"get_class_score

    RGset = read.metharray(newsamplesheet$Basename, force = T, verbose=T);
    GMsetEx <- mapToGenome(RGset);
	estSex  <- getSex(GMsetEx);
    newsamplesheet$Sex <- estSex$predictedSex;

	###Run DKFZ classifier
	cat("Run DKFZ classifier...\n");
	CNSclassifier = lapply(1:length(colnames(RGset)), function(x) classifier(RGset, newsamplesheet$idat[x], newsamplesheet$Material_Type[x]));
	CNSclassifier = data.frame(do.call(rbind, CNSclassifier));
	colnames(CNSclassifier) = c("ID",
								"predFFPE",
								"CNS.MCF",
								"CNS.Subclass", 
								"RFpurity.ABSOLUTE",
								"RFpurity.ESTIMATE",
								"LUMP");
    rm(RGset)
    gc()
   
    newsamplesheet = cbind(newsamplesheet, CNSclassifier);
	newsamplesheet = separate(newsamplesheet, CNS.MCF, into = c("CNS.MCF", "CNS.MCF.score"), sep = ":", remove = TRUE);
   	newsamplesheet = separate(newsamplesheet, CNS.Subclass, into = c("CNS.Subclass", "CNS.Subclass.score"), sep = ":", remove = TRUE);
    newsamplesheet$CNS.MCF.score      = as.numeric(gsub(";.+", "", newsamplesheet$CNS.MCF.score, perl=T));
	newsamplesheet$CNS.Subclass.score = as.numeric(gsub(";.+", "", newsamplesheet$CNS.Subclass.score, perl=T));
    
    tempsheet <- data.frame(
		nn = 1:nrow(newsamplesheet),
		Sample = newsamplesheet$Sample,
		idat_filename = newsamplesheet$idat,
		idat = newsamplesheet$ID,
		Sentrix_ID =  newsamplesheet$Sentrix_ID,
		Sentrix_position =  newsamplesheet$Sentrix_position,
		material_prediction = newsamplesheet$Material_Type,
		Platform_methy = newsamplesheet$Platform_methy,
		Age = rep("", nrow(newsamplesheet)),
		Sex = newsamplesheet$Sex, 
		matched_cases     = rep("", nrow(newsamplesheet)),
		Location_general  = rep("", nrow(newsamplesheet)),
		Location_specific = rep("", nrow(newsamplesheet)),
		Neoplastic = rep("Neuropathology", nrow(newsamplesheet)),
		Primary_category = rep("Case", nrow(newsamplesheet)),
		CNS_study  =  rep("Case", nrow(newsamplesheet)),
		OS_months  =  rep("", nrow(newsamplesheet)),
		OS_status  =  rep("", nrow(newsamplesheet)),
		PFS_months =  rep("", nrow(newsamplesheet)),
		PFS_status =  rep("", nrow(newsamplesheet)),
		Histology  =  rep("", nrow(newsamplesheet)),
		Molecular  =  rep("", nrow(newsamplesheet)),
		Study      =  rep("", nrow(newsamplesheet)),
		predFFPE   =  newsamplesheet$predFFPE,
		CNS.MCF    =  newsamplesheet$CNS.MCF,
		CNS.MCF.score =  newsamplesheet$CNS.MCF.score,
		CNS.Subclass  =  newsamplesheet$CNS.Subclass,
		CNS.Subclass.score= newsamplesheet$CNS.Subclass.score,
		RFpurity.ABSOLUTE = newsamplesheet$RFpurity.ABSOLUTE,
		RFpurity.ESTIMATE = newsamplesheet$RFpurity.ESTIMATE,
		LUMP              = newsamplesheet$LUMP,
		NCI_METRIC        = newsamplesheet$CNS.Subclass,
		Basename          = newsamplesheet$Basename)

      newsamplesheet <- tempsheet;
	  write.csv(newsamplesheet, file.path(batchdirout, "newsamplesheet.csv"));
}

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

new.qc.objects <- meffil.qc(new_anno, 
                          featureset = "common", 
                          cell.type.reference = NA, 
                          verbose = T);
#saveRDS(new.qc.objects, file.path(batchdirout, "new.qc.objects.rds"));

message("\n Mix new and reference data for quantile normalization\n");
qc.objects <- readRDS(file.path(batchdir1, "qc.objects.withoutbads.11K.rds"));
bad.cpgs   <- readRDS(file.path(batchdir1, "bad.cpgs.F.rds"));

combo.qc.objects <- c(new.qc.objects, qc.objects);
combo.norm.objects <- meffil.normalize.quantiles( combo.qc.objects,
                       fixed.effects  = "Array", verbose = TRUE, 
                       number.pcs = 14 ); 

options(mc.cores=12)  
cores <- options()$mc.cores
gc();

message("\n Cores number changed to: ",  cores, "\n Extracting new betas...\n")

new.norm.objects <- combo.norm.objects[names(combo.norm.objects) %in% names(new.qc.objects)]
new.norm.objects <- new.norm.objects[!duplicated(names(new.norm.objects))]
new.norm.beta    <- meffil.normalize.samples(new.norm.objects,
                                       just.beta = T,
                                       cpglist.remove = bad.cpgs,
                                       verbose = FALSE)
rm(new.norm.objects, new.qc.objects, combo.qc.objects, new_anno, bad.cpgs)
gc()

message("\n Blending betas from reference and new datasets.\n");
ref.norm.betas <- fread(file.path(batchdir1,"beta_top027SDprobes_pc14.tsv"), sep = "\t", nThread = cores);
ref.norm.betas <- ref.norm.betas %>% column_to_rownames("V1");
new.norm.beta  <- as.data.frame(t(new.norm.beta))
new.norm.beta  <- new.norm.beta[,colnames(ref.norm.betas)]
beta <- rbind(ref.norm.betas, new.norm.beta)
rm(ref.norm.betas, new.norm.beta)
gc()


#fwrite(beta, file.path(batchdirout,"beta_combined.txt"), row.names=TRUE, sep = "\t", nThread = cores);


message( Sys.time(), "\n Principal component reduction on reduced beta matrix.\n");
## Faster PC function
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

PC <- prcomp_svds(beta, k=150); 
saveRDS(PC, file.path(batchdirout,"combo_150pc_prcomp14.rds"));
PC <- readRDS(file.path(batchdirout,"combo_150pc_prcomp14.rds"));

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

umap <- as.data.frame(umap)
rownames(umap) <- rownames(PC$x)
umap <- cbind(umap, dmap)
umap <- umap %>% rownames_to_column("idat_filename")
names(umap) <- c("idat_filename",  "x1","y1",  "x2","y2" );


# Blending metadata with umap
anno_base <- openxlsx::read.xlsx(samplesheet);
combo_samplesheet <- rbind(anno_base, newsamplesheet[,-33]);
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

outf <- outf[!duplicated(outf$Sample),]
outf[(outf$Sample %in% newsamplesheet$Sample), "Study"] <- "compass"

## sanitizing labels ##
outf$CNS.MCF      <- gsub(" |,|/|__", "_", outf$CNS.MCF) 
outf$CNS.MCF      <- gsub("__", "_", outf$CNS.MCF) 
outf$CNS.MCF      <- gsub("__", "_", outf$CNS.MCF) 
outf$CNS.Subclass <- gsub(" |,|/|__", "_", outf$CNS.Subclass ) 
outf$CNS.Subclass <- gsub("__", "_", outf$CNS.Subclass)
outf$CNS.Subclass <- gsub("__", "_", outf$CNS.Subclass)
outf$CNS.Subclass <- gsub("_$", "",  outf$CNS.Subclass, perl=TRUE)
outf$NCI_METRIC   <- gsub(" |,|/|__", "_", outf$NCI_METRIC)
outf$NCI_METRIC   <- gsub("__", "_", outf$NCI_METRIC)
outf$NCI_METRIC   <- gsub("__", "_", outf$NCI_METRIC)
outf$NCI_METRIC   <- gsub("_$", "",  outf$NCI_METRIC, perl=TRUE)

outf <- arrange(outf, Study);

fwrite(outf, file.path(batchdirout,"umap_2-PCA5.txt"), row.names=TRUE, sep = "\t", nThread = cores);
#outf <-read.csv( file.path(batchdirout,"umap_2-PCA5.txt"), sep="\t", row.names=1) 

## Prepare distances (for individual boxplots) 
X=aggregate(outf[,"x1"],  list(outf$NCI_METRIC), median)
Y=aggregate(outf[,"y1"],  list(outf$NCI_METRIC), median)
centers <- merge(X,Y, by = "Group.1")
names(centers) <- c("LABEL", "x1", "y1")

message("\nComputing distances to the centers of all clusters...");
distances=c();
for (n in 1:nrow(outf)){
 X=outf[n,"x1"];
 Y=outf[n,"y1"];
 distout = c(); 
   for (k in centers$LABEL){
   	 dist <- sqrt((centers[centers$LABEL == k, "x1"] - outf[n,"x1"]) ^ 2 
		    + (centers[centers$LABEL == k, "y1"] - outf[n,"y1"]) ^ 2  )
	 distout <- c(distout, dist);
   }
 distances <- cbind(distances, distout);
}
row.names(distances)<- centers$LABEL;
colnames(distances) <- outf$Sample;
distances<- distances[(row.names(distances)!=""),]

## General html with all labelled (new) samples
#tstamp <- format(Sys.time(), "%a%d%b%Y");
tstamp <- format(Sys.Date(), "%m%d%Y");


rmarkdown::render(file.path(batchdir1, "plot_CNS_umap_all.Rmd"), 
    output_file = file.path(batchdirout,paste0("plot_CNS_umap_",tstamp,".html")));


## save html for individual samples
for(dna in newsamplesheet$Sample){
	print(dna)
    result <- outf %>% filter(Study  != "compass") 
    dnarow <- outf %>% filter(Sample == dna)
	result <- rbind(dnarow[1,], result)
	result$Study[1] <- "compass";

	fwrite(result, file.path(batchdirout,"umap_2-PCA5x1.txt"), 
    row.names=TRUE, sep = "\t", nThread = cores);    
#  
    rmarkdown::render(file.path(batchdir1, "plot_CNS_umap_x1.Rmd"), 
      output_file = paste0(batchdirout,"/", dna ,"_umap_",tstamp,".html"));	
}

message(Sys.time());
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
message(" *", tstamp, ".html files are ready for collection: ");
message(" ",  batchdirout);
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
