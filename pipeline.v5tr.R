#!/usr/bin/Rscript
##Code use DKFZ classifier mnp.v11b6 to predict samples profiled by EPIC/450k array
##Zhichao's code dockerized & upgraded for R version 4.0.3 by Rust [ 6/3/2021 ]

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 1){
      stop("Usage: [idat.folder] \nThis is only centrix barcode from location:\n /data/MDATA/compass/iScan_raw ")
}

centrix = as.character(args[1])
idat_path = paste0("compass/iScan_raw/",centrix)
csv_file_name = "Sample_Sheet.csv";
output = "compass/ClassifierReports";
#csv_file_name = as.character(args[3])
setwd("/mnt");

#####################################################################3
paired.color = RColorBrewer::brewer.pal(12,"Paired")
library(rmarkdown)
library(caret)
library(dplyr)
library(plotly)
library(RANN)
library(knitr)
library(kableExtra)
library(IlluminaHumanMethylationEPICmanifest)          #CNV
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #CNV
suppressMessages(library(mnp.v11b6))
#suppressMessages(library(mnp.v11b4))
library(RFpurify)
library(LUMP)
#source("R/cns_classifier.R")
source("/app/R/cns_classifier.vNCI.R")	         #"classifier" "get_predicted"
source("/app/R/processBeta.pca.v1.2.R")          #"pca.center" "pca.rotation"  "pca.sd" "pca_transform" "getBetas32k" "x.ref" "y.ref"
source("/app/R/predict_umap.R")                  #"predict_umap"
source("/app/R/run_tsne2d.R")                    #"run_tsne"
source("/app/R/cancer_class_color.R")            #Set of cancer classes with associated color 
source("/app/R/queryNeighbors.R")                #"knn.search"
source("/app/R/process_prediction_result.R")	 #"get_class_score"  "get_family_name"  "get_family_score"
rmd_file = "/app/R/Generate_HTMLreport.v5.Rmd"   #"rmd_file" just an adress to the kable template
source("/app/R/MNPcnvplot1.R")                   #"MNPcnvplot1" #the copy number bar
#source("R/dend.sub.R")
#load models
load("/app/R/CNS.KNN.v1.1.Rda")                  #"cns.knn.v1.1" #15-nearest neighbor model
load("/app/R/Ref_distance.rda")                  #"ref_distance.subclass"
#library(ComplexHeatmap)
#####################################################################
#####################################################################


targets <- read.csv(file.path(idat_path, csv_file_name), stringsAsFactors = FALSE, skip = 7)
#targets <- targets[targets$Sample_Plate == 'Clinical'|targets$Sample_Plate == 'Sarcoma Project',]; ## Only known CNS samplses 'Sarcoma Project' later
if(nrow(targets) < 2){
	message(file.path(idat_path, csv_file_name));
	stop("ERROR! There are no appropriate samples for classifier on this chip.\n Or may be just a single sample. Check the input file above.\n");
};

targets$idat <- paste0(targets$Sentrix_ID,"_",targets$Sentrix_Position);
targets$Basename <- file.path(idat_path, targets$idat);
targets$Outname  <- file.path("compass/ClassifierReports", paste0(targets$Sample_Name, "_", targets$idat));

print(head(targets, 2))
RGset = read.metharray.exp(targets = targets, verbose=T, force=TRUE);
print(colnames(RGset))
###Run DKFZ classifier
cat("Run DKFZ classifier...\n");
#CNSclassifier = data.frame(t(sapply(1:nrow(targets), function(x) classifier(targets$Basename[x], targets$Material_Type[x]))))
targets$Material_Type[targets$Material_Type!="FFPE"]<-"Frozen";
CNSclassifier = lapply(1:nrow(targets), function(x) classifier(RGset, targets$idat[x], targets$Material_Type[x]));
CNSclassifier = data.frame(do.call(rbind, CNSclassifier));
#p = c(arrayID, p.ffpe, p1, p2, p.rf.ABSOLUTE, p.rf.ESTIMATE, p.lump)
colnames(CNSclassifier) = c("ID","predFFPE","CNS.MCF","CNS.Subclass", 
							"RFpurity.ABSOLUTE", "RFpurity.ESTIMATE", "LUMP");

CNSclassifier = cbind(CNSclassifier, targets[,c(1:(ncol(targets)))])
#write.csv(CNSclassifier, file = paste0(output,"/KNN.classifier.csv"), row.names=F) ## full table
for(i in 1:nrow(CNSclassifier)) {
   dir.create(CNSclassifier$Outname[i]);
   write.csv(CNSclassifier[i,], file = paste0(CNSclassifier$Outname[i],"/", CNSclassifier$Sample_Name[i], ".KNN.csv"), row.names=F);
   message( CNSclassifier$Sample_Name[i], " KNN report saved");
}
message("Moving to generation TSNE plots in htmls")

#get beta values
betas = getBetas32k(RGset, targets);

#get pcas
cat("Get PCA transformation...\n")
pca.res = pca_transform(betas)
x.ref = pca.res[[1]]
y.ref = pca.res[[2]]
x.test = pca.res[[3]]
rm(pca.res) 
if(sum(!y.ref %in% names(cancer.colour))> 0 ){ stop("Some cancer type not find matched color:", y.ref[!y.ref %in% names(cancer.colour)]); }
if(sum( !names(ref_distance.subclass) %in% y.ref) > 0){ stop("Some cancer distance density not found: ",names(ref_distance.subclass)[! names(ref_distance.subclass) %in% y.ref]); }

####TSNE
cat("TSNE...\n");
tsne_s = run_tsne(x.ref, x.test, y.ref, n.pc = 300)
ref.tsne = tsne_s[[1]] 
test.tsne = tsne_s[[2]]
ref.tsne.center = tsne_s[[3]]
rm(tsne_s)
#save(ref.tsne, test.tsne, ref.tsne.center, file = paste0(output, ".tsne.rda"))


####UMAP
cat("UMAP...\n");
umap_s = predict_umap(x.test[,1:94])
ref.umap = umap_s[[1]]
test.umap = umap_s[[2]]
ref.umap.center = umap_s[[3]]
rm(umap_s)
y.test = targets$Sample_Name
#save(ref.umap, test.umap, ref.umap.center, y.ref, y.test, file = paste0(output, ".umap.rda"))


###
#KNN predict on new test
cat("KNN classifer...\n");
knn.p.raw = predict(cns.knn.v1.1, x.test, type = "prob")

#calibration
knn.p.calibrated = predict(knn.score.calibration$glmnet.fit, newx = knn.p.raw,
						   type = "response", s = knn.score.calibration$lambda.1se)[,,1]
knn.p.calibrated.MCF = get_family_score(knn.p.calibrated) #not very good 
#save(knn.p.raw, knn.p.calibrated, knn.p.calibrated.MCF, file = paste0(output, ".KNN_prediciton.rda"))

####prepare to generate reports
knn.p1 = get_class_score(knn.p.calibrated.MCF)
knn.p2 = get_class_score(knn.p.calibrated)

CNSclassifier$KNN.pred.MCF = names(knn.p1)
CNSclassifier$KNN.pred.MCF = as.vector(knn.p1)
CNSclassifier$KNN.pred.subclass = names(knn.p2)
CNSclassifier$KNN.pred.subclass = as.vector(knn.p2)

CNSclassifier$MCF = sapply( as.character(CNSclassifier$CNS.MCF), function(x) unlist(strsplit(x,";"))[1])
CNSclassifier$MCF1 = sapply(CNSclassifier$MCF, function(x) unlist(strsplit(x,":"))[1])
CNSclassifier$MCF1.score = sapply(CNSclassifier$MCF, function(x) unlist(strsplit(x,":"))[2])

CNSclassifier$Class = sapply(as.character(CNSclassifier$CNS.Subclass), function(x) unlist(strsplit(x,";"))[1])
CNSclassifier$Class1 = sapply(CNSclassifier$Class, function(x) unlist(strsplit(x,":"))[1])
CNSclassifier$Class1.score = sapply(CNSclassifier$Class, function(x) unlist(strsplit(x,":"))[2])

CNSclassifier$RFpurity.ABSOLUTE = as.numeric(as.character(CNSclassifier$RFpurity.ABSOLUTE))
CNSclassifier$RFpurity.ESTIMATE = as.numeric(as.character(CNSclassifier$RFpurity.ESTIMATE))
CNSclassifier$LUMP = as.numeric(as.character(CNSclassifier$LUMP))

setwd("/mnt");
pwd = getwd()

mnpversion = "v11b6"

colnames(targets)[grep("sample.name", colnames(targets),  perl=T,ignore.case=T)] = "Sample_Name";
.. = lapply(1:nrow(targets), function(i){
#    sampleIdat = gsub(paste0("compass/iScan_raw/",centrix,"/"), "", as.character(targets$idat[i]))
    sampleIdat = as.character(targets$idat[i])
    sampleName = as.character(targets$Sample_Name[i])
    sampleDir  = as.character(targets$Outname[i])

    cat(sampleName,sampleIdat,"...")
    #
    cns.p1 = CNSclassifier$MCF1[i]
    cns.p1.score = CNSclassifier$MCF1.score[i]
    cns.p2 = CNSclassifier$Class1[i]
    cns.p2.score = CNSclassifier$Class1.score[i]
    #
    knn.pred.MCFclass = names(knn.p1[i])
    knn.pred.MCFclass.score = round(knn.p1[i],2)
    knn.pred.class = names(knn.p2[i])
    knn.pred.class.score = round(knn.p2[i],2)
    class.table = data.frame(Classifier = c("CNSv11b6"), MCF = c(cns.p1), Score = c(cns.p1.score), 
                             Subclass = c(cns.p2), Score = c(cns.p2.score), check.names = F)
    if(cns.p2 == cns.p1) class.table = data.frame(Classifier = c("CNSv11b6"), 
                                                  Prediction = c(cns.p1), Score = c(cns.p1.score));
    row.names(class.table) = NULL

    histologicDx = paste0("**Diagnosis**: ", targets$Diagnosis[i],". **Age**: ", targets$Age[i], 
                          ", **Gender**:", targets$Gender[i] ,", **Tumor size**: ", 
                          targets$Tumor_site[i],". <br>**Notes**: ", targets$Notes[i],".");
    #prepare data for tsne plot
    plot.col = c(cancer.colour, "#000000");
    names(plot.col)[length(plot.col)] = sampleName;
    x.query = x.test[i, , drop=F]
    x.query.tsne = as.numeric(as.character(test.tsne[i,]))
    tsne = data.frame(rbind(ref.tsne, x.query.tsne));
    tsne$class = c(y.ref, sampleName);
    tsne$size = c(rep(1, length(y.ref)), 4);
    x.query.umap = as.numeric(as.character(test.umap[i,]));
    umap.layout = data.frame(rbind(ref.umap, x.query.umap));
    umap.layout$class = c(y.ref, sampleName);
    umap.layout$size = c(rep(1, length(y.ref)), 4);
    #
    cat("knn.search...")
    res1 = knn.search(x.ref = x.ref, x.query = x.query, y.ref = y.ref, k = 15)
    #
    sex = targets$Gender[i];
    sex = ifelse(sex == "Male" | sex == "M", "Male",ifelse(sex == "Female" | sex == "F", "Female", "unknow"))
    idat_file = targets$Basename[i];
    p.rf.ABSOLUTE = round(CNSclassifier$RFpurity.ABSOLUTE[i], 2);
    p.rf.ESTIMATE = round(CNSclassifier$RFpurity.ESTIMATE[i], 2);
    p.lump = round(CNSclassifier$LUMP[i], 2)
    p.table = data.frame("RFpurity(ABSOLUTE)" = p.rf.ABSOLUTE,
                         "RFpurity(ESTIMATE)" = p.rf.ESTIMATE, "LUMP" = p.lump, check.names=F);
    cat("read array\n");
    RGset <- minfi::read.metharray(idat_file, verbose=FALSE);  Mset <- minfi::preprocessRaw(RGset); rm(RGset);

	cat("generating report\n");
    rmarkdown::render(rmd_file, output_file = paste0( "/mnt/",sampleDir,"/", sampleName,"_MNNv5_",sampleIdat,".html"), 
                      output_dir = paste0("/mnt/",sampleDir), intermediates_dir = paste0("/mnt/",sampleDir), knit_root_dir = paste0("/mnt/",sampleDir),
                      quiet = T);

    nearest.neighbor = NA
    if(sum(res1[[2]]$`N PASS`) > 0 ){
        nearest.neighbor = paste0(as.character(paste0(res1[[2]]$Cancer,"/",
                                  res1[[2]]$`N PASS`)[which(res1[[2]]$`N PASS` > 0)]), collapse = ";");
    }
    nearest.neighbor;
    #break;
})
CNSclassifier$NearestNeighbor = unlist(..)

write.csv(CNSclassifier, file = file.path(idat_path, paste0(centrix, "_KNN.combined.csv")));
for(i in 1:nrow(CNSclassifier)) {
   write.csv(CNSclassifier[i, -c(26,27)], file = paste0(CNSclassifier$Outname[i],"/", CNSclassifier$Sample_Name[i], ".combined.csv"), row.names=F);
   message( CNSclassifier$Sample_Name[i], " combined report saved");
}
message("Centrix: ", centrix , " is done on ", Sys.time()) 
#system(paste0("rm ",output,".cns_tmp.Rda"))
