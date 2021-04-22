#!/usr/bin/env Rscript
library(openxlsx)
args = commandArgs(trailingOnly=TRUE)
sourcedir <- args[1] ; 

if (length(args)==0) {
  stop("Input directory must be supplied", call.=FALSE)
} else if (length(args)==1) {
  print(paste("Folder for data topup: ", sourcedir), quote=F);
}

# identify files anfolders
setwd("/data/MDATA/TRANSFER");
newfolder  <- paste0("/data/MDATA/TRANSFER/", sourcedir );
samplesout <- paste0("/data/MDATA/TRANSFER/SAMPLESHEETS/", sourcedir, "_", Sys.Date(), ".xlsx" );
bigpooldir <- "/data/MDATA/idat";

# find specific samplesheet files with new data, parsing data only with filled Sample.name
samplesheet   <- list.files(newfolder, "Sample_Sheet_batch.xlsx$")
newidatlist   <- list.files(newfolder, pattern = ".idat$|.idat.gz$")
anno_base     <- openxlsx::read.xlsx("/data/MDATA/NormRcode/Sample_sheet_test.xlsx")
newbatch_anno <- openxlsx::read.xlsx(paste0(newfolder, "/", samplesheet))
newbatch_anno <- newbatch_anno[!is.na(newbatch_anno[,1]),1:9]; 


print("~~~~ Raw data transfer ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", quote=F)
# copy idat files to the major pool folder, 
# those have to be compressed for later normalisation script intake

for(sample in newbatch_anno[,1]) {
	  centrix = newbatch_anno[newbatch_anno[,1]==sample, 2];
	  rgfiles = newidatlist[grep(centrix, newidatlist)]; 
	  for(idatfile in rgfiles) {
          if(file.exists(paste0(bigpooldir,"/",idatfile, ".gz"))) {
			  print(paste(idatfile, ": is already in base pool"), quote=F)
		  }else{
			  if (grepl(".gz$", idatfile)) {
				  print(paste(idatfile, ": transfer as it is"), quote=F) 
				  file.copy(paste0("cp -nv ",newfolder,"/",idatfile," ", bigpooldir));
			  }else{
				  print(paste(idatfile, ": GZip and transfer"), quote=F);  
				  system(paste0("gzip -k ",newfolder,"/",idatfile," && mv -v ", newfolder,"/",idatfile, ".gz ", bigpooldir))
			  } 
	     }
	  }
}

print("~~~~ Adding new files metadata to the master spreadsheet ~~~~~~~~~", quote=F);
# Appending master spreadsheet of centrixpool with new batch samples
centrixpool = anno_base$idat_filename ;
nn <- length(centrixpool);

for(sample in newbatch_anno[,1]) {
	 centrix  = newbatch_anno[newbatch_anno[,1]==sample, 2] ;
     if (grepl(centrix, centrixpool)){
		 print(paste(centrix, " has record in master file"), quote=F);
         # add later: append missed fields if any from newbatch_anno --> anno_base 
	 }else{
		 print(paste(centrix, " adding to pool annotation"), quote=F);
		 nn <- nn + 1 ;
	     arraytype   = newbatch_anno[centrix, 3] ;
	     preparation = newbatch_anno[centrix, 4] ;
		 gender      = newbatch_anno[centrix, 5] ;
         diagnosis   = newbatch_anno[centrix, 6] ;
		 location    = newbatch_anno[centrix, 7] ;
		 age         = newbatch_anno[centrix, 8] ;
         notes       = newbatch_anno[centrix, 9] ;
         barcode     = substr(centrix,  1, 12);
		 position    = substr(centrix, 14, 20);
		 idat_basename  = paste0("idat/",centrix);

                                                       #Columns names in Drew's excel file
         anno_base[nn, 1]  <- nn;                       #order
         anno_base[nn, 2]  <- sample;                   #DNA id
         anno_base[nn, 3]  <- centrix;                  #Chip barcode full
         anno_base[nn, 4]  <- idat_basename;            #Basename 
         anno_base[nn, 5]  <- idat_basename;            #Basename AWS
         anno_base[nn, 6]  <- idat_basename;            #basename Ubuntu
         anno_base[nn,27]  <- "NIH";                    #Primary study
         anno_base[nn,28]  <- centrix;                  #idat
         anno_base[nn,29]  <- barcode;                  #Sentrix_ID
         anno_base[nn,30]  <- position;                 #Sentrix_position		 
         anno_base[nn,31]  <- "FALSE";                  #duplicate_idat
         anno_base[nn,32]  <- preparation;              #material_prediction
         anno_base[nn,33]  <- paste0("HumanMethylation",arraytype);   #Platform_methy
         anno_base[nn,34]  <- "NIH";                    #Center_methy
         anno_base[nn,35]  <- newfolder;                  #Accession_methy
         anno_base[nn,36]  <- centrix;                  #Filename_methy
         anno_base[nn,57]  <- age;                      #Age
         anno_base[nn,58]  <- gender;                   #Sex
         anno_base[nn,68]  <- Sys.Date();               #idat_or_published_date
         anno_base[nn,108] <- diagnosis;                #Likely_integrated_diagnosis
         anno_base[nn,110] <- "Case";                   #pan_study
         anno_base[nn,111] <- "Case";                   #panCNS_study
         anno_base[nn,112] <- "case";                   #Copy_number_CNS
	 }
};
write.xlsx(anno_base, samplesout,  colNames = TRUE )
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", quote=F)
print("Idat files transfered and samplesheet saved to:", quote=F); 
print(samplesout, quote=F); 
print(Sys.time(), quote=F);