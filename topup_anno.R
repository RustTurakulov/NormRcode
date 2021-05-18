#!/usr/bin/env Rscript
# 5/1/2021 By Rust
library(openxlsx)
sourcedir = "compass";

#args = commandArgs(trailingOnly=TRUE)
#sourcedir <- args[1] ; 
#
#if (length(args)==0) {
#  stop("Input directory must be supplied", call.=FALSE)
#} else if (length(args)==1) {
#  message("Folder for data topup: ", sourcedir);
#}

# identify files anfolders
setwd("/data/MDATA/TRANSFER");
newfolder  <- paste0("/data/MDATA/TRANSFER/", sourcedir );
samplesout <- paste0("/data/MDATA/TRANSFER/SAMPLESHEETS/", sourcedir, "_", Sys.Date(), ".xlsx" );
bigpooldir <- "/data/MDATA/idat";

# Use concatenated Sample_Sheet.xlsx to find files with new data, parsing data only with filled Sample.name
#samplesheet   <- list.files(newfolder, "Sample_Sheet.xlsx$"); ## This will be list of files
samplesheet   <- "Sample_Sheet.xlsx";
newidatlist   <- list.files(newfolder, pattern = ".idat$|.idat.gz$");
anno_base     <- openxlsx::read.xlsx("/data/MDATA/NormRcode/Sample_sheet.xlsx");
# weeding duplicates in mastersheet (this mess need to be cleaned one day)
anno <- anno_base[!duplicated(anno_base$idat_filename), ]
centrixpool = unique(anno$idat_filename) ;
nn <- length(centrixpool);
#anno$order <- c(1:nn);

newbatch_anno <- openxlsx::read.xlsx(paste0(newfolder, "/", samplesheet));
newbatch_anno <- newbatch_anno[!is.na(newbatch_anno[,1]),1:9];            ## DNA id is must
newbatch_anno <- newbatch_anno[!duplicated(newbatch_anno$Sentrix.id), ];  ## Duplicates are pain
## Leaving only new samples for addition: left join on centrix IDs in new samplesheet
newbatch_anno <- newbatch_anno[newbatch_anno$Sentrix.id %in% setdiff(newbatch_anno$Sentrix.id, centrixpool),];
rownames(newbatch_anno) <- newbatch_anno$Sentrix.id;

message("~~~~ Raw data transfer ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# copy idat files to the major pool folder, 
# those have to be compressed for later normalisation script intake

for(sample in newbatch_anno[,1]) {
	  centrix = newbatch_anno[newbatch_anno[,1]==sample, 2];
	  rgfiles = newidatlist[grep(centrix, newidatlist)]; 
	  for(idatfile in rgfiles) {
          if(file.exists(paste0(bigpooldir,"/",idatfile, ".gz"))) {
			  message(idatfile, ": is already in base pool")
		  }else{
			  if (grepl(".gz$", idatfile)) {
				  message(idatfile, ": transfer as it is"); 
				  file.copy(paste0(newfolder,"/",idatfile), bigpooldir);
			  }else{
				  message(idatfile, ": GZip and transfer");  
				  system(paste0("pigz -k ",newfolder,"/",idatfile," && mv -v ", newfolder,"/",idatfile, ".gz ", bigpooldir))
			  } 
	     }
	  }
}

message("~~~~ Adding new files metadata to the master spreadsheet ~~~~");
# Appending master spreadsheet of centrixpool with new batch samples

for(centrix in newbatch_anno[,2]) {
	 sample = newbatch_anno[newbatch_anno[,2]==centrix, 1] ;
#	 options(warn=-1);
     if (centrix %in% centrixpool){
		 message(centrix, ":        has record in master file");
         # add later: append missed fields if any from newbatch_anno --> anno_base 
	 }else{
		 message(centrix, ":        adding to pool annotation");
		 nn <- nn + 1 ;
	     today       = format(Sys.Date(), "%d/%m/%Y");
		 arraytype   = newbatch_anno[centrix, 3] ;
		 if(grep("EPIC", arraytype)){
             arraytype = "HumanMethylationEPIC";
         }else{
			 arraytype = "HumanMethylation450";
		 };
	     preparation = newbatch_anno[centrix, 4] ;
		 if(is.na(preparation)){
               preparation = "FFPE"
		 }; # No empty fields here otherwise problem with normalization script 
		 gender      = newbatch_anno[centrix, 5] ;
         diagnosis   = newbatch_anno[centrix, 6] ;
		 location    = newbatch_anno[centrix, 7] ;
		 age         = newbatch_anno[centrix, 8] ;
         notes       = newbatch_anno[centrix, 9] ;
         barcode     = substr(centrix,  1, 12);
		 position    = substr(centrix, 14, 20);
#		 idat_fn_alt =	paste0(newfolder, centrix);

                                    #Columns names in New excel file
         anno[nn, 1] <- nn; 	     	#		nn	1
         anno[nn, 2] <- sample;		    #		Sample	DKFZ_INF_335
         anno[nn, 3] <- centrix;    	#		idat_filename	10005771111_R06C02
         anno[nn, 4] <- centrix;		#		idat	10005771111_R06C02
         anno[nn, 5] <- barcode;  		#		Sentrix_ID	10005771111
         anno[nn, 6] <- position;		#		Sentrix_position	R06C02
         anno[nn, 7] <- preparation;	#		material_prediction	Frozen
         anno[nn, 8] <- arraytype;		#		Platform_methy	HumanMethylation450
         anno[nn, 9] <- age;     		#		Age	0
         anno[nn,10] <- gender;  		#		Sex_prediction	M
         anno[nn,11] <- "";     		#		matched_cases	
         anno[nn,12] <- location;		#		Location_general	
         anno[nn,13] <- "Neoplastic";	#		Location_specific	Neoplastic
         anno[nn,14] <- "Neuropathology";#		Neoplastic	Neuropathology
         anno[nn,15] <- "Case"; 		#		Primary_category	Case
         anno[nn,16] <- "Case";  		#		CNS_study       	Case
         anno[nn,17] <- "";		        #		OS_months      #
         anno[nn,18] <- ""; 	    	#		OS_status	1  #  survival 
         anno[nn,19] <- "";	        	#		PFS_months	   #   data
         anno[nn,20] <- "";    		    #		PFS_status	   #
		 anno[nn,21] <- diagnosis;		#		Histology	
         anno[nn,22] <- "";          	#		Molecular	
         anno[nn,23] <- "NCI_LP";    	#		compass folder	
		 anno[nn,24] <- " ";    	   	#		predFFPE	Frozen              #             
		 anno[nn,25] <- " ";    		#		CNS.MCF	IHG                     # 
		 anno[nn,26] <- " ";    		#		CNS.MCF.score	0.953           #
		 anno[nn,27] <- " ";    		#		CNS.Subclass	IHG             #   parse from classifier 
		 anno[nn,28] <- " ";    		#		CNS.Subclass.score	0.953       #   to do later 
		 anno[nn,29] <- " ";    		#		RFpurity.ABSOLUTE	0.6349255   #
		 anno[nn,30] <- " ";    		#		RFpurity.ESTIMATE	0.923143717 #
		 anno[nn,31] <- " ";    		#		LUMP	0.950736459             # 


#                                               #Columns names in Drew's excel file
#         anno[nn, 1]  <- nn;                       #order
#         anno[nn, 2]  <- sample;                   #DNA id
#         anno[nn, 3]  <- centrix;                  #Chip barcode full
#         anno[nn, 4]  <- idat_basename;            #Basename 
#         anno[nn, 5]  <- idat_basename;            #Basename AWS
#         anno[nn, 6]  <- idat_basename;            #basename Ubuntu
#         anno[nn, 7]  <- idat_fn_alt;              #idat_filename_alt
#         anno[nn,27]  <- "NIH";                    #Primary study
#         anno[nn,28]  <- centrix;                  #idat
#         anno[nn,29]  <- barcode;                  #Sentrix_ID
#         anno[nn,30]  <- position;                 #Sentrix_position		 
#         anno[nn,31]  <- "FALSE";                  #duplicate_idat
#         anno[nn,32]  <- preparation;              #material_prediction
#         anno[nn,33]  <- paste0("HumanMethylation",arraytype);   #Platform_methy
#         anno[nn,34]  <- "NIH.NCI.LP";             #Center_methy
#         anno[nn,35]  <- sourcedir;                #Accession_methy
#         anno[nn,36]  <- centrix;                  #Filename_methy
#         anno[nn,57]  <- age;                      #Age
#         anno[nn,58]  <- gender;                   #Sex
#         anno[nn,68]  <- today;                    #idat_or_published_date
#         anno[nn,108] <- diagnosis;                #Likely_integrated_diagnosis
#         anno[nn,110] <- "Case";                   #pan_study
#         anno[nn,111] <- "Case";                   #panCNS_study
#         anno[nn,112] <- "case";                   #Copy_number_CNS
	 }
#     options(warn=0);
};

## There are some toxic chips hanging around and makes life miserable
blacklist = c("203723060097_R03C01",
			  "203219730196_R08C01",
			  "204957740079_R08C01",
			  "202831040010_R06C01",
			  "203723050099_R06C01");

message("~~~~ Are idat files from merged samplesheet are in main pool now?");
updatedpool <- list.files(bigpooldir, pattern = ".idat.gz$");
missed=0;
for(centrix in anno[,4]){
	kk <- length(grep(centrix, updatedpool));
	if(kk < 2){
        message( centrix, " -- has files missed in idat/ ");
		missed = missed + 1;
		blacklist  <- append(blacklist, centrix);
	}
};
anno <- anno[!(anno$idat_filename %in% blacklist), ];

write.xlsx(anno, samplesout,  colNames = TRUE );
newnn <- abs(nn - length(centrixpool));
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
message(paste0("Idat files transfered and new samplesheet [+", newnn ,", -", missed, "] saved to:"));
message(samplesout); 
message(Sys.time());