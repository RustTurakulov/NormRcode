options(scipen = 999, "width"=180)
suppressMessages(library(tidyr))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))

samplesheet   <- "/data/MDATA/NormRcode/CNS13Krefset/Compass_15K.xlsx";

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

compdirs  <- list.files("/data/MDATA/compass/iScan_raw", pattern = "^[0-9]+");
compdirs  <- compdirs[!grepl('[a-zA-Z]', compdirs ,  perl =T)];
newdirs   <- setdiff(compdirs, anno$Slide);

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
   NCI_METRIC = c();
   knnsheet <- paste0("/data/MDATA/compass/iScan_raw/", centrix, "/", centrix,"_KNN.combined.csv");
   if((file.exists(knnsheet))&(file.info(knnsheet)$size)>1){
       kk = kk +1;
         rawsheetfull  <- read.csv(knnsheet, row.names=1);
	     sampletype = rawsheetfull[,grep('Sample_Plate', names(rawsheetfull), ignore.case = T, perl = T)];
	     rawsheet <- rawsheetfull[grepl("Clinical|Brain|CAP", sampletype, ignore.case = TRUE, perl = TRUE), ]; ## Only some rows/types go in
#	     rawsheet <- rawsheetfull[grepl("Clinical|Brain|E3F05|CBTN", sampletype, ignore.case = TRUE, perl = TRUE), ]; ## Only some rows/types go in
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
		 NCI_METRIC           = rawsheet[,which( names(rawsheet) == 'Class1')]; ## Same as CNS.Subclass 
		 Basename             = rawsheet$Basename;

		clms <- list(nnn, Sample,idat_filename,idat,Sentrix_ID,Sentrix_position,material_prediction,Platform_methy,Age,Sex,matched_cases, Location_general, Location_specific, Neoplastic, Primary_category, CNS_study, OS_months, OS_status,PFS_months, PFS_status, Histology,Molecular, Study,predFFPE,CNS.MCF, CNS.MCF.score,CNS.Subclass, CNS.Subclass.score, RFpurity.ABSOLUTE, RFpurity.ESTIMATE, LUMP, NCI_METRIC, Basename)
		names(clms) <- c("nn", 
			"Sample", 
			"idat_filename", 
			"idat", 
			"Sentrix_ID", 
			"Sentrix_position", 
			"material_prediction", 
			"Platform_methy", 
			"Age", 
			"Sex", 
			"matched_cases", 
			"Location_general", 
			"Location_specific", 
			"Neoplastic", 
			"Primary_category", 
			"CNS_study", 
			"OS_months",
			"OS_status", 
			"PFS_months", 
			"PFS_status", 
			"Histology", 
			"Molecular", 
			"Study", 
			"predFFPE",
			"CNS.MCF", 
			"CNS.MCF.score", 
			"CNS.Subclass", 
			"CNS.Subclass.score", 
			"RFpurity.ABSOLUTE", 
			"RFpurity.ESTIMATE", 
			"LUMP", 
			"NCI_METRIC", 
			"Basename");

		clms$CNS.MCF      <- gsub(" |,|/|__", "_", clms$CNS.MCF) 
		clms$CNS.MCF      <- gsub("__", "_",       clms$CNS.MCF) 
		clms$CNS.MCF      <- gsub("__", "_",       clms$CNS.MCF) 
		clms$CNS.Subclass <- gsub(" |,|/|__", "_", clms$CNS.Subclass) 
		clms$CNS.Subclass <- gsub("__", "_",       clms$CNS.Subclass)
		clms$CNS.Subclass <- gsub("__", "_",       clms$CNS.Subclass)
		clms$CNS.Subclass <- gsub("_$", "",        clms$CNS.Subclass, perl=TRUE)
		clms$NCI_METRIC   <- gsub(" |,|/|__", "_", clms$NCI_METRIC) 
		clms$NCI_METRIC   <- gsub("__", "_",       clms$NCI_METRIC)
		clms$NCI_METRIC   <- gsub("__", "_",       clms$NCI_METRIC)
		clms$NCI_METRIC   <- gsub("_$", "",        clms$NCI_METRIC, perl=TRUE)

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

message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\t", 
 nx, "\t-- chips processed\n\t", 
 kk, "\t-- with CNS sample sheets\n\t", 
 nrow(newsamplesheet),"\t-- new samples collected\n",
"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

write.csv(newsamplesheet, "/data/MDATA/NormRcode/CNS13Krefset/newsamples.csv", row.names=F);
