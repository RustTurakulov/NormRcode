# NormRcode

Code for Normalization Of Methylation Array with Meffil

---

#### Sample_sheet_test.xlsx 
 *. Example samplesheet. That is master database which is collated for the major reference dataset.

#### dme2biowulf.sync.pl 
 *.  The first script to run in pipeline. It pull data from DME folder which has extension .idat or Samle_Sheet_Batch.xlsx to the TRANSFER/ folder. 
Excel file will be renamed based on centrix barcode (not to overwrite).
Excel files with samplesheet will be concatenated to single excel file for the pipeline. 

#### methylation.json
*. Configuration file for DME metadata pull for dme2biowulf.sync.pl script. Just to search \[ \*.idat\ ] and \[ SampleSheet_batch.xlsx \]


#### topup_anno.R
*Script for zipping and transfering copy of \*.idat files from TRANSFER folder to major pool folder. This script is also updating master database (thi big excel file) with metadata this excel is used in normalization script late.

#### norm.sh
*Wrapper for* **normalization.R**. * This code will go over Drew's pipeline and generates Umap files for R-Shiny ingestion. 

----

### Script order

*  dme2biowulf.sync.pl (No special requarements on machine configuration)
*  topup_anno.R   (R environment, or singularity with rbox_v0.sif) 
*  normalization.R  ( singularity with rbox_v0.sif otherwise install bunch of R packages)


