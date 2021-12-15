#!/usr/bin/perl 
use strict;
use warnings;
print("Script for sbatching MNN  pipeline");
print("Launching multiple instances of mnn_rbox.sh");

## Wiping slurm logs
chdir "/data/MDATA/compass/slurm_logs";
opendir(DIRHANDLE, ".") or die "Cannot open logdir: $!";
my @logfiles = grep { -e and /^slurm-\d+.out/ } readdir DIRHANDLE;
unlink for @logfiles;

 
my $dir = "/data/MDATA/compass/iScan_raw";
chdir $dir;
opendir DIR, $dir or die "Could not open '$dir' for reading '$!'\n";
my @centrix = grep(/^\d+$/gx, readdir DIR);
closedir DIR;

my $k=0;
foreach my $folder (@centrix) {
	 chdir $dir;
	 my $results = "$folder/$folder\_KNN.combined.csv";
	 print("\n $folder >> ");
	 if (-e $results) {
		 print "$folder\_KNN.combined.csv -- KNN report exists";
	 }else{
		 open INFILE, "<$folder/Sample_Sheet.csv" || die "No samplesheet in $folder/\n$!\n"; 
		 my @samplesheet = <INFILE> ; 
		 close INFILE ;
         my @brain = grep { /,Clinical|,Sarcoma Project,|,CBTN_|,RES_E3F05/ } @samplesheet;
#         my @brain = grep { /,Clinical|,Sarcoma Project,/ } @samplesheet;
		 if ((scalar @brain)>1) {
			 $k++;
			 chdir "/data/MDATA/compass/slurm_logs";
			 my $command = "sbatch /data/MDATA/NormRcode/mnn_rbox.sh $folder";
			 print "Submit job $k to squeue:\n   $command\n ";
			 system($command);
		 }else{
			 print " No [Clinical: Brain] category in $folder samplesheet";
		 };
     };
};
print "\n$k jobs were sent to sbatch\n";