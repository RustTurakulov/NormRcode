#!/usr/bin/perl
## run this one after you top up reference set samplesheet 
## Y:\NormRcode\Compass_13K.xlsx


print "Script to move and zip idat files from compass folder to MainRefernce pool\n";
$n=0;
$convert = `dos2unix topup_list.txt`;
print $convert;
@barcodes=`cat topup_list.txt`;
foreach $sample (@barcodes) {
	chomp $sample;
	$n++;
	$centrix = $sample;
    $centrix =~ s/_R0\dC0\d//gx;
    $centrix =~ s/^.{_R0\dC0\d//gx;
	@centrix = split "_", $centrix;
    $centrix = pop @centrix;
	$zippedf =  "/data/MDATA/idat/$sample\_Grn.idat.gz";
	if(-e $zippedf){
		print "$sample --> Exists: no action is needed\n";
	}else{
		print "$n\t$sample -->folder: $centrix copying.. \n";
#	    $commad = "cp /data/MDATA/compass/iScan_raw/$centrix/$sample* /data/MDATA/idat";
	    $commad = "cp /data/MDATA/pancancer/iScan_raw/$centrix/$sample* /data/MDATA/idat";
 	    system($commad); 
	}
};
print "\n\nZipping now!\n";
chdir "/data/MDATA/idat";
$pigzipping = "pigz -v *idat";
system($pigzipping); 
print"\nDone transfer of $n sampes!\n";