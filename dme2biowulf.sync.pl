#!/usr/bin/perl
### Script for idat files  pull from DME to Biowulf. 4/23/2021
### Just new batches which are not in Excel master file
use lib "/data/MDATA/NormRcode/perl_modules/lib/perl5/";
use warnings;
use strict;
use Spreadsheet::Read qw(ReadData);
use Excel::Writer::XLSX;

### Setting environment and proxies details for DME API. There are other files with constant parameters 
### to tune straight after API install from GIT [ git clone https://github.com/CBIIT/HPC_DME_APIs.git ]
### Major config file is here:
### /home/turakulovr2/HPCAPI/HPC_DME_APIs/utils/hpcdme.properties
##       #Production server settings
##		hpc.server.url=https://hpcdmeapi.nci.nih.gov:8080
##		hpc.ssl.keystore.path=hpc-client/keystore/keystore-prod.jks
##		hpc.ssl.keystore.password=changeit
##		#Proxy Settings
##		hpc.server.proxy.url=10.1.200.240
##		hpc.server.proxy.port=3128
##		#Globus settings
##		#default globus endpoint to be used in registration and download
##		hpc.default.globus.endpoint=7fd394ac-a128-11eb-8a91-d70d98a40c8d
##		#HPC DME Login token file location
##		hpc.login.token=hpcdme-auth.txt
## More documentation is on wiki page:
##      https://wiki.nci.nih.gov/display/DMEdoc
## Install perl modules locally [note the path line N4 at the top: (use lib)]
##      cpanm --local-lib=/data/MDATA/NormRcode/perl_modules Excel::Writer::XLSX 
##      cpanm --local-lib=/data/MDATA/NormRcode/perl_modules Spreadsheet::Read qw(ReadData)

$ENV{HPC_HOME}     = "/home/turakulovr2/HPCAPI/HPC_DME_APIs";
$ENV{HPC_DM_UTILS} = "/home/turakulovr2/HPCAPI/HPC_DME_APIs/utils";
$ENV{PATH}        .= ":/home/turakulovr2/HPCAPI/HPC_DME_APIs/utils/scripts";
$ENV{http_proxy}   = "http://dtn01-e0:3128";
$ENV{ftp_proxy}    = "http://dtn01-e0:3128";
$ENV{RSYNC_PROXY}  = "dtn01-e0:3128";
$ENV{https_proxy}  = "http://dtn01-e0:3128";


## Setting output and config files
my $dmepath  = "LP_COMPASS_Archive/PI_Kenneth_Aldape/Project_Methylation/Sentrix_"; #6.5K files 
#my $dmepath = "LP_COMPASS_Archive/PI_Kenneth_Aldape/Project_Methylation/Sentrix_202827620074";
my $fnmlist = "TRANSFER/dme_methylation_list.txt";
my $jsndump = "TRANSFER/dme_methylation_list.json";
my $jsncnfg = "NormRcode/methylation.ss.json";        # pull only files from 2021 
my $destination  = "TRANSFER/compass" ;               # ad as vaiable later
my $bigexcelfile = "$destination/Sample_Sheet.xlsx";  # concatenated samplesheets

### Stage 1.  Collect and parse metadata from DME with API
### Parse first page of files (there is limit of 100 files) to save.   
print "\nLaunching metadata dump from DME:\n\n";
print("dm_query_dataobject -D $fnmlist -o $jsndump $jsncnfg $dmepath\n\n");
system("dm_query_dataobject -D $fnmlist -o $jsndump $jsncnfg $dmepath");

my @ALLFILES;
my @JASON = `cat NormRcode/methylation.json`;    
my @PAGE  = `cat TRANSFER/dme_methylation_list.json`;
my $limit = 0;
my $totalCount =0;
foreach my $line (@PAGE) {
	chomp $line;
    if ($line =~ m/limit/) {
       $limit = $line;
       $limit =~ s/\D//gx;
    };
    if ($line =~ m/totalCount/) {
       $totalCount = $line;
       $totalCount =~ s/\D//gx;
	};
	if ($line =~ m/idat|xlsx/) {
	    $line =~ s/\s|\"|,//gx;
		push @ALLFILES, $line;
	};
};
print "LIMIT = $limit\n";
print "TOTAL = $totalCount\n\n";

# Now repeating query for all pages by updating page number in JASON file
my $lastpage = int($totalCount / $limit); 
my $curentpage = 1;
for ($curentpage .. $lastpage) {
	$curentpage++;
	open JASON, ">NormRcode/methylation.json";
	foreach my $jline (@JASON) {
		if ($jline =~ m/ \"page\":/g) {
            print JASON "   \"page\": $curentpage,\n";
		}else{
		    print JASON $jline;
		}
	};
	close JASON;
    print "Collecting page $curentpage out of $lastpage\n";
    system("dm_query_dataobject -D TRANSFER/dme_methylation_list.txt -o TRANSFER/dme_methylation_list.json NormRcode/methylation.json $dmepath");
    @PAGE  = `cat TRANSFER/dme_methylation_list.json`;
    foreach my $line (@PAGE) {
		if ($line =~ m/idat|xlsx/) {
		    $line =~ s/\s|\"|,//gx;
		    push @ALLFILES, $line;
	    };
    };
};
#Putting page counter back to page 1
open JASON, ">NormRcode/methylation.ss.json";
foreach my $jline (@JASON) {
	if ($jline =~ m/ \"page\":/g) {
		print JASON "   \"page\": $curentpage,\n";
	}else{
		print JASON $jline;
	}
};
close JASON;
print "\nTotal files collected: ".scalar(@ALLFILES)."\n";
print "\nStart sifting with master sheet\n";

## Dump full list to output file and create hash
my %DMEIDAT =();
my %DMESSHT =();
open FULLIST, ">>$fnmlist" ;
foreach my $line (@ALLFILES) {
   print FULLIST "$line\n";
   my @LINE = split "/", $line;
   my $file = pop @LINE;
   if ($file =~ m/idat$/) {
         $file =~ s/gz$//gx;
	     $file =~ s/_Grn.idat$|_Red.idat$//gx;
         $DMEIDAT{$file} .= " $line"; 
   };
   if ($file =~ m/xlsx$/) {
	   my $cxbarcode = $line;
	   $cxbarcode =~ s/\/Raw_Data\/Sample_Sheet_batch.xlsx$//g;
	   $cxbarcode =~ s/\///gx;
	   $cxbarcode =~ s/^\w+Sentrix_//gx;
       $DMESSHT{$cxbarcode}=$line;
   };
};
close FULLIST;

### Stage 2. Collect hash of centrix existed in mastersheet 
my %bookcentrix = ();
my $book = ReadData('NormRcode/Sample_sheet.xlsx');
my @rows = Spreadsheet::Read::rows($book->[1]);
foreach my $i (1 .. $#rows) {
	my $dnaid   = $rows[$i][1];
   my $centrix = $rows[$i][2];
	$bookcentrix{$centrix} .= "$dnaid ";
};

foreach my $centrix (sort keys %DMEIDAT){
	if (exists $bookcentrix{$centrix}) {
        print  "\n\n $centrix -- $bookcentrix{$centrix} is already in master sheet\n";
	}else{
        print  "\n\n $centrix -- DNA is NOT in master sample sheet\n";
		my @FILES = split " ", $DMEIDAT{$centrix};
		my @unique = uniq( @FILES ); 
		foreach my $idat (@unique) {
		   my @PATH = split "/", $idat;
		   my $fname = pop @PATH;
           my $idattarg = $destination ."/". $fname; 
           if (-e $idattarg) {
			   print "\n\t$fname -- transfered before";
		   }else{
			   my $transfer = "dm_download_dataobject $idat $idattarg";
  			   print "\n\t$fname -- transferring now\n";
			   system ($transfer) ;
          };
		};
   
        ## Checking  samplesheet for each centrix sample 
		my $cxbarcode = $centrix;
		   $cxbarcode =~ s/_R\d+C\d+$//gx;
		if (exists $DMESSHT{$cxbarcode}) {
			print "\nHave samplesheet file for this centrix barcode on DME\n";
			my $sshtname = $DMESSHT{$cxbarcode};
			my $sshttarg = $destination ."/". $cxbarcode . "_Sample_Sheet.xlsx"; 
			   if (-e $sshttarg) {
				   print "\t$cxbarcode -- samplesheet was transfered before\n";
			   }else{
				   my $transfer = "dm_download_dataobject $sshtname $sshttarg";
				   print "\n\t".$cxbarcode . "_Sample_Sheet.xlsx -- transferring now\n";
				   system ($transfer) ;
			  };
		}else{
			print "\nProblem: have NO samplesheet file for this centrix barcode on DME\n";
		}
	}
};
print "\n\n  DONE\n\n";


### Stage 3. Concatenate excel files (merge samplesheets)
my @excelfiles = ();
my @data = ("X");
opendir(DIR, $destination) or die "couldn't open $destination: $!\n";
@excelfiles = grep { $_ =~ /\d+_Sample_Sheet.xlsx/ } readdir DIR;
closedir DIR;

## Stuff array
for my $xlf (@excelfiles)
{
 my $excelfile = "$destination/$xlf"; 
 print "\n$excelfile";  
 my $book = ReadData($excelfile);
 my @extradata = Spreadsheet::Read::rows($book->[1]) ;
	foreach my $i (2 .. scalar @extradata) {                     ## 2 <- No header is needed
		my @row=();
		foreach my $j (1 .. 9) {   ## scalar @{$extradata[$i-1]} ## if full line is needed: there is some useless stuff after column 10  
			push @row, $extradata[$i-1][$j-1] ;
		}
#		my $line = join "\t", @row;  this generates warning messages.
		my $line = join( "\t", map { defined ? $_ : '' } @row );
        push @data, $line;
	}
};

## Dump concatenated array to Excel file
my $workbook = Excel::Writer::XLSX->new("$bigexcelfile");
my $worksheet = $workbook->add_worksheet();
my $firstline = 1;
my ($x,$y) = (1,0);
foreach my $line (@data)
{
  chomp $line;
  if ($firstline eq 1) # Header Lines replace "X"
  {
    $firstline++;
    $worksheet->write( 0, 0, "Sample name" );
    $worksheet->write( 0, 1, "Sentrix id" );
    $worksheet->write( 0, 2, "Array" );
    $worksheet->write( 0, 3, "Preparation type");
    $worksheet->write( 0, 4, "Gender");
    $worksheet->write( 0, 5, "Diagnosis");
    $worksheet->write( 0, 6, "Location");
    $worksheet->write( 0, 7, "Age (in years)");
    $worksheet->write( 0, 8, "Notes");
  }
  else
  {    $firstline++;
    my @formatline = split( '\t' , $line );
    foreach my $cell (@formatline)
    {
      $worksheet->write($x, $y++, $cell);
    }
    $x++;$y=0;
  }
} 
$workbook->close;

print "\n~~~~~~~~~~~~~~~\n   The  end.  \n~~~~~~~~~~~~~~~\n"; 


## Not very useful piece
sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}