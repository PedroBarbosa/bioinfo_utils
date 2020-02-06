use strict;
use warnings;

############
#Main start#
############
# It executes eqclean from an sff file.
# You must pass the sff file and the type of data 0 for transcriptomic data; 1 for genomic data.
 

my $single=$ARGV[0];
my $dataType=$ARGV[1]; # 0 - transcriptomic; 1 - genomic
my $logFile=$ARGV[2]; # logFile
my $index=rindex ($single, "/");
my $genomic_directory="/mnt/msa/BIOCANT/genomic-data/SFF_genom/";
my $transcriptomic_directory="/mnt/msa/BIOCANT/transcriptomic-data/SFF_trans/";
my $nameFile=substr $single, $index+1;
my $outputFasta="";
my $outputSeqClean="";
my $outputCln="";
my $outputNewIDS="";
my $outputIDS="";
my $outputTrimPoints="";
my $outputTRIM="";
my $out="";
do {
  local *STDOUT;
  local *STDERR;
  # redirect STDOUT to log.txt
  open (STDOUT, '>>', $logFile);
  # redirect STDOUT to log.txt
  open (STDERR, '>>', $logFile);

if ($dataType == 0 ){
  $outputFasta=$transcriptomic_directory."FASTA_FILES/".$nameFile.".fna";
  $outputSeqClean=$transcriptomic_directory."SEQCLEAN/".$nameFile.".cleaned";
  $outputCln=$nameFile.".fna.cln"; # current directory
  $outputNewIDS=$transcriptomic_directory."SEQCLEAN/TRIM/".$nameFile.".newIDS";
  $outputIDS=$transcriptomic_directory."SEQCLEAN/TRIM/".$nameFile.".IDS";
  $outputTrimPoints=$transcriptomic_directory."SEQCLEAN/TRIM/".$nameFile.".newtrimpoints";
  $outputTRIM=$transcriptomic_directory."SEQCLEAN/TRIM/".$nameFile.".trim";
  }
elsif ($dataType == 1){
  $outputFasta=$genomic_directory."FASTA_FILES/".$nameFile.".fna";
  $outputSeqClean=$genomic_directory."SEQCLEAN/".$nameFile.".cleaned";
  $outputCln=$nameFile.".cln"; # current directory
  $outputNewIDS=$genomic_directory."SEQCLEAN/TRIM/".$nameFile.".newIDS";
  $outputIDS=$genomic_directory."SEQCLEAN/TRIM/".$nameFile.".IDS";
  $outputTrimPoints=$genomic_directory."SEQCLEAN/TRIM/".$nameFile.".newtrimpoints";
  $outputTRIM=$genomic_directory."SEQCLEAN/TRIM/".$nameFile.".trim";
  }
else{
  print "Error in dataType\n";
  exit(-1);
  }

# Step 1: Convert SFF to FASTA
my $sff2fasta="sffinfo -seq -notrim $single > $outputFasta";
print "COMMAND: $sff2fasta \n";
#system($sff2fasta);
print qx/$sff2fasta/;
#print "$out\n";

# Step 2: Execute SeqClean

my $seqclean="seqclean $outputFasta -n 10000 -o $outputSeqClean";
print "COMMAND: ".$seqclean."\n";
#system($seqclean);
print qx/$seqclean/;
# print "$out \n";

# TRIMMING POINTS

# Step 3: Get the IDS of valids and trimmed reads
my $getIDS="cat $outputCln | grep -vE \"(short|low_qual|dust)\" | cut -f 1 > $outputNewIDS"; 
print "COMMAND: $getIDS \n";
#system($getIDS);
$out=qx/$getIDS/;
print "$out \n";
#Step 4: Get the reads in the SFF with the IDS of valids and trimmed reads
my $getReadsIDS="sfffile -o $outputIDS -i $outputNewIDS $single";
print "COMMAND: $getReadsIDS \n";
#system($getReadsIDS);
$out=qx/$getReadsIDS/;
print "$out \n";
# Step 5: Get the coordinates for the new trimpoints
my $getTrimPoints="cat $outputCln | grep -vE \"(short|low_qual|dust)\" | cut -f 1,3,4 > $outputTrimPoints";
print "COMMAND: $getTrimPoints \n";
#system($getTrimPoints);
$out=qx/$getTrimPoints/;
print "$out\n";
#Step 6: Adjust the trim points for the reads in the .sff file
my $getReadsTRIM="sfffile -o $outputTRIM -t $outputTrimPoints $outputIDS";
print "COMMNAD: $getReadsTRIM \n";
#system($getReadsTRIM);
$out=qx/$getReadsTRIM/;
print "$out \n";
};
print "DONE\n";
exit;
