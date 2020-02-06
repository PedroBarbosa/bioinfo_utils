use strict;
use warnings;

############
#Main start#
############

# Creating the names of the output files from the input file (the first PE of each pair files), using regular expressions to replace the endings/extensions.
# It will create the corresponding WithN and NoN files for each pair end, together with an output file for the record of some statistics. 
# It also will apply the trimming of the NoN-reads executing sickle.
# It will save the results  in the defined directory.

## Assigning parameters for WithN and NoN reads
my $pe_1=$ARGV[0];
my $directory=$ARGV[1];

## Assigning parameters for Sickle
my $qual=$ARGV[2];
my $len=$ARGV[3];

#Create the directory
unless(-e $directory or mkdir($directory, 0775)) {
  die "Unable to create $directory ". $!."\n";
  }
## Defining name file for PE-2
my $findPE=quotemeta "_1.fq";
my $replacePE="_2.fq";
my $pe_2=$pe_1;
$pe_2=~ s/$findPE/$replacePE/g; # replace the name file

## Defining output files for WithN and NoN reads
my $findN= quotemeta ".fq";
my $replaceWithN=".WithN.fq";
my $replaceNoN=".NoN.fq";
my $replace_out=".output_trimmingN.txt";

my $WithNpe_1=$pe_1;
my $NoNpe_1=$pe_1;
my $WithNpe_2=$pe_2;
my $NoNpe_2=$pe_2;
my $outputWithN_NoN = $pe_1;

$WithNpe_1=~ s/$findN/$replaceWithN/g; # replace the name file
$NoNpe_1=~ s/$findN/$replaceNoN/g; # replace the name file
$WithNpe_2=~ s/$findN/$replaceWithN/g; # replace the name file
$NoNpe_2=~ s/$findN/$replaceNoN/g; # replace the name file
$outputWithN_NoN=~ s/$findPE/$replace_out/g;

#print "$pe_1 $pe_2 $WithNpe_1 $WithNpe_2 $NoNpe_1 $NoNpe_2\n";
#exit;
# Removing the PATH of PE files and defining the output directory.
my $index=rindex ($pe_1, "/");
$WithNpe_1=$directory."/".substr $WithNpe_1, $index+1;
$NoNpe_1=$directory."/".substr $NoNpe_1, $index+1;
$WithNpe_2=$directory."/".substr $WithNpe_2, $index+1;
$NoNpe_2=$directory."/".substr $NoNpe_2, $index+1;
$outputWithN_NoN=$directory."/".substr $outputWithN_NoN, $index+1;

## Defining output files for sickle
my $find= quotemeta ".fq"; 
my $findNoNPE= quotemeta "_1.NoN.fq"; 
my $replacePE_Sickle=".Q".$qual."L".$len.".fq";
my $replace_Single=".NoN.Single".$replacePE_Sickle;
my $replace_Sickle=".NoN.output_Q".$qual."L".$len.".txt";

my $sickle_pe_1=$NoNpe_1;
my $sickle_pe_2=$NoNpe_2;
my $single_out=$NoNpe_1; 
my $output_Sickle=$NoNpe_1;

$sickle_pe_1=~ s/$find/$replacePE_Sickle/g;
$sickle_pe_2=~ s/$find/$replacePE_Sickle/g;
$single_out=~ s/$findNoNPE/$replace_Single/g;
$output_Sickle=~ s/$findNoNPE/$replace_Sickle/g;

##### WithN and NoN reads ####
my $totalReadN_p1=0;
my $totalReadN_p2=0;
my $totalN_p1=0;
my $totalN_p2=0;

open (F1, ">".$WithNpe_1) or die "cannot open for writting! " . $!;
open (F2, ">".$NoNpe_1) or die "cannot open for writtinq! " . $!;
open (F3, ">".$WithNpe_2) or die "cannot open for writting! " . $!;
open (F4, ">".$NoNpe_2) or die "cannot open for writtinq! " . $!;
open (F5, ">".$outputWithN_NoN) or die "cannot open for writtinq! " . $!;

open (H1,  $pe_1) or die "cannot open for reading!" . $!;
open (H2,  $pe_2) or die "cannot open for reading!" . $!;

while (my $h1_p1 = <H1>) { # read first header of PE_1
  my $s_p1  = <H1>;   # read sequence of PE_1
  my $h2_p1 = <H1>;   # read second header of PE_1
  my $q_p1  = <H1>;   # read quality scores of PE_1
  my $h1_p2 = <H2>;   # read first header of PE_2
  my $s_p2 = <H2>;    # read sequence of PE_2
  my $h2_p2 = <H2>;   # read second header of PE_2
  my $q_p2 = <H2>;    # read quality scores of PE_2

  my $n_p1 = $s_p1 =~ tr/N//; # count the number of N's in read 1
  my $n_p2 = $s_p2 =~ tr/N//; # count the number of N's in read 2

  $totalN_p1=$totalN_p1+$n_p1;
  $totalN_p2=$totalN_p2+$n_p2;

  if($n_p1 > 0){
    $totalReadN_p1++;
    }
  if($n_p2 > 0){
    $totalReadN_p2++;
    }
  if($n_p1 > 0 || $n_p2 >0){## WIthN
    print F1 $h1_p1.$s_p1.$h2_p1.$q_p1;
    print F3 $h1_p2.$s_p2.$h2_p2.$q_p2;
    }
  else{ ## NoN
    print F2 $h1_p1.$s_p1.$h2_p1.$q_p1;
    print F4 $h1_p2.$s_p2.$h2_p2.$q_p2;
    }
  }

print F5 "Total N (read 1): ".$totalN_p1."\n";
print F5 "Total reads with at least an N (read 1): ".$totalReadN_p1."\n";
print F5 "Total N (read 2): ".$totalN_p2."\n";
print F5 "Total reads with at least an N (read 2): ".$totalReadN_p2."\n";

close H1;
close H2;
close F1;
close F2;
close F3;
close F4;
close F5;

#### Execute Sickle ####
my $SICKLE="sickle pe -f $NoNpe_1 -r $NoNpe_2 -t sanger -o $sickle_pe_1 -p $sickle_pe_2 -s $single_out -q $qual -l $len>> $output_Sickle";
print $SICKLE."\n";
system($SICKLE);
exit;
