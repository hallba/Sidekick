#!/usr/bin/perl
# Please use the -h option for help\n";
# This is a perl script that generates coarse-grained topology for proteins.
# The coarse-grained model is based on SJM's lipid/aminoacid MARTINI forcefield
# Author: Senthil K. Kandasamy
#
# 27.3.2008
#   Added support for using elastic networks instead of dihedrals for extended 
#   segments. -- Martti Louhivuori (m.j.louhivuori@rug.nl)
#
#
use warnings;
#use Time::local;
$force_field = "ffsjm";
$this_file_name = "seq2itp_martini2.1-elastic.pl";
$version = "2.1tryout";
$author = "Senthil K. Kandasamy";
$protein_name = "Protein";



####Define Options hash#######
%options_hash = (
"-seq" => "sequence.seq",
"-ssd" => "ssdump.ssd",
"-par" => "changepar.par",
"-itp" => "protein_cgmartini.itp",
"-cys" => "cysteine.cys",
"-out" => "outpar.out",
"-pro" => "no",
"-elastic" => 0,
);

###########################
# The following section reads and validates the input arguments
$numargs = scalar @ARGV;
#print ("The number of arguments is $numargs \n");
if ($numargs >14){
    print (" Sorry!! Too many arguments.\n");
    print (" Use the -h option for help\n");
    print"\nTOPOLOGY NOT GENERATED\n";
    exit ;
    }

# Check if help is requested
if ((grep $_ eq -h, @ARGV) >=1)
{
print("\n################################\nProgram : seq2cgtop Version $version \n################################\n");
print ("This is a Perl script. This assumes that you are familiar with  GROMACS \nand the MARTINI CG force field. This script generates a topology file for \ncoarse-grained proteins using SJM's MARTINI force field.\n");
print("\n#########\n# Input #\n#########\n");
print("The script REQUIRES two input files.\n
1. The aminoacid sequence (with a .seq extension)  of the protein, in FASTA \nformat. This can be downloaded from the protein data bank.\n
2. The secondary structure of the aminoacid residues in SSDUMP format (with a \n.ssd extension). The SSDUMP file can be generated using the do_dssp program \nof the GROMACS Suite\n");
print("\n##########\n# Output #\n##########\n");

print("By default, two output files will be generated\n
1. The topology file (default name is protein_cgmartini.itp. This can be changed \nby using the -itp option.) \n
2. A file containing all the parameters used (outpar.out).\n");
print("\n############\n# Optional #\n############\n");

print("1. The default parameters can be changed by using the -par option. This will \nrequire an input file with the list of parameters that need to be changed \n(default name is changepar.par). The outpar.out file can be used as a template \nfor the changepar.par file.\n\n");
print ("2. Disulfide bridges can be built using the -cys option \n\n");
print("------------------------------------------------------------\n");
print (" Option             Filename  Type         Description\n");
print ("------------------------------------------------------------\n");
print (" -seq             sequence.seq  Input        FASTA file\n");
print (" -ssd               ssdump.ssd  Input        SSDUMP file \n");
print (" -itp    protein_cgmartini.itp  Output       Topolgy file \n");
print (" -out               outpar.out  Output       Parameter file \n");
print (" -par            changepar.par  Input, Opt.  Parameter file \n");
print (" -cys             cysteine.cys  Input, Opt.  Disulfide file\n");
print (" -elastic                                    Use elastic networks\n");
print (" -h                                          Help!!\n");
 print ("------------------------------------------------------------\n");


print ("seq2cgtop -seq sequence.seq -ssd  ssdump.ssd -itp topology.itp\n\n");
exit;
}

### First, make sure that only one file is provided per input option
if ((@ARGV != 0)){
if  (!(exists $options_hash{$ARGV[0]})){
        print "\nInput Error 1: Unknown option $ARGV[0]\n";
        print "\nCheck your options again\n";
        print "\nUse -h for help\n";
        print"\nTOPOLOGY NOT GENERATED\n";
        exit;
        }
for ( my $i=1; $i< @ARGV; $i++){
        if ((!(exists $options_hash{$ARGV[$i]})) && (!(exists $options_hash{$ARGV[$i-1]}))){
                print "\nInput Error : Unknown option $ARGV[$i]\n";
                print "\nUse -h for help\n";
                print"\nTOPOLOGY NOT GENERATED\n";
                exit;
                }
        }
}

#### Now check to make sure that each option is passed only once
while ( (my $key, my $value) = each (%options_hash)){
        my $number_of_times = 0;
        for (my $i=0; $i< @ARGV; $i++){
                if ($key eq $ARGV[$i]){
                        $number_of_times++;
                        }
                }
        if ($number_of_times >1){
                print"\nFATAL error: You used option $key $number_of_times times. \nPlease use $key once only.";
                print "\nUse -h for help.\n";
                print"\nTOPOLOGY NOT GENERATED\n";
                exit;
                }
        }

#### check for sequence file
$index=undef;
if((grep $_ eq "-seq", @ARGV)<1){
        print "\n-seq option not specified. Will use the default file   $options_hash{-seq}\n";
        $filename_seq = $options_hash{-seq};
        }
elsif((grep $_ eq "-seq", @ARGV)==1){
        for(my $i = 0; $i < @ARGV; $i++ ){
                if( $ARGV[$i] eq "-seq"){
                        $index = $i;
                        last;
                        }
                }
        if ( $ARGV[-1] eq "-seq"){
                $filename_seq = $options_hash{-seq};
                print"Sequence file not specified. Will use the default file $options_hash{-seq}\n";
                }
        else{
                if ((exists $options_hash{$ARGV[$index+1]})){
                        $filename_seq = $options_hash{-seq};
                        print"Sequence file not specified. I will use the default file $options_hash{-seq}\n";
                        }
                else{
                        if (substr($ARGV[$index+1], -4,4) eq ".seq"){
                                $filename_seq = $ARGV[$index+1];
                                }
                        else{
                                $filename_seq =${ARGV[$index+1]}.".seq";
                                }
                        print"I will use the sequence file $filename_seq\n";
                        }
                }
        }
if(!open FILESEQ, $filename_seq){
        die "Fatal Error: Sequence file $filename_seq not found ($!)";
        }


##### Check for ssd file
$index=undef;
if((grep $_ eq "-ssd", @ARGV)<1){
        print "\n-ssd option not specified. Will use the default file $options_hash{-ssd}\n";
        $filename_ssd = $options_hash{-ssd};
        }
elsif((grep $_ eq "-ssd", @ARGV)==1){
        for(my $i = 0; $i < @ARGV; $i++ ){
                if( $ARGV[$i] eq "-ssd"){
                        $index = $i;
                        last;
                        }
                }
        if ( $ARGV[-1] eq "-ssd"){
                $filename_ssd = $options_hash{-ssd};
                print"Secondary Structure file not specified. Will use the default file $options_hash{-ssd}\n";
                }
        else{
                if ((exists $options_hash{$ARGV[$index+1]})){
                        $filename_ssd = $options_hash{-ssd};
                        print"Secondary Structure file not specified. I will use the default file $options_hash{-ssd}\n";
                                }
                        else{
                                if (substr($ARGV[$index+1], -4,4) eq ".ssd"){
                                $filename_ssd = $ARGV[$index+1];
                                        }
                                else{
                                        $filename_ssd =${ARGV[$index+1]}.".ssd";
                                        }
                                print"I will use the secondary structure file $filename_ssd\n";

                                }
                        }
                }
        if(!open FILESSD, $filename_ssd){
        die "Fatal Error: Secondary Structure file $filename_ssd not found ($!)";
        }

#### Check for -par option
$index=undef;
if((grep $_ eq "-par", @ARGV)<1){
        print "\n-par option not specified. Will use the default parameters\n";
        $option_par ="no";
        }
elsif((grep $_ eq "-par", @ARGV)==1){
        for(my $i = 0; $i < @ARGV; $i++ ){
                if( $ARGV[$i] eq "-par"){
                        $index = $i;
                        last;
                        }
                }
        if ( $ARGV[-1] eq "-par"){
                print"\nError: -par option is used, but parameter file not specified\n";
                print"\nPlease provide a parameter file or turn off the -par option\n";
                print"\nTOPOLOGY NOT GENERATED\n";
                exit;
                }
        else{
                if ((exists $options_hash{$ARGV[$index+1]})){
                        print"\nFATAL Error: -par option used, but parameter file not specified\n";
                print"\nPlease provide a parameter file or turn off the -par option\n";
                print"\nTOPOLOGY NOT GENERATED\n";
                exit;
                        }
                        else{
                                if (substr($ARGV[$index+1], -4,4) eq ".par"){
                                        $filename_par = $ARGV[$index+1];
                                        $option_par="yes";
                                        }
                                else{
                                        $filename_par =${ARGV[$index+1]}.".par";
                                        $option_par="yes";
                                        }
                                print"I will use the parameter file $filename_par\n";
                            print"\n\nCAUTION : SOME PARAMETERS WILL BE CHANGED\n\n";
                                }
                        }
                        if(!open FILEPAR, $filename_par){
        die "Fatal Error: Parameter file $filename_par not found ($!)";
        }
                }

#### Now, the -cys option
        $index=undef;
if((grep $_ eq "-cys", @ARGV)<1){
        print "\n-cys option not specified. Will not build any bridges\n";
        $option_cys ="no";
        }
elsif((grep $_ eq "-cys", @ARGV)==1){
        for(my $i = 0; $i < @ARGV; $i++ ){
                if( $ARGV[$i] eq "-cys"){
                        $index = $i;
                        last;
                        }
                }
        if ( $ARGV[-1] eq "-cys"){
                print"\nError: -cys option is used, but a file specifying cysteine bridge residues not provided\n";
                print"\nPlease provide a .cys file or turn off the -cys option\n";
                print"\nTOPOLOGY NOT GENERATED\n";
                exit;
                }
        else{
                if ((exists $options_hash{$ARGV[$index+1]})){
                        print"\nFATAL Error: -cys option used, but .cys file not provided\n";
                print"\nPlease provide a file specifying cysteine bridges or turn off the -cys option\n";
                print"\nTOPOLOGY NOT GENERATED\n";
                exit;
                        }
                        else{
                                if (substr($ARGV[$index+1], -4,4) eq ".cys"){
                                        $filename_cys = $ARGV[$index+1];
                                        $option_cys="yes";
                                        }
                                else{
                                        $filename_cys =${ARGV[$index+1]}.".cys";
                                        $option_cys="yes";
                                        }
                                print"I will use the cysteine-bridge file $filename_cys\n";
                                if(!open FILECYS, $filename_cys){
        die "Fatal Error: File specifying cysteine bridges, $filename_cys not found ($!)";
                                }
                        }
                }

        }

### Elastic network?
if ((grep $_ eq "-elastic", @ARGV)==1) {
	$ELASTIC_NETWORK = 1;
	print "\n-elastic option specified. Using elastic networks.\n";
} else {
	$ELASTIC_NETWORK = 0;
	print "\n-elastic option not specified. Using dihedrals.\n";
}

#### Now, the -itp option
$index=undef;
if((grep $_ eq "-itp", @ARGV)<1){
        print "\n-itp option not specified.\n";
        print"\n####################################################################\n";
        print"Will write the output to the default filename $options_hash{-itp}";
        $filename_itp=$options_hash{-itp};
        print"\n####################################################################\n";
        }
elsif((grep $_ eq "-itp", @ARGV)==1){
        for(my $i = 0; $i < @ARGV; $i++ ){
                if( $ARGV[$i] eq "-itp"){
                        $index = $i;
                        last;
                        }
                }
        if ( $ARGV[-1] eq "-itp"){
                print"\nWarning: -itp option is used, but a filename is not provided\n";
                print"\nI will write the output to the default filename $options_hash{-itp}\n";
                $filename_itp=$options_hash{-itp};
                }
        else{
                if ((exists $options_hash{$ARGV[$index+1]})){
                        print"\nWarning: -itp option is used, but a filename is not specified\n";
                print"\nI will write the output to the default filename $options_hash{-itp}\n\n";
                $filename_itp=$options_hash{-itp};
                        }
                        else{
                                if (substr($ARGV[$index+1], -4,4) eq ".itp"){
                                        $filename_itp = $ARGV[$index+1];
                                        }
                                else{
                                        $filename_itp =${ARGV[$index+1]}.".itp";
                                        }
                        print"\n####################################################\n";
                                print"\nI will write the output to  $filename_itp\n";
                                print"\n####################################################\n";
                        }
                }

        }
#### finally, the -out option
$index=undef;
if((grep $_ eq "-out", @ARGV)<1){
        print "\n-out option not specified. Will write all the parameters to the default filename $options_hash{-out}\n\n";
        $filename_out=$options_hash{-out};
        }
elsif((grep $_ eq "-out", @ARGV)==1){
        for(my $i = 0; $i < @ARGV; $i++ ){
                if( $ARGV[$i] eq "-out"){
                        $index = $i;
                        last;
                        }
                }
        if ( $ARGV[-1] eq "-out"){
                print"\nWarning: -out option is used, but a filename is not provided\n";
                print"I will write all the parameters to the default filename $options_hash{-out}\n";
                $filename_out=$options_hash{-out};
                }
        else{
                if ((exists $options_hash{$ARGV[$index+1]})){
                        print"\nWarning: -out option is used, but a filename is not specified";
                print"\nI will write the parameters to the default filename $options_hash{-out}\n\n";
                        $filename_out=$options_hash{-out};

                        }
                        else{
                                if (substr($ARGV[$index+1], -4,4) eq ".out"){
                                        $filename_out = $ARGV[$index+1];
                                        }
                                else{
                                        $filename_out =${ARGV[$index+1]}.".out";
                                        }
                                print"\nI will write the parameters to  $filename_out\n";
                        }
                }

        }


# Now, check if any of the protonation states needs to be changed
if((grep $_ eq "-pro", @ARGV)==1){
print"-pro option specified. How many residues need to be changed?\n";
chomp (my $protonation_changes = <STDIN>);
        }
###########################################################################
###########################################################################
########Hopefully all the arguments have been correctly parsed
##########################################################################
##########################################################################
# complete the help file later...
##########################

### IF the seq and ssd files do not exist, the program must have quit by now
### Now, read the contents of the seq and ssd files
### Sequence File First

open (FILESEQ,$filename_seq);
@lines_seq = <FILESEQ>;
close (FILESEQ);

$numlines_seq = scalar @lines_seq;
$temp_seq ="";
$complete_seq="";
for ($i=0; $i<$numlines_seq; $i++)
        {
        if(!( $lines_seq[$i] =~ /^>/))
                {
                chomp ($lines_seq[$i]);
                $complete_seq = $temp_seq.$lines_seq[$i];
                }
                $temp_seq=$complete_seq;
        }
#print $complete_seq;
@complete_sequence=split (//, $complete_seq);
$total_residues = scalar @complete_sequence;
print"\n";

##### The sequence (FASTA) file has been read , the comment lines have been
#### ignored, the different chains merged into one array, the array has been #### split . The sequence is finally stored in @complete_sequence

##### NOW check if there are any CYS residues and prompt for bridges.
my $number_of_cys=0;
for ($i=0;$i< $total_residues; $i++){
        if ($complete_sequence[$i] eq "C"){
                $number_of_cys++;
                }
        }
print "I found $number_of_cys CYS residues \n";
if (($number_of_cys >1) && ($option_cys eq "no")){
        print "\nPlease use the -cys option to build CYS bridges, if necessary\n\n";
        }

# Build Cysteine Bridges. Read the filename_cys
if($option_cys eq "yes"){
        open (FILECYS,$filename_cys);
        chomp(@lines_cys = <FILECYS>);
        close (FILECYS);

$number_of_bridges = scalar @lines_cys;
print"I will generate $number_of_bridges CYS bridges\n";
print "\nPLEASE MAKE SURE THAT YOU PROVIDE THE CORRECT CYSTEINE RESIDUES IN THE $filename_cys FILE OR ELSE THE TOPOLOGY WILL BE GARBAGE\n";
for (my $i =0; $i<$number_of_bridges; $i++){
        my @temp2arr =split (/\s+/,$lines_cys[$i]);
        if (@temp2arr!=2 ){
        print"\nSomething is wrong with the $filename_cys file. Please make sure that there are two and only two residues per line @temp2arr\n";
        print"\nTOPOLOGY NOT GENERATED\n";
        exit;
        }
        $begin_bridge[$i]=$temp2arr[0];
        $end_bridge[$i]=$temp2arr[1];
        if ($complete_sequence[$begin_bridge[$i]-1] ne "C"){
                print "\nHey, residue $begin_bridge[$i] is not a cysteine. Check your .cys file\n";
                print"TOPOLOGY NOT WRITTEN\n";
                exit;
                }
        if ($complete_sequence[$end_bridge[$i]-1] ne "C"){
                print "\nHey, residue $end_bridge[$i] is not a cysteine. Check your .cys file\n";
                print"TOPOLOGY NOT WRITTEN\n";
                exit;
                }
        print"CYS Bridge between residues $begin_bridge[$i] and $end_bridge[$i]\n";
        }
}
#### NOW, parse the SSDUMP file

open (FILESSD,$filename_ssd);
@lines_ssd = <FILESSD>;
close (FILESSD);
chomp ($lines_ssd[1]);
@complete_ssdump=split(//,$lines_ssd[1]);
my $total_ss = scalar @complete_ssdump;
for(my $i=0;$i< $total_ss; $i++){
        if($complete_ssdump[$i] eq "~"){
                $complete_ssdump[$i]="C";
                }
        }
for(my $i=0;$i< $total_ss; $i++){
        if($complete_ssdump[$i] eq "G"){
                $complete_ssdump[$i]="H";
                }
        }
for(my $i=0;$i< $total_ss; $i++){
        if($complete_ssdump[$i] eq "I"){
                $complete_ssdump[$i]="H";
                }
        }
        for(my $i=0;$i< $total_ss; $i++){
        if($complete_ssdump[$i] eq "B"){
                $complete_ssdump[$i]="E";
                }
        }
### The ssdump file is split and stored in @complete_ssdump

### Now Check to see if the number of residues in the seq and ssd files match
if ($total_residues != scalar @complete_ssdump){
        printf ("%s%d%s%d\n", "The number of residues in the sequence file, ", $total_residues, ", does not match the number of residues in the ssdump file, ", $total_ss);
        print"\nTOPOLOGY NOT GENERATED\n";
        exit;
        }



### Now, split into chains ...

open (FILECHN,$filename_seq);
@lines_chn = <FILECHN>;
close (FILECHN);

$numlines_chn = scalar @lines_chn;
$number_of_chains =0;

for(my $i= 0; $i< $numlines_chn; $i++){
        if ($lines_chn[$i] =~ /^>/){
                $number_of_chains++;
                }
        }
 $current_chain=0;
 $chain_length=0;
for(my $i= 0; $i< $numlines_chn; $i++){
        if ($lines_chn[$i] =~ /^>/){
                $each_chain[$current_chain]=$chain_length;
                $current_chain++;
                $chain_length=0;
                }
        else{
                #$current_chain_name = "chain$current_chain";
                chomp ($lines_chn[$i]);
                @current_line=split (//, $lines_chn[$i]);
                $chain_length+= scalar @current_line;
                }
        }
        $each_chain[$current_chain]=$chain_length;

print "Number of chains = $number_of_chains \n";
$current_length=0;
for($i=1;$i<=$current_chain;$i++) {
        $begin_chain[$i]= $current_length;
        $end_chain[$i]=$begin_chain[$i]+$each_chain[$i]-1;;
        $current_length=$end_chain[$i]+1;
        }
for ($i=1;$i<=$current_chain;$i++){
        printf ("%s%d%s%d%s%d%s%d\n","chain ",$i," has ",$each_chain[$i]," residues : from ",${begin_chain[$i]}+1," to ",$end_chain[$i]+1);
        }

########################### Consistency check done #################

#Define some hashes
#This is the main parameter listing hash. All the parameters are defined here. Edit this if a parameter is finalized. If just testing , use the -par option.
%parameter_hash = (
BEAD_BBN_HLX => N0,
BEAD_BBN_COI => P5,
BEAD_BBN_EXT => Nda,
BEAD_BBN_TRN => P1,
BEAD_BBN_BET => P1,
BEAD_BBN_BND => P1,
BEAD_BBN_5HX => N0,
BEAD_BBN_3HX => N0,
BEAD_BBN_NTR => Qd,
BEAD_BBN_CTR => Qa,
BEAD_BBH_CN1 => Nd,
BEAD_BBH_CN2 => Na,
BEAD_BBH_CN3 => Nda,
BEAD_ALA_HLX => C5,
BEAD_ALA_EXT => N0,
BEAD_ALA_COI => P4,
BEAD_ALA_TRN => C5,
BEAD_ALA_BND => C5,
BEAD_GLY_HLX => N0,
BEAD_GLY_EXT => Nda,
BEAD_GLY_COI => P5,
BEAD_GLY_TRN => P5,
BEAD_GLY_BND => P5,
BEAD_SC1_ASN => P5,
BEAD_SC1_ASP => Qa,
BEAD_SC1_GLU => Qa,
BEAD_SC1_GLN => P4,
BEAD_SC1_VAL => AC2,
BEAD_SC1_LEU => AC1,
BEAD_SC1_ILE => AC1,
BEAD_SC1_MET => C5,
BEAD_SC1_THR => P1,
BEAD_SC1_SER => P1,
BEAD_SC1_CYS => C5,
BEAD_SC1_PRO => AC2,
BEAD_SC1_LYS => C3,
BEAD_SC2_LYS => Qd,
BEAD_SC1_ARG => Nd,
BEAD_SC2_ARG => Qd,
BEAD_SC1_HIS => SC4,
BEAD_SC2_HIS => SP1,
BEAD_SC3_HIS => SP1,
BEAD_SC1_TYR => SC4,
BEAD_SC2_TYR => SC4,
BEAD_SC3_TYR => SP1,
BEAD_SC1_PHE => SC4,
BEAD_SC2_PHE => SC4,
BEAD_SC3_PHE => SC4,
BEAD_SC1_TRP => SNd,
BEAD_SC2_TRP => SC4,
BEAD_SC3_TRP => SC4,
BEAD_SC4_TRP => SC4,
BNLN_BBN_HLX => 0.35,
BNLN_BBN_COI => 0.35,
BNLN_BBN_EXT => 0.35,
BNLN_BBN_TRN => 0.35,
BNLN_BBN_BET => 0.35,
BNLN_BBN_BND => 0.35,
BNLN_BBN_5HX => 0.35,
BNLN_BBN_3HX => 0.35,
BNLN_SC1_ASN => 0.32,
BNLN_SC1_ASP => 0.32,
BNLN_SC1_GLU => 0.40,
BNLN_SC1_GLN => 0.40,
BNLN_SC1_VAL => 0.26,
BNLN_SC1_LEU => 0.33,
BNLN_SC1_ILE => 0.31,
BNLN_SC1_MET => 0.40,
BNLN_SC1_THR => 0.26,
BNLN_SC1_SER => 0.25,
BNLN_SC1_CYS => 0.31,
BNLN_SC1_PRO => 0.26,
BNLN_SC1_LYS => 0.33,
BNLN_SC2_LYS => 0.28,
BNLN_SC1_ARG => 0.33,
BNLN_SC2_ARG => 0.31,
BNLN_SC1_HIS => 0.32,
BNLN_TR1_HIS => 0.27,
BNLN_TR2_HIS => 0.27,
BNLN_TR3_HIS => 0.27,
BNLN_SC1_TYR => 0.32,
BNLN_TR1_TYR => 0.27,
BNLN_TR2_TYR => 0.27,
BNLN_TR3_TYR => 0.27,
BNLN_SC1_PHE => 0.31,
BNLN_TR1_PHE => 0.27,
BNLN_TR2_PHE => 0.27,
BNLN_TR3_PHE => 0.27,
BNLN_SC1_TRP => 0.3,
BNLN_TR1_TRP => 0.27,
BNLN_TR2_TRP => 0.27,
BNLN_TR3_TRP => 0.27,
BNLN_TR4_TRP => 0.27,
BNLN_TR5_TRP => 0.27,
BNLN_BRD_CYS => 0.4,
BNLN_HBN_HLX => 0.61,
BNKB_BBN_HLX => 1250,
BNKB_BBN_COI => 200,
BNKB_BBN_EXT => 1250,
BNKB_BBN_TRN => 500,
BNKB_BBN_BET => 400,
BNKB_BBN_BND => 500,
BNKB_BBN_5HX => 1250,
BNKB_BBN_3HX => 1250,
BNKB_SC1_ASN => 5000,
BNKB_SC1_ASP => 7500,
BNKB_SC1_GLU => 5000,
BNKB_SC1_GLN => 5000,
BNKB_SC1_VAL => 7500,
BNKB_SC1_LEU => 7500,
BNKB_SC1_ILE => 7500,
BNKB_SC1_MET => 2500,
BNKB_SC1_THR => 7500,
BNKB_SC1_SER => 7500,
BNKB_SC1_CYS => 7500,
BNKB_SC1_PRO => 7500,
BNKB_SC1_LYS => 5000,
BNKB_SC2_LYS => 5000,
BNKB_SC1_ARG => 5000,
BNKB_SC2_ARG => 5000,
BNKB_SC1_HIS => 7500,
BNKB_TR1_HIS => 1000,
BNKB_TR2_HIS => 1000,
BNKB_TR3_HIS => 1000,
BNKB_SC1_TYR => 5000,
BNKB_TR1_TYR => 1000,
BNKB_TR2_TYR => 1000,
BNKB_TR3_TYR => 1000,
BNKB_SC1_PHE => 7500,
BNKB_TR1_PHE => 1000,
BNKB_TR2_PHE => 1000,
BNKB_TR3_PHE => 1000,
BNKB_SC1_TRP => 5000,
BNKB_TR1_TRP => 1000,
BNKB_TR2_TRP => 1000,
BNKB_TR3_TRP => 1000,
BNKB_TR4_TRP => 1000,
BNKB_BRD_CYS => 5000,
BNKB_HBN_HLX => 1250,
BNKB_RNG_EQL => 1000,
ANGL_BBN_HLX => 96,
ANGL_BBN_COI => 127,
ANGL_BBN_EXT => 134,
ANGL_BBN_TRN => 120,
ANGL_BBN_BET => 134,
ANGL_BBN_BND => 120,
ANGL_BBN_5HX => 96,
ANGL_BBN_3HX => 96,
ANGL_SC1_ASN => 100,
ANGL_SC1_ASP => 100,
ANGL_SC1_GLU => 100,
ANGL_SC1_GLN => 100,
ANGL_SC1_VAL => 100,
ANGL_SC1_LEU => 100,
ANGL_SC1_ILE => 100,
ANGL_SC1_MET => 100,
ANGL_SC1_THR => 100,
ANGL_SC1_SER => 100,
ANGL_SC1_CYS => 100,
ANGL_SC1_PRO => 100,
ANGL_SC1_LYS => 100,
ANGL_SC2_LYS => 180,
ANGL_SC1_ARG => 100,
ANGL_SC2_ARG => 180,
ANGL_SC1_HIS => 100,
ANGL_SC2_HIS => 150,
ANGL_SC3_HIS => 150,
ANGL_SC1_TYR => 100,
ANGL_SC2_TYR => 150,
ANGL_SC3_TYR => 150,
ANGL_SC1_PHE => 100,
ANGL_SC2_PHE => 150,
ANGL_SC3_PHE => 150,
ANGL_SC1_TRP => 100,
ANGL_SC2_TRP => 90,
ANGL_SC3_TRP => 210,
ANKB_BBN_HLX => 700,
ANKB_BBN_COI => 25,
ANKB_BBN_EXT => 25,
ANKB_BBN_TRN => 25,
ANKB_BBN_BET => 25,
ANKB_BBN_BND => 25,
ANKB_BBN_5HX => 25,
ANKB_BBN_3HX => 25,
ANKB_SC1_ASN => 25,
ANKB_SC1_ASP => 25,
ANKB_SC1_GLN => 25,
ANKB_SC1_GLU => 25,
ANKB_SC1_VAL => 25,
ANKB_SC1_LEU => 25,
ANKB_SC1_ILE => 25,
ANKB_SC1_MET => 25,
ANKB_SC1_THR => 25,
ANKB_SC1_SER => 25,
ANKB_SC1_CYS => 25,
ANKB_SC1_LYS => 25,
ANKB_SC2_LYS => 25,
ANKB_SC1_PRO => 25,
ANKB_SC1_ARG => 25,
ANKB_SC2_ARG => 25,
ANKB_SC1_HIS => 25,
ANKB_SC2_HIS => 50,
ANKB_SC3_HIS => 50,
ANKB_SC1_TYR => 25,
ANKB_SC2_TYR => 50,
ANKB_SC3_TYR => 50,
ANKB_SC1_PHE => 25,
ANKB_SC2_PHE => 50,
ANKB_SC3_PHE => 50,
ANKB_SC1_TRP => 25,
ANKB_SC2_TRP => 50,
ANKB_SC3_TRP => 50,
DIAN_BBN_HLX => -120,
DIAN_BBN_EXT => 0,
DIAN_BBN_COI => 0,
DIAN_RG1_PHE => 0,
DIAN_RG1_HIS => 0,
DIAN_RG1_TYR => 0,
DIAN_RG1_TRP => 0,
DIAN_RG2_TRP => 0,
DIKB_BBN_HLX => 400,
DIKB_BBN_EXT => 10,
DIKB_BBN_COI => 10,
DIKB_RG1_PHE => 50,
DIKB_RG1_HIS => 50,
DIKB_RG1_TYR => 50,
DIKB_RG1_TRP => 50,
DIKB_RG2_TRP => 200,
DIKB_HLX_PRO => 100,
YES_BBN_CONST15 => "no",
YES_BBN_PROPDIH => "yes",
ELASTIC_SHORT => 0.64,
ELASTIC_LONG => 0.97,
ELASTIC_FORCE => 2500,
ITV_BONDS => "constraints",
RING_BONDS => "constraints",
Q_TERMINI => "CHARGED",
);

%aa_hash = ( G=>"GLY", A=>"ALA", D=>"ASP", N=>"ASN", E=>"GLU", Q=>"GLN", V=>"VAL", L=>"LEU",I=>"ILE", M=>"MET", T=>"THR", S=>"SER", C=>"CYS", K=>"LYS", R=>"ARG", H=>"HIS", F=> "PHE", P=>"PRO", W=>"TRP", Y=>"TYR" );

%ss_hash = (C=>"COI", H=>"HLX", S=>"BND", T=>"TRN", B=>"BET", G=>"3HX", I=>"5HX", E=>"EXT");

%q_hash = (G=>0, A=>0, D=>-1.0, N=>0, E=>-1.0, Q=>0, V=>0, L=>0,I=>0, M=>0, T=>0, S=>0, C=>0, K=>1.0, R=>1.0, H=>0, F=>0, P=>0, W=>0, Y=>0);

%numbeads_hash = (G=>1, A=>1, D=>2, N=>2, E=>2, Q=>2, V=>2, L=>2,I=>2, M=>2, T=>2, S=>2, C=>2, K=>3, R=>3, H=>4, F=>4, P=>2, W=>5, Y=>4);

%ssprecedence_hash = ("3HX"=>1, "HLX"=>2, "5HX"=>3, "TRN"=>4, "BND"=>5, "BET"=>6, "EXT"=>8, "COI"=>7);



##### if -par option is passed, read the parameters and change the values in the hash.

if ($option_par eq "yes"){

        open (FILEPAR, $filename_par);
        @lines_par = <FILEPAR>;
        close (FILEPAR);
        #print @lines_par;

        $number_of_parameters_changed = scalar @lines_par;
        print "\nWARNING : Parameters will be changed. Make sure that the new parameter values are meaningful. Otherwise, the topology file generated will be meaningless.\n\n";
        print "The number of paramters to be changed: $number_of_parameters_changed\n\n";
        for($i=0;$i<$number_of_parameters_changed;$i++){
                chomp($lines_par[$i]);
                my @this_line=split (/\s+/,$lines_par[$i]);
                if (exists $parameter_hash{$this_line[0]}){
                        my $old_parameter_value = $parameter_hash{$this_line[0]};
                        print "Parameter $this_line[0] will be changed from $old_parameter_value to $this_line[1]\n";
                        $parameter_hash{$this_line[0]}= $this_line[1];
                        }
                else{
                        print "\nFATAL ERROR:  The parameter, $this_line[0] does not exist. Please check for Typos. Use a parout.par file as a template.\n ";
                        print"\nTOPOLOGY NOT GENERATED\n";exit;
                        }
                }


        }
        if (($parameter_hash{ITV_BONDS} ne "constraints") && ($parameter_hash{ITV_BONDS} ne "bonds")) {
        print "\nFATAL ERROR:  The parameter ITV_BONDS should be either bonds or constraints. Please check for Typos.\n ";
                        print"\nTOPOLOGY NOT GENERATED\n";exit
                }
        if (($parameter_hash{RING_BONDS} ne "constraints") && ($parameter_hash{RING_BONDS} ne "bonds")) {
        print "\nFATAL ERROR:  The parameter RING_BONDS should be either bonds or constraints. Please check for Typos.\n ";
                        print"\nTOPOLOGY NOT GENERATED\n";exit
                }

######## DEFINE Subroutines

sub atoms_single_bead {
        $atoms1[$_[0]]=$_[0]+1;
        $atoms2[$_[0]]=         $parameter_hash{"BEAD_$aa_hash{${complete_sequence[($_[1])]}}_$ss_hash{${complete_ssdump[($_[1])]}}"};
        $atoms3[$_[0]]=$_[1]+1;
        $atoms4[$_[0]]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]]= "B${complete_ssdump[$_[1]]}$atoms2[$_[0]]";
        $atoms6[$_[0]]= $_[0]+1;
        $atoms7[$_[0]]= 0;
        $atoms8[$_[0]]=$ss_hash{${complete_ssdump[$_[1]]}};
        $backbone_atoms[$_[1]]=$_[0]+1;
        }

sub atoms_two_beads {
        $atoms1[$_[0]]=$_[0]+1;
        $atoms2[$_[0]]=         $parameter_hash{"BEAD_BBN_$ss_hash{${complete_ssdump[($_[1])]}}"};
        $atoms3[$_[0]]=$_[1]+1;
        $atoms4[$_[0]]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]]= "B${complete_ssdump[$_[1]]}$atoms2[$_[0]]";
        $atoms6[$_[0]]= $_[0]+1;
        $atoms7[$_[0]]= 0;
        $atoms8[$_[0]]=$ss_hash{${complete_ssdump[$_[1]]}};
        $backbone_atoms[$_[1]]=$_[0]+1;

        $atoms1[$_[0]+1]=$_[0]+2;
        $atoms2[$_[0]+1]=       $parameter_hash{"BEAD_SC1_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+1]=$_[1]+1;
        $atoms4[$_[0]+1]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+1]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+1]";
        $atoms6[$_[0]+1]= $_[0]+2;
        $atoms7[$_[0]+1]=$q_hash{${complete_sequence[($_[1])]}};
        $atoms8[$_[0]+1]=$ss_hash{${complete_ssdump[($_[1])]}};
}


sub atoms_three_beads {
        $atoms1[$_[0]]=$_[0]+1;
        $atoms2[$_[0]]=         $parameter_hash{"BEAD_BBN_$ss_hash{${complete_ssdump[($_[1])]}}"};
        $atoms3[$_[0]]=$_[1]+1;
        $atoms4[$_[0]]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]]= "B${complete_ssdump[$_[1]]}$atoms2[$_[0]]";
        $atoms6[$_[0]]= $_[0]+1;
        $atoms7[$_[0]]= 0;
        $atoms8[$_[0]]=$ss_hash{${complete_ssdump[$_[1]]}};
        $backbone_atoms[$_[1]]=$_[0]+1;

        $atoms1[$_[0]+1]=$_[0]+2;
        $atoms2[$_[0]+1]=       $parameter_hash{"BEAD_SC1_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+1]=$_[1]+1;
        $atoms4[$_[0]+1]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+1]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+1]";
        $atoms6[$_[0]+1]= $_[0]+2;
        $atoms7[$_[0]+1]=0;
        $atoms8[$_[0]+1]=$ss_hash{${complete_ssdump[($_[1])]}};

        $atoms1[$_[0]+2]=$_[0]+3;
        $atoms2[$_[0]+2]=       $parameter_hash{"BEAD_SC2_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+2]=$_[1]+1;
        $atoms4[$_[0]+2]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+2]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+2]";
        $atoms6[$_[0]+2]= $_[0]+3;
        $atoms7[$_[0]+2]=$q_hash{${complete_sequence[($_[1])]}};
        $atoms8[$_[0]+2]=$ss_hash{${complete_ssdump[($_[1])]}};
}

sub atoms_four_beads {
        $atoms1[$_[0]]=$_[0]+1;
        $atoms2[$_[0]]=         $parameter_hash{"BEAD_BBN_$ss_hash{${complete_ssdump[($_[1])]}}"};
        $atoms3[$_[0]]=$_[1]+1;
        $atoms4[$_[0]]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]]= "B${complete_ssdump[$_[1]]}$atoms2[$_[0]]";
        $atoms6[$_[0]]= $_[0]+1;
        $atoms7[$_[0]]= 0;
        $atoms8[$_[0]]=$ss_hash{${complete_ssdump[$_[1]]}};
        $backbone_atoms[$_[1]]=$_[0]+1;

        $atoms1[$_[0]+1]=$_[0]+2;
        $atoms2[$_[0]+1]=       $parameter_hash{"BEAD_SC1_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+1]=$_[1]+1;
        $atoms4[$_[0]+1]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+1]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+1]";
        $atoms6[$_[0]+1]= $_[0]+2;
        $atoms7[$_[0]+1]=0;
        $atoms8[$_[0]+1]=$ss_hash{${complete_ssdump[($_[1])]}};

        $atoms1[$_[0]+2]=$_[0]+3;
        $atoms2[$_[0]+2]=       $parameter_hash{"BEAD_SC2_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+2]=$_[1]+1;
        $atoms4[$_[0]+2]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+2]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+2]";
        $atoms6[$_[0]+2]= $_[0]+3;
        $atoms7[$_[0]+2]=0;
        $atoms8[$_[0]+2]=$ss_hash{${complete_ssdump[($_[1])]}};

        $atoms1[$_[0]+3]=$_[0]+4;
        $atoms2[$_[0]+3]=       $parameter_hash{"BEAD_SC3_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+3]=$_[1]+1;
        $atoms4[$_[0]+3]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+3]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+3]";
        $atoms6[$_[0]+3]= $_[0]+4;
        $atoms7[$_[0]+3]=$q_hash{${complete_sequence[($_[1])]}};
        $atoms8[$_[0]+3]=$ss_hash{${complete_ssdump[($_[1])]}};
}


sub atoms_five_beads {
        $atoms1[$_[0]]=$_[0]+1;
        $atoms2[$_[0]]=         $parameter_hash{"BEAD_BBN_$ss_hash{${complete_ssdump[($_[1])]}}"};
        $atoms3[$_[0]]=$_[1]+1;
        $atoms4[$_[0]]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]]= "B${complete_ssdump[$_[1]]}$atoms2[$_[0]]";
        $atoms6[$_[0]]= $_[0]+1;
        $atoms7[$_[0]]= 0;
        $atoms8[$_[0]]=$ss_hash{${complete_ssdump[$_[1]]}};
        $backbone_atoms[$_[1]]=$_[0]+1;

        $atoms1[$_[0]+1]=$_[0]+2;
        $atoms2[$_[0]+1]=       $parameter_hash{"BEAD_SC1_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+1]=$_[1]+1;
        $atoms4[$_[0]+1]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+1]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+1]";
        $atoms6[$_[0]+1]= $_[0]+2;
        $atoms7[$_[0]+1]=0;
        $atoms8[$_[0]+1]=$ss_hash{${complete_ssdump[($_[1])]}};

        $atoms1[$_[0]+2]=$_[0]+3;
        $atoms2[$_[0]+2]=       $parameter_hash{"BEAD_SC2_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+2]=$_[1]+1;
        $atoms4[$_[0]+2]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+2]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+2]";
        $atoms6[$_[0]+2]= $_[0]+3;
        $atoms7[$_[0]+2]=0;
        $atoms8[$_[0]+2]=$ss_hash{${complete_ssdump[($_[1])]}};

        $atoms1[$_[0]+3]=$_[0]+4;
        $atoms2[$_[0]+3]=       $parameter_hash{"BEAD_SC3_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+3]=$_[1]+1;
        $atoms4[$_[0]+3]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+3]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+3]";
        $atoms6[$_[0]+3]= $_[0]+4;
        $atoms7[$_[0]+3]=0;
        $atoms8[$_[0]+3]=$ss_hash{${complete_ssdump[($_[1])]}};

        $atoms1[$_[0]+4]=$_[0]+5;
        $atoms2[$_[0]+4]=       $parameter_hash{"BEAD_SC4_$aa_hash{${complete_sequence[($_[1])]}}"};
        $atoms3[$_[0]+4]=$_[1]+1;
        $atoms4[$_[0]+4]=$aa_hash{$complete_sequence[($_[1])]};
        $atoms5[$_[0]+4]= "S${complete_ssdump[$_[1]]}$atoms2[$_[0]+4]";
        $atoms6[$_[0]+4]= $_[0]+5;
        $atoms7[$_[0]+4]=$q_hash{${complete_sequence[($_[1])]}};
        $atoms8[$_[0]+4]=$ss_hash{${complete_ssdump[($_[1])]}};
}


sub ss_precedence {
        if ( ${ssprecedence_hash{$ss_hash{$complete_ssdump[$_[0]]}}} >
                ${ssprecedence_hash{$ss_hash{$complete_ssdump[$_[1]]}}}) {
                $ss_hash{$complete_ssdump[$_[0]]};
                }
        else {
                $ss_hash{$complete_ssdump[$_[1]]};
                }
        }

###########################
####### Print Stuff#########
### Now, start printing everything to the itp file

### First, dump the parameters used to the parout.par file

open FPAROUT, ">$filename_out";
foreach $key (sort keys %parameter_hash){
        printf FPAROUT ("%s\t%s\n",$key, $parameter_hash{$key});
        }
close FPAROUT;

open FITP,"> $filename_itp";
printf FITP ("%s%s\n", ";;; SJM topology ", $force_field);
printf FITP ("%s%s%s%s\n\n", ";;; Generated by ",$this_file_name," written by ", $author);
printf FITP ("%s", ";;; SEQ: ");
printf ("\n%s\n\n", "I will generate a topology file for the following sequence: ");

for ($i=0;$i< $total_residues; $i++)
        {
        printf FITP ("%s", $complete_sequence[$i]);
        printf ("%s", $complete_sequence[$i]);

        $count1++;
        if ($count1%10==0)
                {
                printf FITP ("%s", " ");
                printf ("%s", " ");

                }
        if ($count1%50==0)
                {
                printf FITP ("%s%d"," ",$count1);
                printf ("%s%d\n"," ",$count1);

                printf FITP ("\n%s",";;; SEQ: ");

                }
        }
printf FITP ("\n");
printf ("\n");
printf FITP ("%s%d\n",";;; Total Number of Amino Acid Residues: ",$total_residues);
printf("\n%s%d\n\n","Total number of amino acid residues: ",$total_residues);
printf FITP ("\n\n");
printf FITP ("%s\n","[moleculetype]");
printf FITP ("%s\t%s\n",";molname","exclusions");
printf FITP ("%s\t%d\n", $protein_name,1);
printf FITP ("\n\n");
my $atom_number=0;
for (my $i=0;$i<$total_residues; $i++){
        if ($numbeads_hash{$complete_sequence[$i]}==1){
                &atoms_single_bead($atom_number,$i);
                $atom_number++;
                }
        elsif ($numbeads_hash{$complete_sequence[$i]}==2){
                &atoms_two_beads($atom_number,$i);
                $atom_number+=2;
                }
        elsif ($numbeads_hash{$complete_sequence[$i]}==3){
                &atoms_three_beads($atom_number,$i);
                $atom_number+=3;
                }
        elsif ($numbeads_hash{$complete_sequence[$i]}==4){
                &atoms_four_beads($atom_number,$i);
                $atom_number+=4;
                }
        elsif ($numbeads_hash{$complete_sequence[$i]}==5){
                &atoms_five_beads($atom_number,$i);
                $atom_number+=5;
                }
        else {
                print "FATAL ERROR # 897 . Something is wrong\n";
                }
        }
# Now,check for chains, helical stretches and change the appropriate backbone bead types

for( my $i=1;$i<=$number_of_chains;$i++){
        $number_of_helices=0;
        @start_helix=undef;
        @end_helix=undef;
        if($complete_ssdump[$begin_chain[$i]] eq "H") {
                $number_of_helices++;
                $start_helix[$number_of_helices]=$begin_chain[$i];
                if($complete_ssdump[${begin_chain[$i]}+1] ne "H") {
                        $end_helix[$number_of_helices]=$begin_chain[$i];
                        }
                }
        for (my $j=${begin_chain[$i]}+1; $j<=${end_chain[$i]}-1; $j++){
                if (($complete_ssdump[$j] eq "H") &&($complete_ssdump[$j-1] ne "H")){
                        $number_of_helices++;
                        $start_helix[$number_of_helices]=$j;
                        }
                if (($complete_ssdump[$j] eq "H") &&($complete_ssdump[$j+1] ne "H")){
                        $end_helix[$number_of_helices]=$j;
                        }
                }
        if($complete_ssdump[$end_chain[$i]] eq "H") {
                if ($complete_ssdump[${end_chain[$i]}-1] eq "H"){
                        $end_helix[$number_of_helices]=$end_chain[$i];
                        }
            else{
                $number_of_helices++;
                        $end_helix[$number_of_helices]=$end_chain[$i];
                        $start_helix[$number_of_helices]=$end_chain[$i];;
                        }
                }
        for (my $k=1;$k<=$number_of_helices; $k++){
                $helix_length[$k]=$end_helix[$k]-$start_helix[$k]+1;
                #print "helix length is $helix_length[$k] from $start_helix[$k] to $end_helix[$k];\n";
                }
        # At this point, for this particular chain, all the helices have been evaluated. Now , go on and change the bead types : stored in @atoms2
         for (my $x=1;$x<=$number_of_helices;$x++){
                $first_helical_residue=$start_helix[$x];
                $last_helical_residue=$end_helix[$x];

                if($helix_length[$x]>=8){
                        $atoms2[${backbone_atoms[$first_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN1"};
                        $atoms2[${backbone_atoms[$first_helical_residue+1]}-1]= $parameter_hash{"BEAD_BBH_CN1"};
                        $atoms2[${backbone_atoms[$first_helical_residue+2]}-1]= $parameter_hash{"BEAD_BBH_CN1"};
                        $atoms2[${backbone_atoms[$first_helical_residue+3]}-1]= $parameter_hash{"BEAD_BBH_CN1"};
                        $atoms2[${backbone_atoms[$last_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        $atoms2[${backbone_atoms[$last_helical_residue-1]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        $atoms2[${backbone_atoms[$last_helical_residue-2]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        $atoms2[${backbone_atoms[$last_helical_residue-3]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        }
                elsif($helix_length[$x]==7){
                        $atoms2[${backbone_atoms[$first_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN1"};
                        $atoms2[${backbone_atoms[$first_helical_residue+1]}-1]= $parameter_hash{"BEAD_BBH_CN1"};
                        $atoms2[${backbone_atoms[$first_helical_residue+2]}-1]= $parameter_hash{"BEAD_BBH_CN1"};
                        $atoms2[${backbone_atoms[$first_helical_residue+3]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        $atoms2[${backbone_atoms[$last_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        $atoms2[${backbone_atoms[$last_helical_residue-1]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        $atoms2[${backbone_atoms[$last_helical_residue-2]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        }
                elsif($helix_length[$x]==6){
                        $atoms2[${backbone_atoms[$first_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN1"};
                        #$atoms5[${backbone_atoms[$first_helical_residue]}-1]= "B".${complete_ssdump[$first_helical_residue]}.$parameter_hash{BEAD_BBH_CN1};
                        $atoms2[${backbone_atoms[$first_helical_residue+1]}-1]= $parameter_hash{"BEAD_BBH_CN1"};
                        #$atoms5[${backbone_atoms[$first_helical_residue+1]}-1]= "B".${complete_ssdump[$first_helical_residue+1]}.$parameter_hash{BEAD_BBH_CN1};
                        $atoms2[${backbone_atoms[$first_helical_residue+2]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        #$atoms5[${backbone_atoms[$first_helical_residue+2]}-1]= "B".${complete_ssdump[$first_helical_residue+2]}.$parameter_hash{BEAD_BBH_CN3};
                        $atoms2[${backbone_atoms[$last_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        #$atoms5[${backbone_atoms[$last_helical_residue]}-1]= "B".${complete_ssdump[$last_helical_residue]}.$parameter_hash{BEAD_BBH_CN2};
                        $atoms2[${backbone_atoms[$last_helical_residue-1]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        #$atoms5[${backbone_atoms[$last_helical_residue-1]}-1]= "B".${complete_ssdump[$last_helical_residue-1]}.$parameter_hash{BEAD_BBH_CN2};
                        $atoms2[${backbone_atoms[$last_helical_residue-2]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        #$atoms5[${backbone_atoms[$last_helical_residue-2]}-1]= "B".${complete_ssdump[$last_helical_residue-2]}.$parameter_hash{BEAD_BBH_CN3};
                        }
                elsif($helix_length[$x]==5){
                        $atoms2[${backbone_atoms[$first_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN1"};


                        $atoms2[${backbone_atoms[$first_helical_residue+1]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        $atoms2[${backbone_atoms[$first_helical_residue+2]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        $atoms2[${backbone_atoms[$last_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN2"};
                        $atoms2[${backbone_atoms[$last_helical_residue-1]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        }
                # The rest are probably redundant and insignificant and probably impossible according to dssp definitions, but for sake of completeness.....
                elsif($helix_length[$x]==4){
                        $atoms2[${backbone_atoms[$first_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        $atoms2[${backbone_atoms[$first_helical_residue+1]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        $atoms2[${backbone_atoms[$last_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        $atoms2[${backbone_atoms[$last_helical_residue-1]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        }
                elsif($helix_length[$x]==3){
                        $atoms2[${backbone_atoms[$first_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        $atoms2[${backbone_atoms[$first_helical_residue+1]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        $atoms2[${backbone_atoms[$last_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        }
                elsif($helix_length[$x]==2){
                        $atoms2[${backbone_atoms[$first_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        $atoms2[${backbone_atoms[$last_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        }
                elsif($helix_length[$x]==1){
                        $atoms2[${backbone_atoms[$first_helical_residue]}-1]= $parameter_hash{"BEAD_BBH_CN3"};
                        }
                }
         }

### Nd, Na, Nda assignment for peripheral helical backbones is complete
### Now, change the Na and Nda back to N0 if it is an alanine residue
for (my $i=0;$i<$atom_number;$i++){
        if (($atoms4[$i] eq "ALA") &&(($atoms2[$i] eq "Na")||($atoms2[$i] eq "Nd") ||($atoms2[$i] eq "Nda"))){
                #print "wazzup\n";
                $atoms2[$i]="N0";
                }
        }


### NOW CHECK FOR CYS BRIDGES

### Now, change the termini to be charged for each chain
if ($parameter_hash{"Q_TERMINI"} eq "CHARGED"){
        for (my $j=1;$j<=$number_of_chains; $j++){
$atoms2[${backbone_atoms[$begin_chain[$j]]}-1]="Qd";
$atoms2[${backbone_atoms[$end_chain[$j]]}-1]="Qa";
$atoms7[${backbone_atoms[$begin_chain[$j]]}-1]=1.0;
$atoms7[${backbone_atoms[$end_chain[$j]]}-1]=-1.0;

                }
        }
        else{
        print "\nCAUTION : The termini are not charged. I hope this is what you wanted\n";
        }
### Change   the atoms names for the atoms now. There are probably several unnecessary atoms5[]lines in the atoms subroutines as the following loop will overwrite everything. Clean-up later...

for (my $j=1;$j<=$number_of_chains; $j++){
        for ( my $i=$begin_chain[$j]; $i<=$end_chain[$j];$i++) {
$atoms5[${backbone_atoms[$i]}-1]="B${complete_ssdump[$i]}$atoms2[${backbone_atoms[$i]}-1]";
                }
        }


##### Now that all the  necessary changes are made, write the [atoms] directive
#####################################ATOMS#####################################
printf FITP ("%s\n","[atoms]");

for (my $i=0;$i<$atom_number;$i++){
        printf FITP ("%d\t%s\t%d\t%s\t%s\t%d\t%1.3f\t%s\t%s\n",$atoms1[$i],$atoms2[$i],$atoms3[$i],$atoms4[$i],$atoms5[$i],$atoms6[$i],$atoms7[$i],";",$atoms8[$i]);
        }
##### [atoms] done

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#####################################ATOMS DONE#################################

#####################################BONDS#####################################
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#####write the bonds [bonds]  directive

printf FITP ("\n%s\n","[ bonds ]");
printf FITP ("%s\n",";backbone-backbone bonds");

for (my $j=1;$j<=$number_of_chains; $j++){
        for ( my $i=$begin_chain[$j]; $i<$end_chain[$j];$i++) {
                $this_bond_ss = &ss_precedence ($i, $i+1);
                #printf "$this_bond_ss\n";
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s\n",$backbone_atoms[$i],$backbone_atoms[$i+1],1,$parameter_hash{"BNLN_BBN_$this_bond_ss"},$parameter_hash{"BNKB_BBN_$ss_hash{${complete_ssdump[$i]}}"},";","$ss_hash{${complete_ssdump[$i]}}"."-"."$ss_hash{${complete_ssdump[$i+1]}}");

                }
        }
printf FITP ("%s\n",";bb-sc bonds");

if ($parameter_hash{"ITV_BONDS"} eq "constraints"){
for ( my $i=0; $i<$total_residues;$i++){
        if (($numbeads_hash{$complete_sequence[$i]} != 1) &&  ($complete_sequence[$i] ne "I") && ($complete_sequence[$i] ne "T") && ($complete_sequence[$i] ne "V"))
        {
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%s\n",$backbone_atoms[$i],${backbone_atoms[$i]}+1,1,$parameter_hash{"BNLN_SC1_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_SC1_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                }
        }
}

elsif ($parameter_hash{"ITV_BONDS"} eq "bonds"){

printf ("%s\n","Generating bb-sc bonds instead of constraints for T I, P\n");
for ( my $i=0; $i<$total_residues;$i++){
        if ($numbeads_hash{$complete_sequence[$i]} != 1)
        {
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%s\n",$backbone_atoms[$i],${backbone_atoms[$i]}+1,1,$parameter_hash{"BNLN_SC1_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_SC1_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                }
        }
}
printf FITP ("%s\n",";sc1-sc2 bonds (ARG, LYS)");
for ( my $i=0; $i<$total_residues;$i++){
        if ($numbeads_hash{$complete_sequence[$i]} ==3 ){
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+2,1,$parameter_hash{"BNLN_SC2_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_SC2_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                }
        }


### NOW, if RINGS_BONDS is yes, generate bonds instead of constraints for ring structures. This option will be useful for equilibration simulations.
if ($parameter_hash{"RING_BONDS"} eq "bonds"){
printf FITP ("%s\n",";sc-sc bonds for RINGS (TRP, TYR, PHE, HIS)");
printf ("%s\n","Generating sc-sc bonds instead of constraints for RINGS (TRP, TYR, PHE, HIS)");

############################
for ( my $i=0; $i<$total_residues;$i++){

        # FIRST, HIS, PHE and TYR
        if ($numbeads_hash{$complete_sequence[$i]} ==4 ){
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+2,1,$parameter_hash{"BNLN_TR1_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_RNG_EQL"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+3,1,$parameter_hash{"BNLN_TR2_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_RNG_EQL"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+2,${backbone_atoms[$i]}+3,1,$parameter_hash{"BNLN_TR3_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_RNG_EQL"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                }

        #       Now, TRP
        if ($numbeads_hash{$complete_sequence[$i]} ==5 ){
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+2,1,$parameter_hash{"BNLN_TR1_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_RNG_EQL"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+3,1,$parameter_hash{"BNLN_TR2_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_RNG_EQL"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+2,${backbone_atoms[$i]}+3,1,$parameter_hash{"BNLN_TR3_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_RNG_EQL"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+2,${backbone_atoms[$i]}+4,1,$parameter_hash{"BNLN_TR3_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_RNG_EQL"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+3,${backbone_atoms[$i]}+4,1,$parameter_hash{"BNLN_TR3_$aa_hash{${complete_sequence[$i]}}"},$parameter_hash{"BNKB_RNG_EQL"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                }
        }
############################
}


### now, cysteine bridges, yipeee
if ($option_cys eq "yes"){
        printf FITP ("%s\n",";CYS-CYS bonds");
        for( my $i =0; $i<$number_of_bridges; $i++){
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s%d\n",${backbone_atoms[$begin_bridge[$i]-1]}+1,${backbone_atoms[$end_bridge[$i]-1]}+1,1,$parameter_hash{"BNLN_BRD_$aa_hash{${complete_sequence[$begin_bridge[$i]-1]}}"},$parameter_hash{"BNKB_BRD_$aa_hash{${complete_sequence[$end_bridge[$i]-1]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                }
        }

#Now build the 1-5 bonds, if necessary
if ($parameter_hash{YES_BBN_CONST15} eq "yes"){
        printf FITP ("%s\n",";1-5 HB constraints");
        for (my $j=1;$j<=$number_of_chains; $j++){
        for ( my $i=$begin_chain[$j]; $i<$end_chain[$j]-3;$i++) {
                if(($complete_ssdump[$i] eq "H" ) && ($complete_ssdump[$i+1] eq "H" ) && ($complete_ssdump[$i+2] eq "H" ) && ($complete_ssdump[$i+3] eq "H" ) &&($complete_ssdump[$i+4] eq "H")){
                #printf "$this_bond_ss\n";
                printf FITP ("%d\t%d\t%d\t%1.3f\t%d\t%s\t%s\n",$backbone_atoms[$i],$backbone_atoms[$i+4],1,$parameter_hash{"BNLN_HBN_HLX"},$parameter_hash{"BNKB_HBN_HLX"},";","$ss_hash{${complete_ssdump[$i]}}"."-"."$ss_hash{${complete_ssdump[$i+4]}}");
                }}
        }
        }

# Elastic network instead of dihedrals?
if ($ELASTIC_NETWORK) {
	my %short_bonds = ();
	my %long_bonds = ();
	for (my $j=1; $j<=$number_of_chains; $j++) {
		for (my $i=$begin_chain[$j]+1; $i<$end_chain[$j]-1; $i++) {
			if ( ($complete_ssdump[$i-1] eq "E") 
					&& ($complete_ssdump[$i] eq "E") 
					&& ($complete_ssdump[$i+1] eq "E") 
					&& ($complete_ssdump[$i+2] eq "E") ) {
				# atom and residue details
				($atom0, $res0) = (sprintf("%08d", $backbone_atoms[$i - 1]), 
					join('', "$aa_hash{${complete_sequence[$i - 1]}}", $i));
				($atom1, $res1) = (sprintf("%08d", $backbone_atoms[$i]), 
					join('', "$aa_hash{${complete_sequence[$i]}}", $i + 1));
				($atom2, $res2) = (sprintf("%08d", $backbone_atoms[$i + 1]), 
					join('', "$aa_hash{${complete_sequence[$i + 1]}}", $i + 2));
				($atom3, $res3) = (sprintf("%08d", $backbone_atoms[$i + 2]), 
					join('', "$aa_hash{${complete_sequence[$i + 2]}}", $i + 3));
				# store first to avoid duplicates
				$short_bonds{"$atom0 $atom2"} = "${res0}-${res2}";
				$short_bonds{"$atom1 $atom3"} = "${res1}-${res3}";
				$long_bonds{"$atom0 $atom3"} = "${res0}-${res3}";
			}
		}
	}
	# print short range elastic bonds
	printf FITP (";short elastic bonds\n");
	foreach $key (sort (keys %short_bonds)) {
		my ($i, $j) = split(/ /, $key);
		printf FITP ("%d\t%d\t1\t%1.2f\t%d ; %s\n", $i, $j, 
			$parameter_hash{'ELASTIC_SHORT'}, $parameter_hash{'ELASTIC_FORCE'}, 
			$short_bonds{$key});
	}
	# print long range elastic bonds
	printf FITP (";long elastic bonds\n");
	foreach $key (sort (keys %long_bonds)) {
		my ($i, $j) = split(/ /, $key);
		printf FITP ("%d\t%d\t1\t%1.2f\t%d ; %s\n", $i, $j, 
			$parameter_hash{'ELASTIC_LONG'}, $parameter_hash{'ELASTIC_FORCE'}, 
			$long_bonds{$key});
	}
}


### NOW, print the triangle constraints
printf FITP ("\n%s\n","[ constraints ]");


if ($parameter_hash{"RING_BONDS"} eq "constraints"){
printf FITP ("%s\n",";sc-sc constraints (Ring Structures)");
for ( my $i=0; $i<$total_residues;$i++){

        # FIRST, HIS, PHE and TYR
        if ($numbeads_hash{$complete_sequence[$i]} ==4 ){
                printf FITP ("%d\t%d\t%d\t%1.3f\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+2,1,$parameter_hash{"BNLN_TR1_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+3,1,$parameter_hash{"BNLN_TR2_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%s\t%s%d\n",${backbone_atoms[$i]}+2,${backbone_atoms[$i]}+3,1,$parameter_hash{"BNLN_TR3_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                }

        #       Now, TRP
        if ($numbeads_hash{$complete_sequence[$i]} ==5 ){
                printf FITP ("%d\t%d\t%d\t%1.3f\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+2,1,$parameter_hash{"BNLN_TR1_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+3,1,$parameter_hash{"BNLN_TR2_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%s\t%s%d\n",${backbone_atoms[$i]}+2,${backbone_atoms[$i]}+3,1,$parameter_hash{"BNLN_TR3_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%s\t%s%d\n",${backbone_atoms[$i]}+2,${backbone_atoms[$i]}+4,1,$parameter_hash{"BNLN_TR3_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                printf FITP ("%d\t%d\t%d\t%1.3f\t%s\t%s%d\n",${backbone_atoms[$i]}+3,${backbone_atoms[$i]}+4,1,$parameter_hash{"BNLN_TR3_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                }
        }

}


        if ($parameter_hash{"ITV_BONDS"} eq "constraints"){
                printf FITP ("%s\n",";bc-sc constraints (ITV)");

for ( my $i=0; $i<$total_residues;$i++){
        if (($complete_sequence[$i] eq "I") || ($complete_sequence[$i] eq "V") || ($complete_sequence[$i] eq "T"))
        {
                printf FITP ("%d\t%d\t%d\t%1.3f\t%s\t%s%s\n",$backbone_atoms[$i],${backbone_atoms[$i]}+1,1,$parameter_hash{"BNLN_SC1_$aa_hash{${complete_sequence[$i]}}"},";",$aa_hash{${complete_sequence[$i]}},$i+1);
                }
        }
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#####################################BONDS DONE#################################

#####################################ANGLES#####################################
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

printf FITP ("\n%s\n","[angles]");

### First, let us write the back-bone angles.
### Backbone angles are ss dependent.
printf FITP ("%s\n",";backbone-backbone-backbone angles");

for (my $j=1;$j<=$number_of_chains; $j++){
        for ( my $i=$begin_chain[$j]+1; $i<$end_chain[$j];$i++) {
                printf FITP ("%d\t%d\t%d\t%d\t%1.2f\t%d\t%s\t%s\n",$backbone_atoms[$i-1],$backbone_atoms[$i],$backbone_atoms[$i+1],2,$parameter_hash{"ANGL_BBN_$ss_hash{$complete_ssdump[$i]}"},$parameter_hash{"ANKB_BBN_$ss_hash{${complete_ssdump[$i]}}"},";","$ss_hash{${complete_ssdump[$i-1]}}"."-"."$ss_hash{${complete_ssdump[$i]}}"."-"."$ss_hash{${complete_ssdump[$i+1]}}");

                }
        }

        printf FITP ("%s\n",";backbone-backbone-sidechain angles");

for (my $j=1;$j<=$number_of_chains; $j++){
        for ( my $i=$begin_chain[$j]; $i<=$end_chain[$j];$i++) {
                if($numbeads_hash{$complete_sequence[$i]} != 1 ){
                        if ( $i== $begin_chain[$j]){
                        printf FITP ("%d\t%d\t%d\t%d\t%1.2f\t%d\t%s\t%s%d\n",$backbone_atoms[$i+1],$backbone_atoms[$i],${backbone_atoms[$i]}+1,2,$parameter_hash{"ANGL_SC1_$aa_hash{$complete_sequence[$i]}"},$parameter_hash{"ANKB_SC1_$aa_hash{${complete_sequence[$i]}}"},";","$ss_hash{${complete_ssdump[$i+1]}}"."-"."$aa_hash{${complete_sequence[$i]}}",$i+1);
                        }
                        else{
                                printf FITP ("%d\t%d\t%d\t%d\t%1.2f\t%d\t%s\t%s%d\n",$backbone_atoms[$i-1],$backbone_atoms[$i],${backbone_atoms[$i]}+1,2,$parameter_hash{"ANGL_SC1_$aa_hash{$complete_sequence[$i]}"},$parameter_hash{"ANKB_SC1_$aa_hash{${complete_sequence[$i]}}"},";","$ss_hash{${complete_ssdump[$i-1]}}"."-"."$aa_hash{${complete_sequence[$i]}}",$i+1);
                                }
                        }

                }
        }

        printf FITP ("%s\n",";backbone-sidechain-sidechain angles : ARG and RINGS");

for (my $j=1;$j<=$number_of_chains; $j++){
        for ( my $i=$begin_chain[$j]; $i<=$end_chain[$j];$i++) {
                if($numbeads_hash{$complete_sequence[$i]} == 3 ){
                        printf FITP ("%d\t%d\t%d\t%d\t%1.2f\t%d\t%s\t%s%d\n",$backbone_atoms[$i],${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+2,2,$parameter_hash{"ANGL_SC2_$aa_hash{$complete_sequence[$i]}"},$parameter_hash{"ANKB_SC2_$aa_hash{${complete_sequence[$i]}}"},";","$aa_hash{${complete_sequence[$i]}}",$i+1);
                        }
                        elsif(($numbeads_hash{$complete_sequence[$i]} == 4 )||($numbeads_hash{$complete_sequence[$i]} == 5 )){
                        printf FITP ("%d\t%d\t%d\t%d\t%1.2f\t%d\t%s\t%s%d\n",$backbone_atoms[$i],${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+2,2,$parameter_hash{"ANGL_SC2_$aa_hash{$complete_sequence[$i]}"},$parameter_hash{"ANKB_SC2_$aa_hash{${complete_sequence[$i]}}"},";","$aa_hash{${complete_sequence[$i]}}",$i+1);
                        printf FITP ("%d\t%d\t%d\t%d\t%1.2f\t%d\t%s\t%s%d\n",$backbone_atoms[$i],${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+3,2,$parameter_hash{"ANGL_SC3_$aa_hash{$complete_sequence[$i]}"},$parameter_hash{"ANKB_SC3_$aa_hash{${complete_sequence[$i]}}"},";","$aa_hash{${complete_sequence[$i]}}",$i+1);
                        }

                }
        }



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#####################################ANGLES DONE###############################


#####################################DIHEDRALS##################################
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

printf FITP ("\n%s\n","[dihedrals]");

# First, impropers (GMX Type 2, for rings)
printf FITP ("%s\n",";improper dihedral angles");
for (my $i=0;$i< $total_residues; $i++){
        if ($numbeads_hash{$complete_sequence[$i]} >= 4 ){
                printf FITP ("%d\t%d\t%d\t%d\t%d\t%1.2f\t%d\t%s\t%s%d\n",$backbone_atoms[$i],${backbone_atoms[$i]}+2,${backbone_atoms[$i]}+3,${backbone_atoms[$i]}+1,2,$parameter_hash{"DIAN_RG1_$aa_hash{$complete_sequence[$i]}"},$parameter_hash{"DIKB_RG1_$aa_hash{${complete_sequence[$i]}}"},";","$aa_hash{${complete_sequence[$i]}}",$i+1);
                        }
                if ($numbeads_hash{$complete_sequence[$i]} == 5 ){
                printf FITP ("%d\t%d\t%d\t%d\t%d\t%1.2f\t%d\t%s\t%s%d\n",${backbone_atoms[$i]}+1,${backbone_atoms[$i]}+2,${backbone_atoms[$i]}+4,${backbone_atoms[$i]}+3,2,$parameter_hash{"DIAN_RG2_$aa_hash{$complete_sequence[$i]}"},$parameter_hash{"DIKB_RG2_$aa_hash{${complete_sequence[$i]}}"},";","$aa_hash{${complete_sequence[$i]}}",$i+1);
                        }
                }

printf FITP ("%s\n",";proper dihedral angles");
if ($parameter_hash{YES_BBN_PROPDIH} eq "yes") {
printf FITP ("%s\n",";helix backbone dihedrals");
for (my $j=1; $j<=$number_of_chains; $j++) {
	for (my $i=$begin_chain[$j]+1; $i<$end_chain[$j]-1; $i++) {
		if ( ($complete_ssdump[$i-1] eq "H") && ($complete_ssdump[$i] eq "H") 
				&& ($complete_ssdump[$i+1] eq "H") 
				&& ($complete_ssdump[$i+2] eq "H") ) {
			# if one of the middle residues is a proline
			if ( ($complete_sequence[$i] eq "P") 
					|| ($complete_sequence[$i+1] eq "P") ) {
				printf FITP ("%d\t%d\t%d\t%d\t%d\t%1.2f\t%d\t%d\t%s\t%s%d\n", 
					${backbone_atoms[$i-1]}, $backbone_atoms[$i], 
					${backbone_atoms[$i+1]}, ${backbone_atoms[$i+2]}, 1, 
					$parameter_hash{"DIAN_BBN_$ss_hash{$complete_ssdump[$i]}"}, 
					$parameter_hash{"DIKB_HLX_PRO"}, 1, ";", 
					"$aa_hash{${complete_sequence[$i]}}", $i+1);
			} else {
				printf FITP ("%d\t%d\t%d\t%d\t%d\t%1.2f\t%d\t%d\t%s\t%s%d\n", 
					${backbone_atoms[$i-1]}, $backbone_atoms[$i], 
					${backbone_atoms[$i+1]}, ${backbone_atoms[$i+2]}, 1, 
					$parameter_hash{"DIAN_BBN_$ss_hash{$complete_ssdump[$i]}"}, 
					$parameter_hash{"DIKB_BBN_$ss_hash{${complete_ssdump[$i]}}"}, 
					1, ";", "$aa_hash{${complete_sequence[$i]}}", $i+1);
			}
		}
		# Dihedrals instead of elastic networks?
		if ( !($ELASTIC_NETWORK) 
				&& ($complete_ssdump[$i-1] eq "E") 
				&& ($complete_ssdump[$i] eq "E") 
				&& ($complete_ssdump[$i+1] eq "E") 
				&& ($complete_ssdump[$i+2] eq "E") ) {
			printf FITP ("%d\t%d\t%d\t%d\t%d\t%1.2f\t%d\t%d\t%s\t%s%d\n", 
				${backbone_atoms[$i-1]}, $backbone_atoms[$i], 
				${backbone_atoms[$i+1]}, ${backbone_atoms[$i+2]}, 1, 
				$parameter_hash{"DIAN_BBN_$ss_hash{$complete_ssdump[$i]}"}, 
				$parameter_hash{"DIKB_BBN_$ss_hash{${complete_ssdump[$i]}}"}, 
				1, ";", "$aa_hash{${complete_sequence[$i]}}", $i+1);
		}
	}
}

}

print"\nDONE. Your Topology has been Generated\n";
print"\nPlease check your itp file for bugs/errors...\n";
##########################
# OBLIGATORY FUNNY LINE
print "\n(gcq#007): \"TOPOLOGY FILE FOR MAKE BENEFIT GLORIOUS CG SIMULATION\" \n\n";
#########################
