#!/usr/bin/awk -f
#Converts Atomic PDB file to Coarse-grain PDB file
#The atom types are default values - BN0, SC1, SC2 
#Should be corrected by using correct itp file
#if columns are not separate (e.g. x or y or z > 99) the script doesnt work
#in that case insert an empty column or use PDBCat to separate columns
#d! groningen   29.01.08

{if ($1=="HEADER"||$1=="REMARK"||$1=="CRYST1"||$1=="MODEL"||$1=="TER"||$1=="ENDMDL"||$1=="END")
print $0
if($1=="ATOM" && $4=="ARG" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if($1=="ATOM" && $4=="ARG" && $3=="CG")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if($1=="ATOM" && $4=="ARG" && $3=="NE")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC2", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="ALA" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="ASN" && $3=="CA" )
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="ASN" && $3=="CG")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="ASP" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="ASP" && $3=="CG")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="CYS" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="CYS" && $3=="SG")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="GLN" && $3=="CA" )
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="GLN" && $3=="CB")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="GLU" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="GLU" && $3=="CB")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="GLY" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="HIS" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="HIS" && $3=="CB")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="HIS" && $3=="ND1")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC2", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="HIS" && $3=="NE2")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC3", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="ILE" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="ILE" && $3=="CD")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="LEU" && $3=="CA" )
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="LEU" && $3=="CG")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="LYS" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="LYS" && $3=="CG")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="LYS" && $3=="CE")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC2", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="MET" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="MET" && $3=="CG")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="PHE" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="PHE" && $3=="CG")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="PHE" && $3=="CE1")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC2", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="PHE" && $3=="CE2")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC3", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="PRO" && $3=="CA") 
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="PRO" && $3=="CG")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="SER" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="SER" && $3=="CB")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="THR" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="THR" && $3=="CB")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="TRP" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="TRP" && $3=="CD1")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="TRP" && $3=="CD2")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC2", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="TRP" && $3=="CE2")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC3", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="TRP" && $3=="CH2")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC4", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="TYR" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="TYR" && $3=="CG") 
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="TYR" && $3=="CE1")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC2", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="TYR" && $3=="CE2")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC3", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="VAL" && $3=="CA")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "BN0", $4, $6, $7, $8, $9,$10,$11);
else if ($1=="ATOM" && $4=="VAL" && $3=="CB")
printf("%4s  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",$1, $2, "SC1", $4, $6, $7, $8, $9,$10,$11);
}
