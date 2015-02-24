#!/usr/sbin/perl -w  
#bond - 07/05
#-lots of additional crap
#
#05/06 - new distances, angles (thanks timbo) + triangles for aromatics + harmonic bonds instead of distance restraints
#new force and angle constants
#impoper dihedrals for maintaining directionality of sidechains wrt secondary structure
#tryptophan corrected
#things tided up for better output
#
#06/06 - did sidechain distances (and angleS) properly - for most residues averaged over atoms, then multiplied by ratio of atoms/4
#      - exceptions were Lys, Arg, and rings
#hall - 07/08
#Sidekick specific edits
#

#This line is hard coded, unfortunately.
#$gromacs = "/Users/Shared/Sidekick/gromacs/bin/";
#$ENV{'PATH'} .= ':/Users/Shared/Sidekick/gromacs/bin/:/usr/local/Sidekick/gromacs/bin/:/nfsmount/Sidekick/gromacs/bin/';
$ENV{'PATH'} .= ":" . $ARGV[$0];
$inpdb="atomistic.pdb" ; #input pdb of protein - no ions,cofactors,water, one chain

system ("pdb2gmx -f $inpdb -ignh <<EOD >&log.txt
2
EOD");


system ("editconf -f conf.gro -o inx.pdb -d 1 >&log.txt");



system ("echo 'cpp                      = /usr/bin/cpp\nintegrator               = steep\nnsteps                   = 100' > temp.mdp");
system ("grompp -f temp.mdp -c inx.pdb -p topol.top -r 0.3 -o temp.tpr >&log.txt");

system ("g_hbond -f inx.pdb -s temp.tpr  -hbn hbond.ndx  <<EOD >&log.txt
7
7
EOD");

##########COMMENTED OUT AS HANDLED BY SIDEKICK###########
#/*
#
#system ("do_dssp  -f inx.pdb -s temp.tpr -o ss.xpm << EOD >&log.txt
#1
#EOD");
#
#system ("cat ss.xpm | grep -v '*' | grep -v '1' > structure.txt");
#*/

@structure=();
open (INSTRUCTURE, "structure.txt") || die "nope\n";
while (<INSTRUCTURE>) 	{
			chomp;
			push (@structure,$_);
			}
close(INSTRUCTURE);

@structure = reverse (@structure);


###################################################### process the index group of hydrogen-bonds

@hbondsa=();
@hbondsb=();
$now=0;
open(INNDX, "hbond.ndx") || die "nope\n";
while (<INNDX>) {
		chomp;
		($xoa, $xob, $xoc) = split;
		if ($now == 1) {push(@hbondsa, $xoa);push(@hbondsb, $xoc)}
		if ($xob eq "hbonds_MainChain+H") {$now=1}
		}
close(INNDX);




@hbondordera=();
@hbondorderb=();
@ordera=();
@orderb=();		
open(INPDB, "inx.pdb") || die "nope\n";
while (<INPDB>) {
		chomp;
		($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
		if ($woa eq "ATOM")
			{
			$loop=0;
			while ($loop<@hbondsb)
				{
				if ($wob == $hbondsb[$loop])
					{
					push(@hbondorderb, $woe);
					push(@orderb, $loop);
					}
				if ($wob == $hbondsa[$loop])
					{
					push(@hbondordera, $woe);
					push(@ordera, $loop);
					}

				$loop++
				}
			}
		}
close(INPDB);



@hbondlista=();
@hbondlistb=();
$loop=0;
while ($loop<@orderb)
	{
	$loopinloop=0;
	while ($loopinloop<@orderb)
		{
		if ($orderb[$loopinloop]==$loop)
			{
			$x=$hbondorderb[$loopinloop];
			push(@hbondlistb, $x)
			}
		$loopinloop++
		}
	
	$loop++
	}

$loop=0;
while ($loop<@ordera)
	{
	$loopinloop=0;
	while ($loopinloop<@ordera)
		{
		if ($ordera[$loopinloop]==$loop)
			{
			$x=$hbondordera[$loopinloop];
			push(@hbondlista, $x)
			}
		$loopinloop++
		}
	
	$loop++
	}


#####make the hbond list array

@hbondingdonor=();
@hbondingacceptor=();
@hbondingdonoracceptor=();

$loopa=0;
$check="one";
while ($loopa<@hbondlista)
		{
		$loopb=0;
		while ($loopb<@hbondlistb)
			{
			if ($hbondlistb[$loopb] == $hbondlista[$loopa]) {$check="both"}; $loopb++
			}
		if ($check eq "both") {push(@hbondingdonoracceptor, $hbondlista[$loopa])}
		if ($check eq "one") {push(@hbondingdonor, $hbondlista[$loopa])}
		$check="one";
		$loopa++
		}

$loopa=0;
$check="one";
while ($loopa<@hbondlistb)
		{
		$loopb=0;
		while ($loopb<@hbondlista)
			{
			if ($hbondlistb[$loopa] == $hbondlista[$loopb]) {$check="both"}; $loopb++
			}
		if ($check eq "one") {push(@hbondingacceptor, $hbondlistb[$loopa])}
		$check="one";
		$loopa++
		}



############################all atoms which take part in hydrogen bonds =

@hbondlists=@hbondlista;
push(@hbondlists, @hbondlistb);
push(@hbondlist, $hbondlists[0]);
$loopa=0;
$check="new";
while ($loopa<@hbondlists)
		{
		$loopb=0;
		while ($loopb<@hbondlist)
			{
			if ($hbondlist[$loopb] == $hbondlists[$loopa]) {$check="nonew"}; $loopb++
			}
		if ($check eq "new") {push(@hbondlist, $hbondlists[$loopa])}
		$check="new";
		$loopa++
		}

			
			
######################### @hbondlist = list of residues that take part in hbonds
######################### @hbondlista and @hbondlistb = lists of pairs of h-bonded residues



##############################################now process pdb to write out body of itp with atoms etc.


open(INPDB, "inx.pdb") || die "nope\n";
open(OUTPDB, ">protein-cg.pdb") || die "nope\n";
open(OUTITP, ">atoms.list") || die "nope\n";


$resnum = 1;
$atomnum = 1;

while (<INPDB>) {
		chomp;
		($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
		($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
		if ($woa eq "ATOM")
			{
			
			if (($wod eq ALA) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
	                if (($wod eq ALA) && ($woc eq CB)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }###reinserted Ala CB
			if (($wod eq CYS) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq CYS) && ($woc eq SG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq ASP) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq ASP) && ($woc eq CG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq GLU) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq GLU) && ($woc eq CD)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq PHE) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq PHE) && ($woc eq CE1)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq PHE) && ($woc eq CE2)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq GLY) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq HIS) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq HIS) && ($woc eq CD2)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq HIS) && ($woc eq ND1)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq ILE) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq ILE) && ($woc eq CG1)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq LYS) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq LYS) && ($woc eq CG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq LYS) && ($woc eq NZ)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq LEU) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq LEU) && ($woc eq CG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq MET) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq MET) && ($woc eq CG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq ASN) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq ASN) && ($woc eq CG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq PRO) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq PRO) && ($woc eq CG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq GLN) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq GLN) && ($woc eq CG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq ARG) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq ARG) && ($woc eq CG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq ARG) && ($woc eq CZ)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq SER) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq SER) && ($woc eq OG)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq THR) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq THR) && ($woc eq OG1)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq VAL) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq VAL) && ($woc eq CG1)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq TRP) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq TRP) && ($woc eq CZ3)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq TRP) && ($woc eq NE1)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq TYR) && ($woc eq CA)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq TYR) && ($woc eq CD1)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			if (($wod eq TYR) && ($woc eq CE2)) { print OUTPDB $spa,$woa,$spb,$wob,$spc,$woc,$spd,$wod,$spe,$woe,$spf,$wof,$spg,$wog,$sph,$woh,$spi,$woi,$spj,$woj,"\n" }
			
			
				if ($resnum == $woe) 
					{

					$loop=0;
					$donoracceptor="no";
					while ($loop<@hbondingdonoracceptor) 
						{
						if ($resnum == $hbondingdonoracceptor[$loop]) { $donoracceptor="yes" }
						$loop++
						}
						
					$loop=0;
					$donoronly="no";
					while ($loop<@hbondingdonor) 
						{
						if ($resnum == $hbondingdonor[$loop]) { $donoronly="yes" }
						$loop++
						}
						
					$loop=0;
					$acceptoronly="no";
					while ($loop<@hbondingacceptor) 
						{
						if ($resnum == $hbondingacceptor[$loop]) { $acceptoronly="yes" }
						$loop++
						}
						
					$process=0;
					if ($donoracceptor eq "yes") { print OUTITP $atomnum,"  ","Nm","  ",$resnum,"  ",$wod,"  ","CA","  ",$atomnum,"  ","0.0","\n";$process=1}
					if ($donoronly eq "yes") { print OUTITP $atomnum,"  ","Nda","  ",$resnum,"  ",$wod,"  ","CA","  ",$atomnum,"  ","0.0","\n";$process=1}
					if ($acceptoronly eq "yes") { print OUTITP $atomnum,"  ","Nda","  ",$resnum,"  ",$wod,"  ","CA","  ",$atomnum,"  ","0.0","\n";$process=1}
					if ($process==0) {print OUTITP $atomnum,"  ","P","  ",$resnum,"  ",$wod,"  ","CA","  ",$atomnum,"  ","0.0","\n"}
					
					$atomnum++;
					
			                if ($wod eq ALA) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++} ###reinserted Ala CB
					if ($wod eq CYS) { print OUTITP $atomnum,"  ","N0 ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq ASP) { print OUTITP $atomnum,"  ","Qa ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq GLU) { print OUTITP $atomnum,"  ","Qa ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq PHE) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq PHE) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CG","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq HIS) { print OUTITP $atomnum,"  ","Nd ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq HIS) { print OUTITP $atomnum,"  ","Nda","  ",$resnum,"  ",$wod,"  ","CG","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq ILE) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq LYS) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq LYS) { print OUTITP $atomnum,"  ","Qd ","  ",$resnum,"  ",$wod,"  ","CG","  ",$atomnum,"  ","1.0","\n"; $atomnum++}
					if ($wod eq LEU) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq MET) { print OUTITP $atomnum,"  ","N0 ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq ASN) { print OUTITP $atomnum,"  ","P  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq PRO) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq GLN) { print OUTITP $atomnum,"  ","P  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq ARG) { print OUTITP $atomnum,"  ","Nd ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq ARG) { print OUTITP $atomnum,"  ","Qd ","  ",$resnum,"  ",$wod,"  ","CG","  ",$atomnum,"  ","1.0","\n"; $atomnum++}
					if ($wod eq SER) { print OUTITP $atomnum,"  ","Nda","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq THR) { print OUTITP $atomnum,"  ","Nda","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq VAL) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq TRP) { print OUTITP $atomnum,"  ","Nd ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq TRP) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CG","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq TYR) { print OUTITP $atomnum,"  ","C  ","  ",$resnum,"  ",$wod,"  ","CB","  ",$atomnum,"  ","0.0","\n"; $atomnum++}
					if ($wod eq TYR) { print OUTITP $atomnum,"  ","P  ","  ",$resnum,"  ",$wod,"  ","CG","  ",$atomnum,"  ","0.0","\n"; $atomnum++} ### Changed P to Nd
					
					$resnum++
					}
					
				#print OUTPDB ,"\n";
				
			
			}


			
		
						
	}

close(INPDB);
close(OUTPDB);
close(OUTITP);




open(INATOMS, "atoms.list") || die "nope\n";
open(OUTBONDS, ">bonds.list") || die "nope\n";
open(OUTANGLES, ">angles.list") || die "nope\n";
open(OUTDIHEDRALS, ">dihedrals.list") || die "nope\n";

@sideanglesa=(); ## sidechain angles for Lys, Arg (180)
@sideanglesb=(); ## sidechain angles for Trp
@sideanglesc=(); ## sidechain angles for Phe
@sideanglesd=(); ## sidechain angles for His
@sideanglese=(); ## sidechain angles for Tyr
@hbondordera=();
@hbondorderb=();
@ordera=();
@orderb=();	
@hbondcga=();
@hbondcgb=();
@angles = ();
$resnum = 1;
$atomnum = 1;
#$bondloop = 0;
#$force="weak";
$type="none";
while (<INATOMS>) {
		chomp;
		($yoa, $yob, $yoc, $yod, $yoe, $yof, $yog) = split;
		
			
			if  (($yod eq "LYS") ||($yod eq "ARG"))
				{
				push (@sideanglesa,$yoa)
				}
			if  ($yod eq "TRP")
				{
				push (@sideanglesb,$yoa)
				}
			if  ($yod eq "PHE")
				{
				push (@sideanglesc,$yoa)
				}
			if  ($yod eq "HIS")
				{
				push (@sideanglesd,$yoa)
				}
			if  ($yod eq "TYR")
				{
				push (@sideanglese,$yoa)
				}	
		
			if ($resnum == 1)
				{
				$ca=$yoa;
				push (@angles,$yoa);
				$resnum++;
				
				
						
		
				$loop=0;
	 			while ($loop<@hbondlista)
 					{
 					if ($hbondlista[$loop]==$yoc) {push (@hbondordera,$yoa);push(@ordera, $loop);} #print $yoa," ",$loop,"\n"}
					if ($hbondlistb[$loop]==$yoc) {push (@hbondorderb,$yoa);push(@orderb, $loop);} #print $yoa," ",$loop,"\n"}
					$loop++
					}
				}
			
			if ($resnum == $yoc) 
				{
		
				$loop=0;
	 			while ($loop<@hbondlista)
 					{
					if ($hbondlista[$loop]==$yoc) {push (@hbondordera,$yoa);push(@ordera, $loop);} #print $yoa," ",$loop,"\n"}
					if ($hbondlistb[$loop]==$yoc) {push (@hbondorderb,$yoa);push(@orderb, $loop);} #print $yoa," ",$loop,"\n"}
					$loop++
 					}

				print OUTBONDS $ca,"  ",$yoa,"  1  0.36  5000\n";
				$ca=$yoa;
				push (@angles,$yoa);
				
				if ($resnum>2)
					{
					$resi=$resnum-1;
					$resii=$resnum-2;
					$resiii=$resnum-3;
				#	$resiv=$resnum-4;			
				#	$loopangle=0;
				#	while ($loopangle<@hbondlist)
				#		{
				#		if ($resi==$hbondlist[$loopangle]) {$particles++}
				#		if ($resii==$hbondlist[$loopangle]) {$particles++}	
				#		if ($resiii==$hbondlist[$loopangle]) {$particles++}
				#		if ($particles==3) {$force="strong"}
				#		print "particles=",$particles,"\n";	
				#		$loopangle++
				#		}	

					if (($structure[$resiii] eq '"H",') && ($structure[$resii] eq '"H",') && ($structure[$resi] eq '"H",')) {$type="helix"}
					if (($structure[$resiii] eq '"E",') && ($structure[$resii] eq '"E",') && ($structure[$resi] eq '"E",')) {$type="sheet"}

					if ($type eq "helix") {print OUTANGLES $angles[$resiii],"  ",$angles[$resii],"  ",$angles[$resi],"  2    90  350 ; angle ",$resii,"-",$resi,"-",$resnum,"\n"}
					if ($type eq "sheet") {print OUTANGLES $angles[$resiii],"  ",$angles[$resii],"  ",$angles[$resi],"  2    130 350 ; angle ",$resii,"-",$resi,"-",$resnum,"\n"}
					if ($type eq "none") {print OUTANGLES $angles[$resiii],"  ",$angles[$resii],"  ",$angles[$resi],"  2    120 250  ; angle ",$resii,"-",$resi,"-",$resnum,"\n"}
					
					#18/05/06 - improper dihedrals
					
					if ($prev ne "GLY") # &&  ($prev ne "PHE") &&  ($prev ne "HIS") &&  ($prev ne "TYR") &&  ($prev ne "TRP"))
						{
						$dihside=$angles[$resii]+1;
						print OUTDIHEDRALS $angles[$resii],"  ",$angles[$resiii],"  ",$dihside,"  ",$angles[$resi],"  2    0   0.05  ; improper angle centred on residue ",$resi,"\n"
						}			
					
					
					
					$type="none";
					}
				$resnum++;
				}
				
			$cai=$ca+1;
			$caii=$ca+2;
			if (($yod eq "ALA") && ($yoa==$cai))               ############For the Ala CB particle
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.09  5000\n"
				}
			if (($yod eq "CYS") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.17  5000\n"
				}
			if (($yod eq "ASP") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.26  5000\n"
				}
			if (($yod eq "ILE") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.26  5000\n"
				}
			if (($yod eq "LEU") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.28  5000\n"
				}
			if (($yod eq "ASN") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.26  5000\n"
				}
			if (($yod eq "GLU") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.36  5000\n"
				}
			if (($yod eq "PHE") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.42  5000\n",$cai,"  ",$caii,"  1  0.42  5000\n",$ca,"  ",$caii,"  1  0.42  5000\n"
				}
			if (($yod eq "HIS") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.38  5000\n",$cai,"  ",$caii,"  1  0.38  5000\n",$ca,"  ",$caii,"  1  0.38  5000\n"
				}
			if (($yod eq "LYS") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.31  5000\n",$cai,"  ",$caii,"  1  0.21  5000\n"
				}
			if (($yod eq "MET") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.31  5000\n"
				}
			if (($yod eq "GLN") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.36  5000\n"
				}
			if (($yod eq "PRO") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.19  5000\n"
				}
			if (($yod eq "ARG") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.37  5000\n",$cai,"  ",$caii,"  1  0.25  5000\n"
				}
			if (($yod eq "SER") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.16  5000\n"
				}
			if (($yod eq "THR") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.2  5000\n"
				}
			if (($yod eq "VAL") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.2  5000\n"
				}
			if (($yod eq "TRP") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.47  5000\n",$cai,"  ",$caii,"  1  0.28  5000\n",$ca,"  ",$caii,"  1  0.38  5000\n"
				}
			if (($yod eq "TYR") && ($yoa==$cai))
				{
				print OUTBONDS $ca,"  ",$cai,"  1  0.38  5000\n",$cai,"  ",$caii,"  1  0.38  5000\n",$ca,"  ",$caii,"  1  0.38  5000\n"
				}
			
		$prev=$yod;
		}
					
###arg,lys sidechain angles					
$x=0;
while ($x<@sideanglesa)
	{
	$xx=$x+1;
	$xxx=$x+2;
	print OUTANGLES $sideanglesa[$x],"  ",$sideanglesa[$xx],"  ",$sideanglesa[$xxx],"  2    180  250 ; angle for Arg and Lys sidechains \n";
	$x=$x+3
	}

#tryptophan sidechain angles
$x=0;
while ($x<@sideanglesb)
	{
	$xx=$x+1;
	$xxx=$x+2;
	print OUTANGLES $sideanglesb[$x],"  ",$sideanglesb[$xx],"  ",$sideanglesb[$xxx],"  2    54  250 ; angle for Trp sidechains \n";
	$x=$x+3
	}
		
#phe sidechain angles
$x=0;
while ($x<@sideanglesc)
	{
	$xx=$x+1;
	$xxx=$x+2;
	print OUTANGLES $sideanglesc[$x],"  ",$sideanglesc[$xx],"  ",$sideanglesc[$xxx],"  2    60  250 ; angle for Phe sidechains \n";
	$x=$x+3
	}
	
#his sidechain angles
$x=0;
while ($x<@sideanglesd)
	{
	$xx=$x+1;
	$xxx=$x+2;
	print OUTANGLES $sideanglesd[$x],"  ",$sideanglesd[$xx],"  ",$sideanglesd[$xxx],"  2    60  250 ; angle for His sidechains \n";
	$x=$x+3
	}
		
#tyr sidechain angles
$x=0;
while ($x<@sideanglese)
	{
	$xx=$x+1;
	$xxx=$x+2;
	print OUTANGLES $sideanglese[$x],"  ",$sideanglese[$xx],"  ",$sideanglese[$xxx],"  2    60  250 ; angle for Tyr sidechains \n";
	$x=$x+3
	}
		
		







close(INATOMS);
close(OUTBONDS);
close(OUTANGLES);




#open(INATOMS, "atoms.list") || die "nope\n";
#open(OUTBONDS, ">bonds.list") || die "nope\n";
#open(OUTANGLES, ">angles.list") || die "nope\n";
#open(OUTDIHEDRALS, ">dihedrals.list") || die "nope\n";

#while (<INATOMS>) {
	#	chomp;
		#($yoa, $yob, $yoc, $yod, $yoe, $yof, $yog) = split;
		


#if ($resnum == 1)
	#			{
		#		$ca=$yoa;
			#	push (@angles,$yoa);
				#$resnum++;
				
				#if (($yoe eq "CA") && ($yod ne "GLY")
				#{
				#print OUTDIHEDRALS $yoe,"   2    0   0.05 \n";
				#}
			#}








#close(INATOMS);
#close(OUTBONDS);
#close(OUTANGLES);
#close(OUTDIHEDRALS);


	
open(OUTHBONDS, ">temphbonds") || die "nope\n";




@hbondcga=();
@hbondcgb=();
$loop=0;
while ($loop<@orderb)
	{
	$loopinloop=0;
	while ($loopinloop<@orderb)
		{
		if ($orderb[$loopinloop]==$loop)
			{
			$x=$hbondorderb[$loopinloop];
			push(@hbondcgb, $x)
			}
		$loopinloop++
		}
	
	$loop++
	}

$loop=0;
while ($loop<@ordera)
	{
	$loopinloop=0;
	while ($loopinloop<@ordera)
		{
		if ($ordera[$loopinloop]==$loop)
			{
			$x=$hbondordera[$loopinloop];
			push(@hbondcga, $x)
			}
		$loopinloop++
		}
	
	$loop++
	}









$loop=0;
while ($loop<@hbondcga)
	{
	print OUTHBONDS $hbondcga[$loop],"  ",$hbondcgb[$loop],"  1  0.60  1000 ; hydrogen bonds ",$hbondlista[$loop],"-",$hbondlistb[$loop],"\n";
	$loop++
 	}

			
	
close(OUTHBONDS);

$yom=0;$yoh=0;$yoi=0;$yoj=0;$yok=0;$yol=0;

$already="no";
@checkhbonds=();					
open(INHBONDS, "temphbonds") || die "nope\n";
open(OUTHBONDS, ">hbonds.list") || die "nope\n";
while (<INHBONDS>) 	{
			chomp;
			($yoa, $yob, $yoc, $yod, $yoe, $yof, $yog, $yoh, $yoi, $yoj, $yok, $yol, $yom) = split;
			
			$loop=0;
			while ($loop<@checkhbonds)
				{
				$next=$loop+1;
				if (($checkhbonds[$loop]==$yoa) && ($checkhbonds[$next]==$yob))
					{
					$already="yes"
					}
				$loop=$loop+2
				}
			
					
			if ($yoa != $yob)
				{
				if ($already eq "no")
					{
						print OUTHBONDS $_,"\n";
						push (@checkhbonds, $yob, $yoa);
					}
				$already="no"
				}
			}
close(INHBONDS);
close(OUTHBONDS);
			
	


		
system ("echo '[ moleculetype ]\n; molname       nrexcl\nPROTEIN            1\n' > tempa");
system ("echo '\n[ atoms ]\n;id type resnr residu atom cgnr   charge' > tempb");
system ("echo '\n[ bonds ]' > tempc"); 
#system ("echo '\n[ distance_restraints ]\n;       atom1  atom2  type  index  type2  low   up1   up2   fac' > tempd");
system ("/usr/bin/touch tempd");
system ("echo '\n[ angles ]\n;  i     j    k     funct   angle  force.c.' > tempe");
system ("echo '\n[ dihedrals ]\n;  j  i  k  l   funct   angle  force.c.' > tempf");


#system ("cat tempa tempb atoms.list tempc bonds.list tempd hbonds.list tempe angles.list tempf dihedrals.list> protein-cg.itp");

system ("/bin/cat tempa tempb atoms.list tempc bonds.list tempd hbonds.list tempe angles.list > protein-cg.itp");


$number=1;
open(INAT, "atoms.list") || die "nope\n";
open(OUTPOS, ">posre.itp") || die "nope\n";
print OUTPOS "[ position_restraints ]\n; atom  type      fx      fy      fz\n";
while (<INAT>)	{
           	print OUTPOS "    ",$number,"     1  1000  1000  1000\n";
		$number++
		}
close(INAT);
close(OUTPOS);

print "Successful Finish\n"
 
##include "ff_v1.4_x.itp"
##include "protein-cg.itp"
# 
#[ system ]
#; name
#PROTEIN CG
# 
#[ molecules ]
#; name  number
#PROTEIN 1
