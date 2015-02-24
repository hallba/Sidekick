#!/usr/bin/perl

use warnings;
use strict;


unless ( $ARGV[1] )
{
	die "Usage:analyze.pl <file with list of pdbfiles> <name string>\n";
}

#Read the list of directories from a file
my $dirfile = $ARGV[0];
my $namestring = $ARGV[1];
my @dirs;
open (DIR, "$dirfile" ) || die "Could not open $dirfile\n";
while ( my $line = <DIR>)
{
	chomp ($line);
	my @temp = split ( /\s+/, $line );
	push ( @dirs, $temp[0] );
}
close (DIR);

my $counter = -1;
#Go through each directory and run the analysis
foreach my $dir (@dirs)
{
	$counter++;
	open (FILE, "$dir");
	my $outfile = $namestring."_".$counter.".pdb";
	open (OUT, ">$outfile");
	while (my $line = <FILE> )
	{
		unless ($line =~ m/POP/ )
		{
			printf (OUT "$line");
			next;
		}
		$line =~ s/POP/$namestring/;
		printf (OUT "$line");
	}
	close (FILE);
	close (OUT);
}
