#!/usr/bin/perl

for (my $i = 0; $i < 51; $i++)
{
	my $j = $i;
	my $file = "../pope/pope_".$i.".pdb";
	my $out =  "popc_".$j.".pdb";

	open (FILE, $file );
	open (OUT, ">$out");
	while ( my $line = <FILE>)
	{
		unless ($line =~ m/^ATOM/ )	
		{
			printf(OUT "$line");
			next;
		}
		$line =~ s/NH3/NC3/;
		$line =~ s/POE/POP/;
		printf(OUT "$line");
	}
	close (FILE);
	close (OUT);
	
}

