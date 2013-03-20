#!/usr/local/bin/perl
use Cwd;
my $file_name = $ARGV[0];
my $cwd = getcwd();
my $seperator=0;
my @forces;
#my $status=system("g09 temp.com");

open(INPUT,"<$cwd/temp.log") or die "Cannot open $cwd/temp.log to read: $!\n";
open(OUTPUT,">$cwd/temp.forces")  or die "Cannot open $cwd/temp.forces to read: $!\n";
while ((my $line = <INPUT>) && $seperator<2) {
#  print "$line";
	if($line =~ /Forces/){
		#We have entered the force section fo the output
		while(($line = <INPUT>) && $seperator<2){			
			if($line =~ /-{2,}/){
				$seperator++;	
			}
			elsif($line =~ /(\d+)( +)(\d+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
				print OUTPUT "$5$6   $8$9   $11$12\n";
			}
		}
		
	}
}
close(INPUT);
close(OUTPUT);


