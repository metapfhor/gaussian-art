#!/usr/local/bin/perl
use Cwd;
my $file_name = $ARGV[0];
my $cwd = getcwd();
my $seperator=0;
my @forces;
open(FILE,"<$cwd/$file_name") or die "Cannot open ./$file_name to read: $!\n";
while ((my $line = <FILE>) && $seperator<2) {
	if($line =~ /Forces/){
		#We have entered the force section fo the output
		while(($line = <FILE>) && $seperator<2){			
			if($line =~ /-{2,}/){
				$seperator++;	
			}
			elsif($line =~ /(\d+)( +)(\d+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
				push(@forces,"$5$6","$8$9","$11$12");
			}
		}
		
	}
}
my $i=0;
while(my $force=$forces[$i]){
	print "$force\n";
	$i++;
}
close(FILE);
return @forces;

