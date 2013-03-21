#!/usr/local/bin/perl
use Cwd;

my $file_name = $ARGV[0];
my $cwd = getcwd();
my $seperator=0;
my @forces;
my $status=system("g09 temp.com");
my $fp="%16.10e";
my $fn="%16.9e";
my $force;

open(INPUT,"<$cwd/temp.log") or die "Cannot open $cwd/temp.log to read: $!\n";
open(OUTPUT,">$cwd/temp.forces")  or die "Cannot open $cwd/temp.forces to read: $!\n";
while ((my $line = <INPUT>) && $seperator<2) {
if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
	 if($1>0){
	 	print OUTPUT sprintf($fp,"$1")."\n";
	 }else{
	 	print OUTPUT sprintf($fn,"$1")."\n";
	 }
  } 
	if($line =~ /Forces/){
		#We have entered the force section of the output
		while(($line = <INPUT>) && $seperator<2){			
			if($line =~ /-{2,}/){
				$seperator++;	
			}
			elsif($line =~ /(\d+)( +)(\d+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
				if("$5$6">0){
					$force=sprintf($fp,"$5$6");
				}else{
					$force=sprintf($fn,"$5$6");
				}
				$force=$force."   ";
				if("$8$9">0){
					$force=$force.sprintf($fp,"$8$9");
				}else{
					$force=$force.sprintf($fn,"$8$9");
				}
				$force=$force."   ";
				if("$11$12">0){
					$force=$force.sprintf($fp,"$11$12");
				}else{
					$force=$force.sprintf($fn,"$11$12");
				}
				print OUTPUT "$force\n";
			}
		}
	}

}
close(INPUT);
close(OUTPUT);


