#!/usr/local/bin/perl
use Cwd;

my $file_name = $ARGV[0];
my $cwd = getcwd();
my $seperator=0;

my $status=system("g09 temp.com");
my $fp="%16.10e";
my $fn="%16.9e";
my $force;
my $forces;
my $enegry;


open(INPUT,"<$cwd/temp.log") or die "Cannot open $cwd/temp.log to read: $!\n";
open(OUTPUT,">$cwd/temp.forces")  or die "Cannot open $cwd/temp.forces to read: $!\n";
while ((my $line = <INPUT>) && $seperator<2) {
if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
	$energy=$1*27.211383;
	 if($energy>0){
	 	print OUTPUT sprintf($fp,"$energy")."\n";
	 }else{
	 	print OUTPUT sprintf($fn,"$energy")."\n";
	 }
  } 
	if($line =~ /Forces/){
		#We have entered the force section of the output
		while(($line = <INPUT>) && $seperator<2){			
			if($line =~ /-{2,}/){
				$seperator++;	
			}
			elsif($line =~ /(\d+)( +)(\d+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
				$force="$5$6"*1;
				if($force>0){
					$forces=sprintf($fp,"$force");
				}else{
					$forces=sprintf($fn,"$force");
				}
				$forces=$forces."   ";
				$force="$8$9"*1;
				if($force>0){
					$forces=$forces.sprintf($fp,"$force");
				}else{
					$forces=$forces.sprintf($fn,"$force");
				}
				$forces=$forces."   ";
				$force="$11$12"*1;
				if("$force">0){
					$forces=$forces.sprintf($fp,"$force");
				}else{
					$forces=$forces.sprintf($fn,"$force");
				}
				print OUTPUT "$forces\n";
			}
		}
	}

}
close(INPUT);
close(OUTPUT);


