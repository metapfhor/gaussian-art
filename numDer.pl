#!/usr/local/bin/perl
#Calculates the force on the first atom in the system by numerical differentiation
#for checking the units of gaussian and seeing where the numerical differentiation breaks down
use Cwd;

my $file_name = $ARGV[0];
my $cwd = getcwd();
my $seperator=0;
my $dx=0.1;
my $bohr=0.529;

#my $status=system("g09 ref.com");
my $energy;
my $force;
my $forces;
my $refEnegry;
my $refForces;
my $energy;
my $factor=1; 
my $delta=1;
my $coord;

open(INPUT,"<$cwd/ref.log") or die "Cannot open $cwd/ref.log to read: $!\n";
while (my $line = <INPUT>) {
#print $line;
	if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
		$refEnergy=$1;
		#print "$refEnergy\n";
	} 
	if($delta && $line =~ /(\d+)( +)(\d+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
				$forces="$5$6"*$factor."   "."$8$9"*$factor."   "."$11$12"*$factor;				
				$refForces=$forces;
				print  "Gaussian force (Hartree/bohr): $refForces\n";
				$delta=0;
			}
}
close(INPUT);

$delta=1;

open(INPUT,"<$cwd/ref.com") or die "Cannot open $cwd/ref.com to read: $!\n";
open(OUTPUT,">$cwd/temp.com")  or die "Cannot open $cwd/temp.com to read: $!\n";
while (my $line = <INPUT>) {
	if($delta && $line =~/( +)(\w+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
		$coord="$4$5"+$dx*$bohr;
		print OUTPUT "$1$2$3"."$coord"."$6$7$8$9$10$11\n";
		$delta=0;
	}else{
		print OUTPUT "$line";
	}
}
close(INPUT);
close(OUTPUT);

#$status=system("g09 temp.com");

open(INPUT,"<$cwd/ref.log") or die "Cannot open $cwd/ref.log to read: $!\n";
while (my $line = <INPUT>) {
	if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
		$energy=$1;
	}
	$force=($energy-$refEnergy)/$dx;#Fx
}
close(INPUT);

$forces=$force;


open(INPUT,"<$cwd/ref.com") or die "Cannot open $cwd/ref.com to read: $!\n";
open(OUTPUT,">$cwd/temp.com")  or die "Cannot open $cwd/temp.com to read: $!\n";
while (my $line = <INPUT>) {
	if($delta && $line =~/( +)(\w+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
		$coord="$7$8"+$dx*$bohr;
		print OUTPUT "$1$2$3$4$5$6"."$coord"."$9$10$11\n";
		$delta=0;
	}else{
		print OUTPUT "$line";
	}
}
close(INPUT);
close(OUTPUT);

#$status=system("g09 temp.com");

open(INPUT,"<$cwd/ref.log") or die "Cannot open $cwd/ref.log to read: $!\n";
while (my $line = <INPUT>) {
	if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
		$energy=$1;
	}
	$force=($energy-$refEnergy)/$dx;#Fx
}
close(INPUT);

$forces=$forces."   ".$force;

open(INPUT,"<$cwd/ref.com") or die "Cannot open $cwd/ref.com to read: $!\n";
open(OUTPUT,">$cwd/temp.com")  or die "Cannot open $cwd/temp.com to read: $!\n";
while (my $line = <INPUT>) {
	if($delta && $line =~/( +)(\w+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
		$coord="$10$11"+$dx*$bohr;
		print OUTPUT "$1$2$3$4$5$6$7$8$9"."$coord\n";
		$delta=0;
	}else{
		print OUTPUT "$line";
	}
}
close(INPUT);
close(OUTPUT);

#$status=system("g09 temp.com");

open(INPUT,"<$cwd/ref.log") or die "Cannot open $cwd/ref.log to read: $!\n";
while (my $line = <INPUT>) {
	if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
		$energy=$1;
	}
	$force=($energy-$refEnergy)/$dx;#Fx
}
close(INPUT);

$forces=$forces."   ".$force;

print "Numerical Force (Hartree/bohr): $forces\n"


