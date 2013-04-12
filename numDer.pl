#!/usr/local/bin/perl
#Calculates the force on the first atom in the system by numerical differentiation
#for checking the units of gaussian and seeing where the numerical differentiation breaks down
use Cwd;

my $file_name = $ARGV[0];
my $cwd = getcwd();
my $seperator=0;
my $bohr=0.529;
my $dx=$ARGV[0];


my $status=system("g09 ref.com");
my $energy;
my $force;
my $forces;
my $refEnegry;
my $refForces;
my $energy;
my $factor=1; 
my $delta=1;
my $coord;
my $factor=27.211;

open(INPUT,"<$cwd/ref.log") or die "Cannot open $cwd/ref.log to read: $!\n";
while (my $line = <INPUT>) {
#print $line;
	if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
		$refEnergy=$1;
		print "E0=$refEnergy\n";
	} 
	if($delta && $line =~ /(\d+)( +)(\d+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
				$forces="$5$6"."   "."$8$9"."   "."$11$12";				
				$refForces=$forces;
				#print  "Gaussian force (Hartree/bohr): $refForces\n";
				$delta=0;
			}
}
close(INPUT);

$delta=1;

open(INPUT,"<$cwd/ref.com") or die "Cannot open $cwd/ref.com to read: $!\n";
open(OUTPUT,">$cwd/tempx.com")  or die "Cannot open $cwd/temxp.com to read: $!\n";
while (my $line = <INPUT>) {
	if($delta && $line =~/( +)(\w+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
		$coord="$4$5"+($dx*$bohr);
    print "$4$5 -> $coord\n";
		print OUTPUT "$1$2$3"."$coord"."$6$7$8$9$10$11\n";
		$delta=0;
	}else{
		print OUTPUT "$line";
	}
}
close(INPUT);
close(OUTPUT);

$status=system("g09 tempx.com");

open(INPUT,"<$cwd/tempx.log") or die "Cannot open $cwd/ref.log to read: $!\n";
while (my $line = <INPUT>) {
	if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
		$energy=$1;
		print "Ex=$energy\n";
	}
	$force=$factor*($energy-$refEnergy)/$dx;#Fx
}
close(INPUT);

$forces=$force;

$delta=1;

open(INPUT,"<$cwd/ref.com") or die "Cannot open $cwd/ref.com to read: $!\n";
open(OUTPUT,">$cwd/tempy.com")  or die "Cannot open $cwd/temp.com to read: $!\n";
while (my $line = <INPUT>) {
	if($delta && $line =~/( +)(\w+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
		$coord="$7$8"+($dx*$bohr);
		#print "$7$8 -> $coord\n";
		print OUTPUT "$1$2$3$4$5$6"."$coord"."$9$10$11\n";
		$delta=0;
	}else{
		print OUTPUT "$line";
	}
}
close(INPUT);
close(OUTPUT);

$status=system("g09 tempy.com");

open(INPUT,"<$cwd/tempy.log") or die "Cannot open $cwd/ref.log to read: $!\n";
while (my $line = <INPUT>) {
	if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
		$energy=$1;
		print "Ey=$energy\n";
	}
	$force=$factor*($energy-$refEnergy)/$dx;#Fx
}
close(INPUT);

$forces="$forces   $force";

$delta=1;

open(INPUT,"<$cwd/ref.com") or die "Cannot open $cwd/ref.com to read: $!\n";
open(OUTPUT,">$cwd/tempz.com")  or die "Cannot open $cwd/temp.com to read: $!\n";
while (my $line = <INPUT>) {
	if($delta && $line =~/( +)(\w+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
		$coord="$10$11"+($dx*$bohr);
		#print "$10$11 -> $coord\n";
		print OUTPUT "$1$2$3$4$5$6$7$8$9"."$coord\n";
		$delta=0;
	}else{
		print OUTPUT "$line";
	}
}
close(INPUT);
close(OUTPUT);

$status=system("g09 tempz.com");

open(INPUT,"<$cwd/tempz.log") or die "Cannot open $cwd/ref.log to read: $!\n";
while (my $line = <INPUT>) {
	if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
		$energy=$1;
		print "Ez=$energy\n";
	}
	$force=$factor*($energy-$refEnergy)/$dx;#Fx
}
close(INPUT);

$forces="$forces   $force";
print  "Gaussian force (Hartree/bohr): $refForces\n";
print "Numerical Force (Hartree/bohr): $forces\n"


