#!/usr/local/bin/perl
#!/usr/local/bin/perl
use strict;
use File::Find;
use Cwd;

my $regex = $ARGV[0];
my $cwd = getcwd();
my $coords="";
my $line_counter;
my @elements=("","H","He","Li","Be","B","C","N","O","F","Ne");

open(OUTPUT,">$cwd/$regex.xyz")  or die "Cannot open $cwd/$regex.xyz to read: $!\n";
opendir(DIR, $cwd);
my @FILES= readdir(DIR);
closedir(DIR);

foreach(@FILES){
  if($_ =~ /$regex/ && $_!~/\.xyz/){
    $coords="";
    $line_counter=0;
    open(INPUT,"<$cwd/$_") or die "Cannot open $cwd/$$_ to read: $!\n";
    while (my $line = <INPUT>){
      if($line =~ /(\d+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
       $coords=$coords."$elements[$1]   $3   $6   $9 \n";
       $line_counter++;
      }
    }
   print OUTPUT "$line_counter\n";
   print OUTPUT "$_\n";
   print OUTPUT "$coords";
   print  "$line_counter\n";
   print  "$_\n";
   print  "$coords";
   close(INPUT);
  }
}

close(OUTPUT);

