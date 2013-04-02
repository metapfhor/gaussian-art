#!/usr/local/bin/perl
#!/usr/local/bin/perl
use Cwd;

my $file_name = $ARGV[0];
my $cwd = getcwd();
my $coords="";
my $line_counter;
my @elements=("","H","He","Li","Be","B","C","N","O","F","Ne");

open(INPUT,"<$cwd/$file_name") or die "Cannot open $cwd/$file_name to read: $!\n";
open(OUTPUT,">$cwd/$file_name.xyz")  or die "Cannot open $cwd/$file_name.xyz to read: $!\n";
while (my $line = <INPUT>){
    if($line =~ /(\d+)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)( +)(-?\d+\.\d+([Ee][+-]?\d+)?)/){
      $coords=$coords."$elements[$1]   $3   $6   $9 \n";
      $line_counter++;
    }
}
print OUTPUT "$line_counter\n";
print OUTPUT "$file_name\n";
print OUTPUT "$coords";
close(OUTPUT);
close(INPUT);
