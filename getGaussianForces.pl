    #!/usr/local/bin/perl
    use Cwd;
    
    my $file_name = $ARGV[0];
    my $cwd = getcwd();
    my $seperator=0;
    
    
    my $fp="%16.10e";
    my $fn="%16.9e";
    my $force;
    my $forces;
    my $energy;
    my $hartree_to_ev=27.2113838668;
    my $bohr_to_angstrom=0.52917721092;
    
    my $factor=$hartree_to_ev/$bohr_to_angstrom;
    #my $factor=1;
    my $line;
    my $line2;
    my @forceout;
    my @forceread;
    my $i;
    my $checkread=1;
    my $ghost_counter=0;
    my $error_thresh=10**(-14);
    my $mol_counter=0;
    my @mol_numbers=(6,1,1,1,6,1,1,1);
    my $force_counter=0;
    my $n_forces=1;
	
    
    my $status=system("g09 temp.com");
    $status=system("cp temp.log log.temp");
    $status=system("cat temp.log >> gaussian.logs");
    
    printForces("temp.forces",readGaussianForces("temp.log"));

    
    sub readGaussianForces{
	my @force;
	my @forces;
	my $mol_counter=1;
	open(INPUT,"<$_[0]") or die "Cannot open $_[0] to read: $!\n";
	while (($line = <INPUT>) && $seperator<2) {
	if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
	    $energy=$1*$hartree_to_ev;
	    $forces[0]=$energy;
	} 
	if($line =~ /Forces/){
	    #We have entered the force section of the output
	    while(($line = <INPUT>) && $seperator<2){
		
		if($line =~ /-{2,}/){
		    $seperator++;	
		}
		else{
		    @force=split(" ",$line);
		    if($force[1]!=0){
		      @force=($force[2],$force[3],$force[4]);
		      #print "@force\n";
		      $forces[$mol_counter]=[@force];
		      $mol_counter++;
		    }  
		}
	    }
	}
    }
	close(INPUT);
	#print "$#forces\n";
	return @forces;
    }
    
    sub printForces{
	my @force;
	my $mol_counter=2;
	my $i;
	open(OUTPUT,">$_[0]") or die "Cannot open $_[0] to write: $!\n";
	print OUTPUT "$_[1]\n";
	for $ref ( @_ ) {
		if (@$ref!="") {
		    $i=0;
		    for $coord (@$ref){
			print OUTPUT "$coord";
			if ($i==2) {
			    print OUTPUT "\n";
			}else{
			    print  OUTPUT "   "
			}
			$i++;
		    }
		}
		
		
	}	
	close(OUTPUT);
    }

