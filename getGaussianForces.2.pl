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
    if (($checkread==1)&&(-e "temp.forces")&&(-e "forces.read") ){
	
	
	open(INPUT1,"<$cwd/temp.forces")  or die "Cannot open $cwd/temp.forces to read: $!\n";
	open(INPUT2,"<$cwd/forces.read")  or die "Cannot open $cwd/forces.read to read: $!\n";
	open(ERROR,">>$cwd/forces.err")  or die "Cannot open $cwd/forces.error to write: $!\n";
	#    $line=<INPUT1>;
	#    $line2=<INPUT2>;
	#    print "DIFF(PERL PRINTED AND FORTRAN PRINTED FORCES + ENERGY)\n";
	#    print "$line"-"$line2"."\n";
	#    while($line=<INPUT1>){
	#    $line2=<INPUT2>;
	#    @forceout=split(" ",$line);
	#    @forceread=split(" ",$line2);
	#    $i=0;
	#      while($i<=$#forceread){
	#	  my $diff=$forceread[$i]-$forceout[$i];
	#	       if($diff>0){
	#		  print sprintf($fp,"$diff");
	#	       }else{
	#		  print sprintf($fn,"$diff");
	#	       }
	#	print "   ";
	#	$i=$i+1;
	#      }
	#    print "\n";
	#    }
	close(INPUT1);
	open(INPUT3,"<$cwd/temp.log") or die "Cannot open $cwd/temp.log to read: $!\n";
	seek(INPUT2,0,0);
	#	print "DIFF(LOGGED AND FORTRAN PRINTED FORCES + ENERGY)\n";
	$seperator=0;
	$ghost_counter=0;
	$mol_counter=0;
	while (($line = <INPUT3>) && $seperator<2) {
	    if($line =~ /E\(\w+\)\s+=\s+(\S+)/){
		$line2=<INPUT2>;
		$energy=$1*$hartree_to_ev;
	    } 
	    elsif($line =~ /Forces/){
		#We have entered the force section of the output
		#read in extra table header line
		$line = <INPUT3>;
		while(($line = <INPUT3>) && $seperator<2){
		    @forceout=split(" ",$line);
		    if($line =~ /-{2,}/){
			$seperator++;	
		    }elsif($forceout[1]!=0){
			if ($forceout[0]=~/\./ || $forceout[1]=~/\./ || $forceout[2]!~/\./ || $forceout[3]!~/\./ || $forceout[4]!~/\./) {
			    print ERROR "WOOOOOOOOOOOOAHHHHHHHHH SOMETHING WENT WROOOOOOOOOOOOOOONNNNNNNNNG!!!!!!!!!\n\nFORCE LINE DOES NOT HAVE CORRECT FORMAT\n\n";
			}
			
			$line2=<INPUT2>;
			@forceread=split(" ",$line2);
			$i=0;
			while($i<=2){
			    my $diff=$forceread[$i]-$forceout[$i+2]*$factor;
			    if(($diff!=0) && (abs("$diff")>="$error_thresh")){
				print ERROR "$diff\n";
			    }
			    
			    
			    $i=$i+1;
			}
			#print "\n";
			if($forceout[1]!=$mol_numbers[$mol_counter]) {
			    print ERROR "ERROR:MISMATCH IN ATOM TYPES!!!!!!!\nread: $forceout[1] when expecting:  $mol_numbers[$mol_counter]\n";
			}
			if ($mol_counter>$#mol_numbers) {
			    print ERROR "unexpected atom!!!!!!!!\n of type:$forceout[1]\n";
			}
			$mol_counter++;
		    }elsif($forceout[0]!~/\./ && $forceout[1]!~/\./ && $forceout[2]=~/\./ && $forceout[3]=~/\./ && $forceout[4]=~/\./){
			$ghost_counter++;
			#print "ghost center\n";
		    }
		}
		$force_counter++;
		
	    }
	}
	if ($force_counter!=$n_forces) {
	    print ERROR "unexpected number of force tables:  $force_counter ,expected: $n_forces  !!!!!!!!\n";
	}
	seek(INPUT3,0,0);
	close(INPUT2);
	close(INPUT3);
	close(ERROR);
	$seperator=0;
    }
    #end checkread
    
    
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
			print "$coord";
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

