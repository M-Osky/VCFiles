#!/usr/bin/perl
use strict ; use warnings;

# vcf_pop_subseter   			# by M'Óscar 
my $version = "vcf_pop_subseter_v2.3.pl";

############################

# Use this script to make a subset of a VCF file
# As usual the individual tags should include the population code
# For options, usage and other information check the help typing the name of the program version and "help" or "--h" or so...
# vcf_pop_subseter.pl -help


#################################################################################################
##########################          CHANGELOG           #########################################
#################################################################################################
#####		                                                                                #####
#####				Version 2.3 (17/12/2020)                                                #####
#####		debugging monomorphic filtering for "onepop" option                             #####
#####		                                                                                #####
#####				Version 2.2 (01/12/2020)                                                #####
#####		added the option to rermove or not monomorphic loci                             #####
#####		                                                                                #####
#####				Version 2.1 (16/04/2020)                                                #####
#####		Now it works with a popmap                                                      #####
#####		                                                                                #####
#####				Version 2 (02/04/2020)                                                  #####
#####		New option to vcf file in as many files as populations (one pop per file)       #####
#####		                                                                                #####
#################################################################################################
#################################################################################################


#######################   PARAMETERS   #########################
# All of them can and should be set from the command line, check the help.

#my $inputname = "populations.snps.vcf";  		# input file name, should be a string either alphanumeric or alphabetic.
my $inputname = "populations.snps.vcf";  		# input file name, should be a string either alphanumeric or alphabetic.

#populations to subset, change from command line.
my $popstring = "allpairs";
#my $popstring = "allpairs";
#my $popstring = "PK PM";
my @populations = ("no default");  		# list of the populations we want in our subset.
#my @populations = ("PK", "PM");  		# list of the populations we want in our subset.

my $infocols = 9;  		# check how many columns of information has the VCF file before the first sample column
#my $infocols = 9;  		# check how many columns of information has the VCF file before the first sample column

#output file name
my $outfile = "generate_new_outname";  		#if this is not changed from here or from command line a name with all the details (populations, number of samples, number of loci), will be generated.
#my $outfile = "no_default_outname";  		#if this is not changed from here or from command line a name with all the details (populations, number of samples, number of loci), will be generated.

#only if you want all pairs of populations
my $poplength = 3;  		#how many characters at the beginning of the sample names are part of the population code
#my $poplength = 2;  		#how many characters at the beginning of the sample names are part of the population code

# popmap name
my $popmap = "No default";  		
#my $popmap = "popmap";  		
#my $popmap = "No default";  		

# remove monomorphic (0 = purge mono, 1 = keep mono)
my $keepmono = 0;
#my $keepmono = 0;
#my $keepmono = 1;


my $outdirect = "split";  		#directory for the oututs when doing all pairs
#################################################################################################################
#################################################################################################################




#Help!
my %arguments = map { $_ => 1 } @ARGV;
if(exists($arguments{"help"}) || exists($arguments{"--help"}) || exists($arguments{"-help"}) || exists($arguments{"-h"}) || exists($arguments{"--h"})) {
	die "\n\n\t   $version   Help Information\n\t-------------------------------------------------\n
	This program will subset a VCF file selecting only samples from the given populations.
	It can delete from the subset any empty or monomorphic SNP.
	A list with the loci saved in the subset will be saved.\n
	It is designed to work with the VCF file generated by the program \"populations\" (Stacks).
	If not using a popmap, it will recognise if a sample belongs to a given population by reading the sample name
	this means that the beginning of the sample name must match the population code:
	 Pop1_A001 (This will be sample A001 from Pop1), HR42 (sample 42 from HR), etc.\n
	\n\tCommand line arguments and defaults:\n
	--input / --vcf           Name (or path) of the VCF file. Default: $inputname\n
	--infocols                Number of locus information columns before the first sample. Default: $infocols\n
	--pops                    Codes of the populations we want in our subset: \"Pop1 Pop2\". Default: $popstring
	                          if \'--pops allpairs\' then will output all pair-wise combinations of populations.
	                          if \'--pops onepop\' then will output a different file for each population.\n
	--poplength               Only if \'--pops allpairs\' or \'--pops onepop\'. Default: $poplength
	                          number of characters at the beginning of the sample names that are the population code.
	                          \'--poplength 4\' for Pop1_A001; \'--poplength 2\' for HR42; etc. Only if no --popmap\n
	--mono                    Flag. Add this to keep monomorphic loci, otherwise they will be removed.\n
	--popmap                  Only if \'--pops allpairs\' or \'--pops onepop\'. Alternative to --poplength
	                          parse a file with the sample-population information. $popmap
	                          Same format as Stacks, one sample per line, tab separated: sample_ID\\tpopulation_ID\\n\n
	--output                  Output file name. Saved in the same directory as the input file unless a full path is provided.
	\t\t\t  If set up to do multiple combinations, files will be named automatically and saved in a directory named \"$outdirect\"
	\t\t\t  If no output name is provided, one with relevant information will be generated:
	\t\t\t  example: inputname_42x2k_pop1_112x6k.vcf\n
	Command line call example:\n\tvcf_pop_subseter --input /home/refmap/populations.snps.vcf --pops \"LEMUR ATLAN LAPUT\" --infocols 8 --output 3POPsubset.vcf\n\n\n";
}

my $rawname = $inputname;
$rawname =~ s/^(.*?)\.vcf/$1/;

################ PASSING ARGUMENTS


use Getopt::Long;

GetOptions( "input=s" => \$inputname,    #   --input
            "vcf=s" => \$inputname,      #   --vcf
            "output=s" => \$outfile,      #   --output
            "popmap=s" => \$popmap,      #   --popmap
            "infocols=i" => \$infocols,      #   --infocols
            "poplength=i" => \$poplength,      #   --poplength
            "mono" => \$keepmono,      #   --mono
            "pops=s" => \$popstring );   #   --pops


if ($popstring ne "allpairs" && $popstring ne "onepop") { @populations = split(' ' , $popstring); }
elsif ($popstring eq "allpairs") { print "No populations names parsed, $version will output all unique combinations of two populations.\ncheck the help information if needed:   \t $version help\n\n"; }
elsif ($popstring eq "onepop") { print "No populations names parsed, $version will output one vcf file for each population.\ncheck the help information if needed: $version help\n"; }


my $alltogethernow = join ("-", @populations);



############### DIRECTORY PATH

use Cwd qw(cwd);
my $localdir = cwd;

my @directorypath = split('/' , $inputname);
my $pathlength = scalar @directorypath;
my $filename="nofilename";
my $subdirpath="nodirpath";
my $filepath="nofilepath";

if ($pathlength > 1) {
	$filename = $directorypath[-1];
	$filepath = $inputname;
	pop (@directorypath);
	$subdirpath = join ('/' , @directorypath);
}
elsif ($pathlength <= 1) {
	$filename = $inputname;
	$filepath = "$localdir" . "/" . "$inputname";
	$subdirpath = $localdir;
}

#print "\nPopulations to extract: $populations[0], $populations[1], $populations[2]\n";


my %popmapnames = ();

if ($popmap ne "No default") {
	open my $POPMAP, '<', $popmap or die "\nUnable to find or open $popmap: $!\n";
	while (<$POPMAP>) {
		chomp;	#clean "end of line" symbols
		next if /^$/;  		#skip if blank
		next if /^\s*$/;  		#skip if only empty spaces
		my $line = $_;  		#save line
		$line =~ s/\s+$//;  		#clean white tails in lines
		my @wholeline= split('\t', $line);  		#split columns as different elements of an array
		my $idpopmap = $wholeline[0];
		my $poppopmap = $wholeline[1];
		if (exists $popmapnames{$idpopmap}) { print "\n$idpopmap appears more than once in the popmap, only the first one will be recorded\n"; } else { $popmapnames{$idpopmap} = "$poppopmap"; }
	}
	close $POPMAP;
}




print "\nReading file and checking sample names\n";



open my $VCFFILE, '<', $filepath or die "\nUnable to find or open $filepath: $!\n";
#print "$filename read succesfully.\n";




if ($popstring eq "onepop") {
	my @keepdata = ();
	my $infoloci = $infocols - 1;
	my $firstsample = $infocols;
	my @extrainfo =();
	my @lociinfo = (0..$infoloci);
	my @metadata = ();
	my $allsamples = 0;
	my %popindex = ();
	my %uniquepops = ();
	my %popmapvcf = ();
	my $headers="not";
	my $numpop = 0;
	my $rawname = $filename;

	print "Checking populations in the file... ";
	while (<$VCFFILE>) {
		chomp;	#clean "end of line" symbols
		next if /^$/;  		#skip if blank
		next if /^\s*$/;  		#skip if only empty spaces
		next if $headers eq "read";
		my $line = $_;  		#save line
		$line =~ s/\s+$//;  		#clean white tails in lines
		my @wholeline= split('\t', $line);  		#split columns as different elements of an array
		my $numcolumn = scalar @wholeline;
		$allsamples = $numcolumn - $infoloci;
		my $lastsample = $numcolumn -1;
		
		if ($wholeline[0]=~ /^##.*?/) {
			push (@metadata, $line);  		#save the first rows with metadata (they all start with "##")
			#print "Saving metadata $wholeline[0]\n";
		}
		elsif ($headers eq "not") {
			my $count=0;
			my @keepline = ();
			foreach my $col ($firstsample..$lastsample) {
				my $sample = $wholeline[$col];  		#take one column each time
				my $popname = "none";
				
				if ($popmap eq "No default") { $popname = substr($sample, 0, $poplength); } 		#extract the population name from the sample name
				else { $popname = $popmapnames{$sample}; }
				
				$uniquepops{$popname} = 42;  		#save all populations
				$popmapvcf{$sample} = $popname;  		#save each sample-population pair present inn the vcf file
				#save indexes of the columns with samples from the chosen population 
				if (exists $popindex{$popname}) { $popindex{$popname} = "$popindex{$popname}" . "\t" . "$col"; }  else { $popindex{$popname} = "$col"; }
			}
			my @popnames = keys %uniquepops;
			$numpop = scalar @popnames;
			print "\n$numpop populations found: @popnames \n\n";
			
			$headers = "read";
		}
	}
	close $VCFFILE;
	print "Now processing the samples from each population...\n";
	
	if ($keepmono == 0) { print "monomorphic loci will be removed from all populations\n\n"; }
	elsif ($keepmono == 1) { print "monomorphic loci will kept in all populations\n\n"; }
	
	my $newfolder = "$subdirpath" . "/" . "$outdirect";
	unless(-e $newfolder or mkdir $newfolder) {die "Unable to create the directory \"$newfolder\"\nMay be you don't have the rights: $!\n"; }
	
	my $k = 0;
	# Now create a file for each population
	foreach my $popname (keys %uniquepops) {
		my $keeptrack = $k+1;
		$headers = "none";
		
		print "$keeptrack of $numpop ($popname) procesed: ";
		#take the indexes from population
		my $allindexes = $popindex{$popname};
		#if ($k==0) { print "\npairs($key) = $allindexes\n"; }
		my @indexes = split("\t", $allindexes);		#new array with all indexes from this population
		my $samplenum = scalar @indexes;
		my $delete = 0;
		my $locinum = 0;
		my $n = 0;
		my @keepinfo=();
		
		my @extralines = @metadata;
		
		open my $VCFFILE, '<', $filepath or die "\nUnable to find or open $filepath: $!\n";
		while (<$VCFFILE>) {
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			next if /##.*?/;  		#skip if commented
			my $line = $_;  		#save line
			$line =~ s/\s+$//;  		#clean white tails in lines
			my @wholeline= split('\t', $line);  		#split columns as different elements of an array
			#if ($k == 0) { print "\nHeaders have been $headers\n"; }
			
			
			
			
			my @colindex = (@lociinfo,@indexes);  		#all indexes we want
			
			
			
			if ($headers eq "none") {
				my @keepheads=();
				foreach my $index (@colindex) { push (@keepheads, $wholeline[$index]); }  		#save the columns we are interested in
				my $headerskept = join ("\t", @keepheads);
				push (@extralines, $headerskept); #keepdata has now the first "##" lines and the headers
				$headers = "done";
				#if ($k == 0 && $locinum ==1) { print "metadata saved, "; }
			} elsif ($headers eq "done") {
				my @genotypes=();
				$n = 0;
				my @keepcols = ();
				#get all columns with the genotypes we want and the loci metadata
				foreach my $index (@colindex) { push (@keepcols, $wholeline[$index]); $n++; }  		#save the columns we are interested in
				#if ($k == 0 && $locinum ==1) { print "samples from $popname selected, "; }
				my $saveline = join("\t", @keepcols);
				
				if ($keepmono == 0) {
					
					#get only the column with the genotypes we want
					foreach my $cols (@indexes) { push (@genotypes, $wholeline[$cols]); }
					
					
					#get how many genotypes there are so you can loop through them
					my $numgenotypes = scalar @genotypes;
					my $maxindex = $numgenotypes - 1;
					my @indexgeno = (0..$maxindex);
					
					#replace genotypes and metadata with only genotypes
					foreach my $genot (@indexgeno) {
						my @infosample = split (':', $genotypes[$genot]);
						$genotypes[$genot] = $infosample[0];
					}
					
					
					#calculate how many different there are
					my %allgenotypes = map { $_ => 1 } @genotypes;
					my $variability =scalar keys %allgenotypes;
					
					#if ($k == 0 && $locinum ==1) { print "monomorphic pruned, "; }
					#filter out monomorphic
					if($variability==1) { $delete++; }
					elsif(exists($allgenotypes{"./."}) && $variability == 2) { $delete++; }
					else { push (@keepinfo, $saveline); }  		#add row
					
					
					
				} elsif ($keepmono == 1) {
					
					
					
					push (@keepinfo, $saveline);
					
					
				} else { die "\n\nERROR\nSomething failed when checking monomorphic loci, variable assigned to --mono flag has unexpected value\nThis is a bug, sorry.\n\n"; }
				
				#if ($k == 0 && $locinum ==1) { print "loci saved\n"; }
				$locinum++;
				
				
			}
			
		}
		
		
		my $popsamples = $n - $infocols;
		
		
		my $finalloci = $locinum - $delete;
		
		
		my $lociprint = $finalloci;
		if ($lociprint > 999) { $lociprint=~ s/^([0-9]*)[0-9][0-9][0-9]$/$1k/; }
		
		
		my $outname = "$popname$samplenum" . "x" . "$lociprint" . "_$rawname";
		my $outpath = "$newfolder" . "/" . "$outname";
		
		my @printinfo = (@extralines, @keepinfo);
		
		open my $OUT, '>', $outpath or die "\nUnable to create or save \"$outpath\": $!\n";
		foreach (@printinfo) {print $OUT "$_\n";} # Print each entry in our array to the file
		close $OUT; 
		
		print "$samplenum samples and $finalloci SNPs from $popname saved in $outname.\n";
		
		close $VCFFILE;
		$k++;
	}
	print "\nAll populations done!\n\n";
}
elsif ($popstring eq "allpairs") {
	my @keepdata = ();
	my $headers = "not";
	my $infoloci = $infocols - 1;
	my $firstsample = $infocols;
	my %pairs =();
	my %pairheads =();
	my @extrainfo =();
	my @lociinfo = (0..$infoloci);
	my @metadata = ();
	my $allsamples = 0;
	my $pairnum = 0;

	print "Checking populations in the file... ";
	while (<$VCFFILE>) {
		chomp;	#clean "end of line" symbols
		next if /^$/;  		#skip if blank
		next if /^\s*$/;  		#skip if only empty spaces
		next if $headers eq "read";
		my $line = $_;  		#save line
		$line =~ s/\s+$//;  		#clean white tails in lines
		my @wholeline= split('\t', $line);  		#split columns as different elements of an array
		my $numcolumn = scalar @wholeline;
		$allsamples = $numcolumn - $infoloci;
		my $lastsample = $numcolumn -1;
		
		if ($wholeline[0]=~ /^##.*?/) {
			push (@metadata, $line);  		#save the first rows with metadata (they all start with "##")
			#print "Saving metadata $wholeline[0]\n";
		}
		elsif ($headers eq "not") {
			#my @extrainfo = ();
			#my $locusinfo = 0;
			#$allsamples = $numcolumn - $infocols;
			#foreach my $count (0..$infoloci) {
			#	push (@extrainfo, $wholeline[$count]);  		#save the first columns that have information about the loci
			#}
			#$locusinfo = join ("\t", @extrainfo);
			my $count=0;
			#my @keepcol = ();
			my @keepline = ();
			my %uniquepops = ();
			foreach my $col ($firstsample..$lastsample) {
				my $sample = $wholeline[$col];  		#take one column each time
				my $popname = "none";  		#extract the population name from the sample name
				
				if ($popmap eq "No default") { $popname = substr($sample, 0, $poplength); } 		#extract the population name from the sample name
				else { $popname = $popmapnames{$sample}; }
				
				$uniquepops{$popname} = 42;  		#save all populations
			}
			my @popnames = keys %uniquepops;
			my $numpop = scalar @popnames;
			print "\n$numpop populations found: @popnames \n";
			
			my $first = 0;
			my $second = 0;
			my $max1 = $numpop - 1;
			my $max2 = $numpop;
			
			until ($first == $max1) {
				$second = $first + 1;
				until ($second == $max2) {
					my $combi = "$popnames[$first]"."-"."$popnames[$second]";
					$pairs{$combi} = 1;
					#push (@pairs, $combi);
					$second++;
				}
				$first++;
			}
			$pairnum = scalar keys %pairs;
			%pairheads = %pairs;
			my @pairprint = keys %pairs;
			print "$pairnum unique combinations found. ";
			
			
			
			foreach my $key (keys %pairs) {
				my @both = split("-", $key);
				my $pop1 = $both[0];
				my $pop2 = $both[1];
				my $samplelist = "0";
				my $indexlist = "0";
				
				#loop through all the samples
				foreach my $col ($firstsample..$lastsample) {
					my $sample = $wholeline[$col];  		#take one column each time
					if ($sample =~ /^$pop1/ || $sample =~ /^$pop2/) {
						# if there is match with any of the populations, save them in a variable
						if ($samplelist eq "0") {
							$samplelist = $wholeline[$col];
							$indexlist = $col;
						}
						else {
							$samplelist = "$samplelist" . "\t" . "$wholeline[$col]";
							$indexlist = "$indexlist" . "\t" . "$col";
						}
					}
				}
				#once all samples are processed, add them to the hash
				$pairs{$key} = "$indexlist\n";
				$pairheads{$key} = "$samplelist\n";
			}
			$headers = "read";
		}
	}
	close $VCFFILE;
	print "Now processing the samples for each pair...\n";
	
	my $newfolder = "$subdirpath" . "/" . "$outdirect";
	unless(-e $newfolder or mkdir $newfolder) {die "Unable to create the directory \"$newfolder\"\nMay be you don't have the rights: $!\n"; }
	
	my $k = 0;
	# Now create a file for each pair
	foreach my $key (keys %pairs) {
		my $keeptrack = $k+1;
		print "Processing $keeptrack of $pairnum ($key): ";
		#take the indexes from each pair of populations
		my $allindexes = $pairs{$key};
		#if ($k==0) { print "\npairs($key) = $allindexes\n"; }
		my @indexes = split("\t", $allindexes);		#new array with all the pairs from the index
		my $samplenum = scalar @indexes;
		my $delete = 0;
		my $locinum = 0;
		
		open my $VCFFILE, '<', $filepath or die "\nUnable to find or open $filepath: $!\n";
		while (<$VCFFILE>) {
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			next if /##.*?/;  		#skip if commented
			my $line = $_;  		#save line
			$line =~ s/\s+$//;  		#clean white tails in lines
			my @wholeline= split('\t', $line);  		#split columns as different elements of an array
			#if ($k == 0) { print "\nHeaders have been $headers\n"; }
			my @colindex = (@lociinfo,@indexes);  		#all indexes we want
			#if ($k == 0) { print "\nLociinfo indexes: @lociinfo\n"; }
			#if ($k == 0) { print "\nIndexes indexes: @indexes\n"; }
			if ($headers eq "read") {
				my @keepheads=();
				foreach my $index (@colindex) { push (@keepheads, $wholeline[$index]); }  		#save the columns we are interested in
				#if ($k == 0) { print "\nHeaders: @keepheads"; }
				my $headerskept = join ("\t", @keepheads);
				@keepdata = @metadata;  		#keepdata has now the metadata (##)
				push (@keepdata, $headerskept); #keepdata has now the first "##" lines and the headers
				$headers = "done";
			}
			elsif ($headers eq "done") {
				my @genotypes = ();
				my @rowinfo = ();
				foreach my $index (@colindex) { push (@rowinfo, $wholeline[$index]); }  		#save the columns we are interested in
				foreach my $cols (@indexes) { push (@genotypes, $wholeline[$cols]); }  		#save only columns with the genotypes
				#if ($k == 0 && $locinum == 0) { print "\nFirst locus: @genotypes"; }
				my $keptline = join("\t", @rowinfo);
				my $numgenotypes = scalar @genotypes;
				my $maxindex = $numgenotypes - 1;
				my @indexgeno = (0..$maxindex);
				
				foreach my $genot (@indexgeno) {
					my @infosample = split (':', $genotypes[$genot]);
					$genotypes[$genot] = $infosample[0];
				}
				
				#if ($locinum <= 4) {
				#	foreach (@genotypes) {print "$_   ";}
				#	print "\n\n";
				#}
				my %allgenotypes = map { $_ => 1 } @genotypes;
				my $variability =scalar keys %allgenotypes;
				
				
				
				#monomorf delete
				if ($keepmono == 0) {
					if($variability==1) { $delete++; }
					elsif(exists($allgenotypes{"./."}) && $variability == 2) { $delete++; }
					else { push (@keepdata, $keptline); }  		#add row
				} elsif ($keepmono == 1) {
					push (@keepdata, $keptline);
				} else { die "\n\nERROR\nSomething failed when checking monomorphic loci, variable assigned to --mono flag has unexpected value\nThis is a bug, sorry.\n\n"; }
				
				$locinum++;
				#$k++;
				
			}
		}
		
		my $finalloci = $locinum - $delete;
		print "$samplenum samples and $finalloci loci kept. ";
		
		
		$rawname = $filename;
		$rawname =~ s/^(.*?)\.vcf/$1/;
		my @outarg = split('/' , $outfile);
		my $outlength = scalar @outarg;
		my $outname = "no_default_outname";
		my $outpath = "no_default_outpath";
		my $lociprint = $finalloci;
		if ($lociprint > 999) { $lociprint=~ s/^([0-9]*)[0-9][0-9][0-9]$/$1k/; }
		my $alllociprint = $locinum;
		if ($alllociprint > 999) { $alllociprint=~ s/^([0-9]*)[0-9][0-9][0-9]$/$1k/; }
		
		if($rawname=~ /.*?raw.*?/) {
			$rawname =~ s/_raw//;
			$rawname =~ s/raw_//;
			$rawname =~ s/raw//;
			$outname = "$key$samplenum" . "x$lociprint" . "_$rawname" . "_$allsamples" . "x$alllociprint" . "_raw.vcf";
		}
		else {
			$outname = "$key$samplenum" . "x$lociprint" . "_$rawname" . ".vcf";
		}
		$outpath = "$newfolder" . "/" . "$outname";
		
		open my $OUT, '>', $outpath or die "\nUnable to create or save \"$outpath\": $!\n";
		foreach (@keepdata) {print $OUT "$_\n";} # Print each entry in our array to the file
		close $OUT; 
		
		print "Saved as $outname. ";
		
		my @newlist =();
		
		#save loci list
		open my $VCFOUT, '<', $outpath or die "\nUnable to find or open $outpath: $!\n";
		while (<$VCFOUT>) {
			chomp;	#clean "end of line" symbols
			
			next if /^(\s*(#.*)?)?$/;   # skip blank lines and comments
			
			my $line = $_;
			$line =~ s/\s+$//;		#clean white tails in lines
			
			my @newline = split('\t', $line);	#split columns as different elements of an array
			
			my $scaffold = $newline[0];
			$scaffold =~ s/\D//g;
			
			my $ref = $newline[2];
			my @refnum = split(':', $ref);
			my $new_name = "$scaffold" . "_" . "$refnum[0]";
			push (@newlist, $new_name);
		}
		close $VCFOUT;

		

		my $listpath = "$newfolder" . "/" . "$key". "_sublocilist";
		open my $LIST, '>', $listpath or die "\nUnable to create or save \"$listpath\": $!\n";
		foreach (@newlist) {print $LIST "$_\n";} # Print each entry in our array to the file
		print " Loci list saved!\n";
		close $LIST; 
	$headers = "read";
	
	close $VCFFILE;
	$k++;
	#if ($key eq "BD-VT" || $key eq "VT-BD") { die "\n\nCheck the outputs\n\n"; }
	}
	print "\nAll pairs done!\n";
}
elsif ($popstring ne "allpairs" && $popstring ne "onepop") {
	my @keepdata = ();
	my $firstsample = $infocols;
	my $infoloci = $infocols - 1;
	my @lociinfo = (0..$infoloci);
	my $headers = "not";
	my @colindex = ();
	my $samplenum="none";
	my $numpop = scalar @populations;
	my @keepsample=();
	my $locinum = 0;
	my $allsamples = 0;
	my $delete = 0;
	my @keepcol=();
	print "Looking for samples from populations: @populations\n";
	while (<$VCFFILE>) {
		chomp;	#clean "end of line" symbols
		next if /^$/;  		#skip if blank
		next if /^\s*$/;  		#skip if only empty spaces
		my $line = $_;  		#save line
		$line =~ s/\s+$//;  		#clean white tails in lines
		my @wholeline= split('\t', $line);  		#split columns as different elements of an array
		my $numcolumn = scalar @wholeline;
		my $lastsample = $numcolumn -1;
		
		if ($wholeline[0]=~ /^##.*?/) {
			push (@keepdata, $line);  		#save the first rows with metadata (they all start with "##")
			#print "Saving metadata $wholeline[0]\n";
		}
		elsif ($headers eq "not") {
			my @extrainfo = ();
			my $locusinfo = 0;
			$allsamples = $numcolumn - $infocols;
			foreach my $count (0..$infoloci) {
				push (@extrainfo, $wholeline[$count]);  		#save the first columns that have information about the loci
			}
			$locusinfo = join ("\t", @extrainfo);
			my $count=0;
			@keepcol = ();
			my @keepline = ();
			foreach my $col ($firstsample..$lastsample) {
				my $sample = $wholeline[$col];  		#take one column each time
				foreach (@populations) {
					my $checkpop = $_;  		# take one population code each time
					if ($sample =~ /^$checkpop/) {
						push (@keepsample, $wholeline[$col]);  		#if sample names matches popname, keep it
						push (@keepcol, $col);  		# keep also the position of the column with that sample
					}
				}
			}
			$samplenum = scalar @keepsample;  		#count number of samples kept
			@colindex = (@lociinfo,@keepcol);  		#save the positions of the columns to keep
			@keepline = (@extrainfo, @keepsample);
			my $headerline = join ("\t", @keepline);
			#print "\n\nSaving headers:\n$headerline\n";
			push (@keepdata, $headerline);
			$headers="done";
			if ($samplenum == 0) { die "\n\n\tERROR\n\n\tNone of the samples found in the vcf file matched any of the populations: @populations\n\n\n\tIn the mode ran (population list parsed with --pops) populations need to match the beginning of sample name.\n\tCheck the help information!\n\n\n"; }
		}
		elsif ($headers="done") {
			my @genotypes = ();
			my @rowinfo = ();
			foreach (@colindex) {
				push (@rowinfo, $wholeline[$_]);  		#save the columns we are interested in
			}
			foreach (@keepcol) {
				push (@genotypes, $wholeline[$_]);  		#save only columns with the genotypes
			}
			my $keptline = join("\t", @rowinfo);
			my $numgenotypes = scalar @genotypes;
			my $maxindex = $numgenotypes - 1;
			my @indexgeno = (0..$maxindex);
			
			foreach (@indexgeno) {
				my @infosample = split (':', $genotypes[$_]);
				#$genotypes[$_]=~ s/^([0-9]\/[0-9])(:.*?)/$1/;  		#save only the genotype
				$genotypes[$_] = $infosample[0];
			}
			
			#if ($locinum <= 4) {
			#	foreach (@genotypes) {print "$_   ";}
			#	print "\n\n";
			#}
			my %allgenotypes = map { $_ => 1 } @genotypes;
			my $variability = scalar keys %allgenotypes;
			#print "$variability ";
			
			#monomorf delete
			if ($keepmono == 0) {
				if($variability==1) { $delete++; }
				elsif(exists($allgenotypes{"./."}) && $variability == 2) { $delete++; }
				else { push (@keepdata, $keptline); }  		#add row
			} elsif ($keepmono == 1) {
				push (@keepdata, $keptline);
			} else { die "\n\nERROR\nSomething failed when checking monomorphic loci, variable assigned to --mono flag has unexpected value\nThis is a bug, sorry.\n\n"; }
			
			
			$locinum++;
		}
	}

	close $VCFFILE;

	my $finalloci = $locinum - $delete;

	$rawname = $filename;
	$rawname =~ s/^(.*?)\.vcf/$1/;


	my @outarg = split('/' , $outfile);
	my $outlength = scalar @outarg;
	my $outname = "no_default_outname";
	my $outpath = "no_default_outpath";

	my $lociprint = $finalloci;
	if ($lociprint > 999) { $lociprint=~ s/^([0-9]*)[0-9][0-9][0-9]$/$1k/; }
	my $alllociprint = $locinum;
	if ($alllociprint > 999) { $alllociprint=~ s/^([0-9]*)[0-9][0-9][0-9]$/$1k/; }



	if ($outfile eq "generate_new_outname") {
		if($rawname=~ /.*?raw.*?/) {
			$rawname =~ s/_raw//;
			$rawname =~ s/raw_//;
			$rawname =~ s/raw//;
			$outname = "$alltogethernow$samplenum" . "x$lociprint" . "_$rawname" . "_$allsamples" . "x$alllociprint" . "_raw.vcf";
		}
		else {
			$outname = "$alltogethernow$samplenum" . "x$lociprint" . "_$rawname" . ".vcf";
		}
		$outpath = "$subdirpath" . "/" . "$outname";
	}
	elsif ($outlength > 1) { $outpath = $outname; }
	elsif ($outlength <= 1) { $outpath = "$subdirpath" . "/" . "$outfile" }



	print "All samples processed!\nOriginal file had $locinum loci and $allsamples samples... $delete monomorphic or \"empty\" loci loci deleted.\nA subset of $finalloci loci and $samplenum samples from $numpop populations (@populations) was done.\nSaving $outpath\n";
	open my $OUT, '>', $outpath or die "\nUnable to create or save \"$outpath\": $!\n";
	foreach (@keepdata) {print $OUT "$_\n";} # Print each entry in our array to the file
	close $OUT; 

	my @newlist =();

	print "Saving list of loci...";
	open my $VCFOUT, '<', $outpath or die "\nUnable to find or open $outpath: $!\n";
	while (<$VCFOUT>) {
		chomp;	#clean "end of line" symbols
		
		next if /^(\s*(#.*)?)?$/;   # skip blank lines and comments
		
		my $line = $_;
		$line =~ s/\s+$//;		#clean white tails in lines
		
		my @newline = split('\t', $line);	#split columns as different elements of an array
		
		my $scaffold = $newline[0];
		$scaffold =~ s/\D//g;
		
		my $ref = $newline[2];
		my @refnum = split(':', $ref);
		my $new_name = "$scaffold" . "_" . "$refnum[0]";
		push (@newlist, $new_name);
	}
	close $VCFOUT;

	print "  done.\n";

	my $listpath = "$subdirpath" . "/" . "$alltogethernow". "_sublocilist";
	open my $LIST, '>', $listpath or die "\nUnable to create or save \"$listpath\": $!\n";
	foreach (@newlist) {print $LIST "$_\n";} # Print each entry in our array to the file
	close $LIST; 

}

print "$version finished!\n";


