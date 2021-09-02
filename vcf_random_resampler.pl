#!/usr/bin/perl
use strict ; use warnings;

# vcf_pop_subseter   			# by M'Óscar 
my $version = "vcf_random_resampler.pl";

############################

# Use to re-sample a vcf file, will downsample every population to an specific number of samples 
# For options, usage and other information check help information by typing the name of the program version and "help" or "--h" or so...
# vcf_random_resampler.pl -help


################################################################################
################################   PARAMETERS   ################################
# All of them can and should be set from the command line, check the help.

#my $inputname = "populations.snps.vcf";  		# input file name, should be a string either alphanumeric or alphabetic.
my $inputname = "no default";  		# input file name, should be a string either alphanumeric or alphabetic.

#number of samples per population
my $cutval = "min";

my $infoloci = 9;  		# check how many columns of information has the VCF file before the first sample column
#my $infocols = 9;  		# check how many columns of information has the VCF file before the first sample column

#output file name
my $outfile = "generate_new_outname";  		#if this is not changed from here or from command line a name with all the details (populations, number of samples, number of loci), will be generated.
#my $outfile = "no_default_outname";  		#if this is not changed from here or from command line a name with all the details (populations, number of samples, number of loci), will be generated.

#how many characters of the sample name are the population name
my $poplength = 2;  		#how many characters at the beginning of the sample names are part of the population code
#my $poplength = 3;  		#how many characters at the beginning of the sample names are part of the population code

# popmap name
my $popmap = "No default";  		
#my $popmap = "popmap";  		
#my $popmap = "No default";  		

################################################################################
################################################################################



#Help!
my %arguments = map { $_ => 1 } @ARGV;
if(exists($arguments{"help"}) || exists($arguments{"--help"}) || exists($arguments{"-help"}) || exists($arguments{"-h"}) || exists($arguments{"--h"})) {
	die "\n\n\t   $version   Help Information\n\t-------------------------------------------------\n
	This program will generate a vcf file with an specific number of samples per population.
	Number of samples kept can correspond with that from the smallest population or with a parsed value.
	It will delete from the subset any empty or monomorphic SNP.
	It is designed to work with the VCF file generated by the program \"populations\" (Stacks).\n
	\n\tCommand line arguments and defaults:
	--input / --vcf           Name (or path) of the VCF file. If none parsed will process the first vcf found at local directory.
	                          If a directory is parsed will open a vcf there and save outputs there.\n
	--infocols                [int] Number of locus information columns before the first sample. Default: $infoloci\n
	--samples                 maximum number of samples to keep per population. Samples will be chosen randomly.
	                          if no value is parsed, will take the value from the population with fewer samples.\n
	--popmap                  Alternative to --poplength; file with the sample-population information. $popmap
	                          Same format as Stacks, one sample per line, tab separated: sample_ID\\tpopulation_ID\\n\n
	--poplength               If no --popmap; number of characters at the beginning of sample names belonging to the population code.
	                          \'--poplength 4\' for Pop1_A001; \'--poplength 2\' for HR42; etc. Default: $poplength\n
	--output                  Output file name. Saved in the same directory as the input file unless a full path is provided.
	                          If no output name is provided, one will be generated.\n
	A new popmap \"popmap_new\" (not very original, I know) will be generated only with the chosen samples.\n\n\n";
}



################ PASSING ARGUMENTS


my $parsednum = 0;

use Getopt::Long;

GetOptions( "input=s" => \$inputname,    #   --input
            "vcf=s" => \$inputname,      #   --vcf
            "output=s" => \$outfile,      #   --output
            "popmap=s" => \$popmap,      #   --popmap
            "infocols=i" => \$infoloci,      #   --infocols
            "poplength=i" => \$poplength,      #   --poplength
            "samples=i" => \$cutval );   #   --pops

my $printnum = "default";
if ($cutval eq "min") { $printnum = " will be determined by the population with fewer samples."; }
else { $printnum = ": $cutval"; }

print "$version is running, check help information for more options: vcf_pop_subseter.pl -h\n\nMaximum number of samples to keep per population$printnum\n\n";





#deal with population names if there is a popmap

my %popmapnames = ();

if ($popmap ne "No default") {
	print "Reading popmap $popmap...\n";
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
		if (exists $popmapnames{$idpopmap}) { print "$idpopmap appears more than once in popmap, only the first one will be recorded\n"; } else { $popmapnames{$idpopmap} = "$poppopmap"; }
	}
	close $POPMAP;
}






#sort out if working with a vcf file or a path

use Cwd qw(cwd);
my $localdir = cwd;


use File::Basename;

my $singlefile = "no default";
my $workpath = "no default";
my $inputfile = "no default";

if ($inputname =~ /.*?\.vcf$/ || $inputname =~ /.*?\.VCF$/) { $singlefile = "yes"; $workpath = dirname($inputname); $inputfile = basename($inputname);}
elsif ($inputname eq "no default") { $workpath = $localdir; $singlefile = "no"; }
else { $workpath = $inputname; $singlefile = "no"; }

#what to print according to path found
my $printpath = "";
if($workpath eq ".") { $workpath = $localdir; } else { $printpath = " from $workpath/" }

#define vcf filename, vcf file path to work with and directory
if ($singlefile eq "no") {
	#read files
	print "No input file parsed, checking $workpath/\n";
	opendir(DIR, $workpath);						#open the directory 
	my @infiles = readdir(DIR);					#extract filenames
	closedir(DIR);
	
	#filter files
	my @vcffiles = ();
	
	foreach my $infile (@infiles) {
		next if ($infile =~ /^\.$/);				#don't use any hidden file
		next if ($infile =~ /^\.\.$/);			
		if ($infile =~ /.*?\.vcf$/) { push (@vcffiles, $infile); }		#save the vcf files
		elsif ($infile =~ /.*?\.VCF$/) { push (@vcffiles, $infile); }		#save the vcf files
	}
	
	my $filenum = scalar @vcffiles;
	if ($filenum == 0) { die "\n\n\tERROR!\n\tNo \".vcf\" or \".VCF\" files found at $workpath\nCan't work without a valid input file. Check the help information and try again: vcf_random_resampler.pl --help\n\n"; }
	else { $inputfile = $vcffiles[0]; } 
	
	if ($filenum > 1) { print "\n\nWARNING! More than one vcf file found, only \"$inputfile\" will be analysed this time.\nReading samples...\n"; }
	else { print "Reading samples from $inputfile...\n\n"; }
}
elsif ($singlefile eq "yes") { print "Reading samples from $inputfile$printpath...\n"; }

my $filepath = "$workpath/$inputfile"; 
$filepath =~ s/\/\//\//g;
$filepath =~ s/\/\//\//g;
#print "$filepath\n";








# open vcf file

open my $VCFFILE, '<', $filepath or die "\nUnable to find or open $filepath: $!\n";


my $allsamples = 0;
my %uniquepops = ();
my %popindex = ();

while (<$VCFFILE>) {
	chomp;	#clean "end of line" symbols
	next if /^$/;  		#skip if blank
	next if /^\s*$/;  		#skip if only empty spaces
	my $line = $_;  		#save line
	$line =~ s/\s+$//;  		#clean white tails in lines
	my @wholeline= split('\t', $line);  		#split columns as different elements of an array
	my $numcolumn = scalar @wholeline;
	$allsamples = $numcolumn - $infoloci;
	my $lastsample = $numcolumn -1;
	my $firstsample = $infoloci;
	
	next if /^##.*?/; #skip first rows with metadata 

	if ($wholeline[0]=~ /^#.*?/) {
		
		foreach my $col ($firstsample..$lastsample) {
			my $sample = $wholeline[$col];  		#take one column each time
			my $popname = "none";
			
			# either extract the population name from the sample name
			if ($popmap eq "No default") {
				$popname = substr($sample, 0, $poplength);
				$popmapnames{$sample}="$popname";	#create a popmap
				$uniquepops{$popname} = 42;  		#save all populations
				#save al column indexes from each population
				if (exists $popindex{$popname}) { $popindex{$popname} = "$popindex{$popname}" . "\t" . "$col"; }  else { $popindex{$popname} = "$col"; }
			} # or use the popmap
			else { 
				if (exists $popmapnames{$sample}) {
					$popname = $popmapnames{$sample}; 
					$uniquepops{$popname} = 42;  		#save all populations
					#save al column indexes from each population
					if (exists $popindex{$popname}) { $popindex{$popname} = "$popindex{$popname}" . "\t" . "$col"; }  else { $popindex{$popname} = "$col"; }
				}
				else { print "\n\tERROR! Sample $sample does not exist in popmap and will be ignored!\n"; }
			}
			
		}
	}
	
	last if ($wholeline[0] !~ /^#.*?/)
}
close $VCFFILE;

#check number of samples per population and extract minimum
my @popnames = keys %uniquepops;
my $numpop = scalar @popnames;
my $min = 999999999;
my $max = 0;
my $bign = 0;

foreach my $pop (keys %popindex) {
	my $indexes = $popindex{$pop};
	my @allindex= split('\t', $indexes);
	my $samplecount = scalar @allindex;
	$bign = $bign + $samplecount;
	#save number of samples/cols per pop
	if ($samplecount < $min) { $min = $samplecount; }
	if ($samplecount > $max) { $max = $samplecount; }
}

if ($cutval eq "min") { $cutval = $min; }

if ($bign != $allsamples) { die "\n\n\tERROR!\n\tSomething unexpected happended, it appears to be $allsamples columns with samples in the file, but when taken by population comparing with the popmap there are $bign\n\tCheck the input file and popmap, and the help information, and report me the error if everything seems alright.\n\n\n"; }

print "A total of $bign samples in $numpop populations found: @popnames \nSaples per population range from $min to $max; downsampling bigger populations to $cutval\n\nSelecting random samples per population...\n";






#select random col indexes if above number

my %selectedsamples = ();

foreach my $pop (keys %popindex) {
	
	my $indexes = $popindex{$pop};
	my @allindex= split('\t', $indexes);
	my $samplecount = scalar @allindex;
	
	#if population has too many samples
	if ($samplecount > $cutval) {
		my %rndindex = ();
		my $size = 0;
		my $null=0;
		#select random positions from the array of columns until you get the number you want
		while ($size < $cutval) {
			my $rndm = int(rand($samplecount)); #get one random position
			
			if (exists $rndindex{$rndm}) { $null++; }  else { $rndindex{$rndm} = 42; } #save it if it is new
			$size = scalar keys %rndindex; #count how many are already saved
		}
		my @prunedindexes = keys %rndindex; #get all the positions saved
		#my $k=0;
		my $keptcols = 0;
		#save the column indexes from all the random positions saved
		foreach my $index (@prunedindexes) {
			my $colnum = $allindex[$index];
			if (exists $selectedsamples{$colnum}) {my $realcol = $colnum+1; print "\n\n\tERROR!\n\tSomething unexpected happened, it seems that column $realcol appears is duplicated in the indexed reference.\nAt populations: $selectedsamples{$colnum} and $pop.\nCheck your input file and popmap, and the help information, and report me the error if everything seems alright.\n\n"; }
			else { $selectedsamples{$colnum} = $pop; }
			
			#if ($k == 0) {my $cols = $allindex[$index]; } else { $cols ="$cols" . "\t" . "$allindex[$index]" }
			#$keptcols = $cols;
			#$k++;
		}
		#add list of indexes to the hash
		#$selectedsamples{$pop} = $keptcols;
	}
	else {
		foreach my $colnum (@allindex) {
			if (exists $selectedsamples{$colnum}) {my $realcol = $colnum+1; print "\n\n\tERROR!\n\tSomething unexpected happended, it seems that column $realcol appears is duplicated in the indexed reference.\nAt populations: $selectedsamples{$colnum} and $pop.\nCheck the input file and popmap, and the help information, and report me the error if everything seems alright.\n\n"; }
			else { $selectedsamples{$colnum} = $pop; }
		}
	}
}

my $samplesfin = scalar keys %selectedsamples;

print "Done! $samplesfin samples will remain in the final dataset (maximum $cutval per population).\n\nSelecting genotypes from each loci for the selected samples, this may take some time...\n";








#open file again and prune

open my $VCFFILE2, '<', $filepath or die "\nUnable to find or open $filepath: $!\n";


my %print_file = ();
my $row = 0;
my %newpopmap = ();

while (<$VCFFILE2>) {
	chomp;	#clean "end of line" symbols
	next if /^$/;  		#skip if blank
	next if /^\s*$/;  		#skip if only empty spaces
	my $line = $_;  		#save line
	$line =~ s/\s+$//;  		#clean white tails in lines
	my @wholeline= split('\t', $line);  		#split columns as different elements of an array
	my $numcolumn = scalar @wholeline;
	$allsamples = $numcolumn - $infoloci;
	my $lastsample = $numcolumn -1;
	my $firstsample = $infoloci;
	
	#save first rows with metadata
	if ($wholeline[0]=~ /^##.*?/) { $print_file{$row} = $line; $row++; }
	else {
		my $num = 0;
		my $rowdata = "NEW";
		#add metadata from first columns
		while ($num < $infoloci) {
			my $info = $wholeline[$num];
			if ($rowdata eq "NEW") { $rowdata = "$info"; }  else { $rowdata = "$rowdata" . "\t" . "$info"; }
			$num++;
		}
		
		#sort indexes and add genotypes
		foreach my $col (sort {$a <=> $b} keys %selectedsamples) {
			my $info = $wholeline[$col];
			$rowdata = "$rowdata" . "\t" . "$info";
			if ($wholeline[0]=~ /^#.*?/) { my $newpop = $popmapnames{$info}; $newpopmap{$info} = "$newpop"; }
		}
		
		#add selected data to a hash indexed according to the row number
		$print_file{$row} = $rowdata;
		$row++;
	
	}
}

close $VCFFILE2;

print "Done! ";






# GET RID OF MONOMORPHIC LOCI
print "Now checking monomorphic loci...\n";

my %final_file = ();
my %discarded = ();
my $rownum = 0;
my $locicount=0;
my $totalloci=0;

foreach my $line (sort {$a <=> $b} keys %print_file) { 
	
	my $checkrow = $print_file{$line};
	$checkrow =~ s/\s+$//;
	
	if ($checkrow=~ /^#.*?/) { $final_file{$rownum} = $checkrow; $rownum++; }
	else {
		my @wholerow = split('\t', $checkrow);
		my $numcolumn = scalar @wholerow;
		my $firstsample = $infoloci;
		my $lastsample = $numcolumn -1;
		
		my %genotypelist = ();
		my $rep = 1;
		
		foreach my $count ($firstsample..$lastsample) {
			
			my $allinfo = $wholerow[$count];
			my @holdallele = split (':', $allinfo);  		#split the alleles from the rest of the info
			my $genot = $holdallele[0];  		#save the allele information
			$genotypelist{$genot} = "$rep";
			$rep++;
		}
		
		my @genlist = keys %genotypelist;
		my $divers = scalar @genlist;
		my $locusname = $wholerow[2];
		my $locusinfo = "$wholerow[0]" ."\t" . "$wholerow[1]" . "\t" . "$wholerow[2]" . "@genlist";
		#print "genotypes from loci $locusname: @genlist  ";
		if (exists $genotypelist{'./.'} && $divers < 3 ) {
			$discarded{$locusname} = "$locusinfo";
			$totalloci++;
			#print "\tREMOVED\n";
		}
		elsif ($divers < 2 ) {
			$discarded{$locusname} = "$locusinfo";
			$totalloci++;
			#print "\tREMOVED\n";
		}
		else {
			$final_file{$rownum} = $checkrow;
			$totalloci++;
			$rownum++;
			$locicount++;
			#print "\tkept\n";
		}
		
	}
}
my $nummono = scalar keys %discarded;
if ($nummono == 0) { print "None. All $locicount kept.\n\n"; }  else { print "Done. $locicount loci kept from $totalloci.\n\n"; }







#save file

#file name and path
my $outname = "no default";
my $outpath = "no default";

if ($outfile eq "generate_new_outname") { 
	my $cleaname = $inputfile;
	$cleaname =~ s/^(.*)\.vcf/$1/;
	$cleaname =~ s/^(.*)\.VCF/$1/;
	
	$outname = "$cleaname" . "_max$cutval" . "perpop.vcf";
	$outpath = $workpath;
}
else {
	$outname = basename($outfile);
	$outpath = dirname($outfile);
	if ($outpath eq ".") { $outpath = $workpath; }
}
my $savepath = "$outpath/$outname";
$savepath =~ s/\/\//\//g;
$savepath =~ s/\/\//\//g;



#print file
print "Saving file as $outname...";

open my $OUT, '>', $savepath or die "\nUnable to create or save \"$savepath\": $!\n";

#sort and print
foreach my $line (sort {$a <=> $b} keys %final_file) { print $OUT "$final_file{$line}\n"; }
close $OUT;

#now popmap
my $poppath = "$outpath/popmap_new";
$poppath =~ s/\/\//\//g;
$poppath =~ s/\/\//\//g;

open my $POPOUT, '>', $poppath or die "\nUnable to create or save \"$poppath\": $!\n";
#sort and print
foreach my $ind (sort keys %newpopmap) { print $POPOUT "$ind\t$newpopmap{$ind}\n"; }
close $POPOUT;

#now monomorphic
my $monopath = "$outpath/monomorphic_loci";
$monopath =~ s/\/\//\//g;
$monopath =~ s/\/\//\//g;

#my $nummono = scalar keys %discarded;

if ($nummono > 0) {
	open my $MONOOUT, '>', $monopath or die "\nUnable to create or save \"$monopath\": $!\n";
	#sort and print
	foreach my $locus (sort keys %discarded) { print $MONOOUT "$discarded{$locus}\n"; }
	close $MONOOUT;
}


print "\n$version is done!\n\n";












