#!/usr/bin/perl
use strict ; use warnings;

# vcf_loci_subseter   			# by M'Óskar 
my $version = "vcf_loci_subseter_v2.6.pl";

############################

# Use this script to subset a number of loci from a vcf filefield
# check the help information (--h; help; etc...) for more details: vcf_loci_subseter.pl --help




#############################################################################################
#####################		Change log	vcf_loci_subseter		#############################
#############################################################################################
#########			Version 2.6 20/07/2020									#################
#########	Fixed some bugs from when locinames in list are missing ":+" or ":-"		#####
#########																				#####
#########			Version 2.5 22/06/2020									#################
#########	Fixed some bugs, now loci-SNP separator is kept as in the original file		#####
#########	Adding metadata to the output vcf file										#####
#########																				#####
#############################################################################################








######################   PARAMETERS   #########################
# All of them can and should be set from the command line, check the help.

my $inputname = "populations.snps.vcf";  		# input file name, either a vcf file or a directory with vcf files.
#my $inputname = "populations.snps.vcf";  		# input file name, either a vcf file or a directory with vcf files.


my $random = 1;  		# you want a subset of loci to be random? then write how many random loci you want, otherwise just "0";
#my $random = 0;  		# you want a subset of loci to be random? then write how many random loci you want, otherwise just "0";
#my $random = 1;  		# you want a subset of loci to be random? then write how many random loci you want, otherwise just "0";
#my $random = 223;  		# you want the subset of loci to be random? then write how many random loci you want

my $locilist = "no default";  		# if you want a specific loci set this is the file name
#my $locilist = "no default";  		# if you do not want a specific loci set this is the file name
#my $locilist = "candidate_loci_list.txt";  		# if you want a specific loci set this is the file name

my $exclude = "default";  		#any set of loci to exclude from the random list?
#my $exclude = "default";  		#any set of loci to exclude from the random list?
#my $exclude = "no";  		#any set of loci to exclude from the random list?
#my $exclude = "bad_loci_list.txt";  		#any set of loci to exclude from the random list?

my $sep = "\n";  		# how are the loci names separated in the list?
#my $sep = "\n";  		# how are the loci names separated in the list?
#my $sep = "\t";  		# how are the loci names separated in the list?
#my $sep = " ";  		# how are the loci names separated in the list?

my $outdir = "out";		#output subdirectory for multiple files


#################################################################################################################
#################################################################################################################



my $printsep = "not defined";
if ($sep eq "\n") { $printsep = "one per line: \"\\n\""; }
elsif ($sep eq "\t") { $printsep = "tabulated: \"\\t\""; }
elsif ($sep eq " ") { $printsep = "space: \" \""; }
else { $printsep = $sep; }

#Help!
my %arguments = map { $_ => 1 } @ARGV;
if(exists($arguments{"help"}) || exists($arguments{"--help"}) || exists($arguments{"-help"}) || exists($arguments{"-h"}) || exists($arguments{"--h"})) {
	die "\n\n\t   $version   Help Information\n\t-------------------------------------------------\n
\tThis program will subset one .vcf file selecting specific loci listed in a file and/or a determinated number of random loci.
\tAlternatively instead of parsing a single file name, a directory with \'.vcf\' files can be given; will subset all of them.
\tIt is designed to work with the VCF files generated by the program \"populations\" (Stacks).
\t\n\tCommand line arguments and defaults:
\t--input / --vcf           name/path of a single .vcf file or to a directory with \'.vcf\' files. Default: $inputname
\t--locilist                name of the file with the loci names. If none is provided it will only output random loci.
\t                          loci names must match one of these formats:\n\t                             \"scaffold number _ loci position\" (1_2)\n\t                             \"loci position : SNP position\" (2:44)
\t--random                  [integer] Number of random loci to subset. No random subset produced if \'--random 0\'\n\t                           \'--random 1\' to subset the same number of loci than in the loci list; default = $random
\t                          chosen random loci names will be saved in a different file as a one name per row file.
\t--sep                     how are the loci names separated in the loci list (if any); default = $printsep
\t                          if the file is not one loci per row, the loci names must be in the first row.
\t--exclude                 list of loci to be excluded from the random subset; by default loci from locilist.
\t                           If \'--exclude no\' then no loci will be excluded.\n\t                           If a new file is provided, column separator must be the same than for locilist.
\t--out                     only if input is a directory with multiple files. Name or path for the output directory.
\t                          Default: output files and names list will be saved at the loci list file location,\n\t                          if none provided then at the vcf file location.
\t                          When multiples files default subdirectory for the outputs = $outdir.
\n\n\tCommand line call example:\n\t$version --input /usr/home/refmap/populations.snps.vcf --random 1 --locilist candidate100loci.txt --sep \"\\t\"
\n\tWARNING!\n\tLoci are selected in random order and I did not implemented yet to sort the final file\n\tyou may want to use this vcftools module to sort them by chromosome:
\n\tvcf-sort -c output_file.vcf > sorted_file.vcf\n\n\n";
}



################ PARSING ARGUMENTS

my $popstring = "not_defined";
use Getopt::Long;

GetOptions( "input=s" => \$inputname,    #   --input
            "vcf=s" => \$inputname,      #   --vcf
            "random=i" => \$random,      #   --random
            "sep=s" => \$sep,      #   --sep
            "exclude=s" => \$exclude,      #   --exclude
            "out=s" => \$outdir,      #   --out
            "locilist=s" => \$locilist );   #   --locilist

if ($random == 0 && $locilist eq "no default") { die "\n\nERROR!\nNumber of random loci to subset = $random\nList of specific loci to subset = not provided\nNot sure what do you want this program to do then... I'm not a clever program, I need some command line arguments\nCheck the help information if needed\n\t$version --help\n\t$version -h\n etc...\n\n" }




############### DIRECTORY PATH

use Cwd qw(cwd);
my $localdir = cwd;

my $filename="nofilename";
my $subdirpath="nodirpath";
my $filepath="nofilepath";
my $locifile = "no locilist";
my $locipath = "no loci list";

my @inputpath = split('/' , $inputname);
my $inputlength = scalar @inputpath;

my $singlefile = "no idea";
if ($inputname =~ /.*?\.vcf$/ || $inputname =~ /.*?\.VCF$/) { $singlefile = "yes"; }
else { $singlefile = "no"; }



if ($locilist eq "no default") {
	if ($inputlength > 1) {
		$filename = $inputpath[-1];
		$filepath = $inputname;
		pop (@inputpath);
		$subdirpath = join ('/' , @inputpath);
	}
	elsif ($inputlength <= 1) {
		$filename = $inputname;
		$filepath = "$localdir" . "/" . "$inputname";
		$subdirpath = $localdir;
	}
}
else {
	my @directorypath = split('/' , $locilist);
	my $pathlength = scalar @directorypath;

	if ($pathlength > 1) {
		$locifile = $directorypath[-1];
		$locipath = $locilist;
		pop (@directorypath);
		$subdirpath = join ('/' , @directorypath);
		if ($inputlength > 1) {
			$filename = $inputpath[-1];
			$filepath = $inputname;
		}
		elsif ($inputlength <= 1) {
			$filename = $inputname;
			$filepath = "$localdir" . "/" . "$inputname";
		}
	}
	elsif ($pathlength <= 1) {
		$locifile = $locilist;
		$locipath = "$localdir" . "/" . "$locifile";
		$subdirpath = $localdir;
		if ($inputlength > 1) {
			$filename = $inputpath[-1];
			$filepath = $inputname;
		}
		elsif ($inputlength <= 1) {
			$filename = $inputname;
			$filepath = "$localdir" . "/" . "$inputname";
		}
	}
}






###### START THE ANALYSIS


my %goodloci = ();
my $numgood=0;
my $sample=0;

# LIST OF GOOD LOCI

if ($locilist ne "no default") {
	print "\nReading loci list $locifile" . "... ";
	open my $LOCIFILE, '<', $locipath or die "\nUnable to find or open $locipath: $!\n";
	if ($sep eq "\n") {
		while (<$LOCIFILE>) {
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			next if /^#.*$/;  		#skip comments
			my $line = $_;  		#save line
			$line =~ s/\s+$//;  		#clean white tails in lines
			$goodloci {$line} = "42";
			$sample = $line;
		}
		$numgood = scalar keys %goodloci;
		print "$numgood loci names saved from list.\n";
	}
	else {
		my $gl=0;
		while (<$LOCIFILE>) {
			exit if $gl > 0;
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			my $line = $_;  		#save line
			$line =~ s/\s+$//;  		#clean white tails in lines
			my @listedloci = split ($sep, $line);
			$sample = $listedloci[-1];
			%goodloci = map {$_ => 42} @listedloci;
			$gl++;
		}
		$numgood = scalar keys %goodloci;
		print "$numgood loci names saved from list.\n";
	}
	close $LOCIFILE;
}



# LIST OF BAD LOCI

my %badloci = ();
my $numbad=0;

if ($exclude eq "default") { %badloci = %goodloci; }
elsif ($exclude ne "no") {
	print "\nReading loci to exclude $exclude" . "... ";
	open my $EXCLUDEFILE, '<', $exclude or die "\nUnable to find or open $exclude: $!\n";
	if ($sep eq "\n") {
		while (<$EXCLUDEFILE>) {
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			next if /^#.*$/;  		#skip comments
			my $line = $_;  		#save line
			$line =~ s/\s+$//;  		#clean white tails in lines
			$badloci {$line} = "36";
		}
		$numbad = scalar keys %badloci;
		print "$numbad loci names saved from list.\n";
	}
	else {
		my $bl=0;
		while (<$EXCLUDEFILE>) {
			exit if $bl > 0;
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			my $line = $_;  		#save line
			$line =~ s/\s+$//;  		#clean white tails in lines
			my @noloci = split ($sep, $line);
			%badloci = map {$_ => 36} @noloci;
			$bl++;
		}
		$numbad = scalar keys %badloci;
		print "$numbad loci names saved from list.\n";
	}
	close $EXCLUDEFILE;
}

#print "\nLast locus $sample\n";


if ($singlefile eq "yes") {
	
	my @headers = ();
	
	# SUBSETING
	
	print "One file parsed, processing $filename... \n";
	open my $VCFFILE, '<', $filepath or die "\nUnable to find or open $filepath: $!\n";



	# READ VCF FILE
	my @rawfile = ();
	my @wholefile = ();

	while (<$VCFFILE>) {
		chomp;	#clean "end of line" symbols
		next if /^$/;  		#skip if blank
		next if /^\s*$/;  		#skip if only empty spaces
		my $row = $_;  		#save line
		$row =~ s/\s+$//;  		#clean white tails in lines
		push (@rawfile, $row);
	}
	
	close $VCFFILE;
	#save headers
	foreach (@rawfile) {
	my $line = $_;
	my @wholeline= split('\t', $line);  		#split columns as different elements of an array
		if ($wholeline[0]=~ /^#.*?/) {
		push (@headers, $line);  		#save the first rows with metadata (they all start with "#")
		#print "Saving metadata $wholeline[0]\n";
		}
		else { push (@wholefile, $line); }  #save the real data; new version 2.6
	}



	# EXTRACT LOCI FROM LIST



	my @indexes = ();
	my $savedloci;
	my @namelist =();
	my $sepsymbol = 0;
	
	# FIRST LETS GET THE LOCI LIST
	if ($locilist ne "no default") {
		my @subsetlisted = ();
		my @keptlisted = ();
		my $index=0;
		foreach my $line (@wholefile) {
			my @wholeline= split('\t', $line);  		#split columns as different elements of an array
			#save the headers
			#if ($wholeline[0]=~ /^#.*?/) {
			#	push (@headers, $line);  		#save the first rows with metadata (they all start with "#")
			#	#print "Saving metadata $wholeline[0]\n";
			#}
			# save new loci names listed
			#elsif ($sample =~ /^[0-9]*_[0-9]*$/) {
			if ($sample =~ /^[0-9]*_[0-9]*$/) {
				# generate the "new locinames" from vcf
				my $scaffold = $wholeline[0];
				$scaffold =~ s/\D//g;
				my $ref = $wholeline[2];
				my @refnum = split(':', $ref);
				my $new_name = "$scaffold" . "_" . "$refnum[0]";
				# if they exist in the list keep them
				if (exists $goodloci{$new_name}) { 
					push (@keptlisted, $line);
					push (@indexes, $index);
					push (@namelist, $new_name);
				}
				$sepsymbol = "_";
			}
			elsif ($sample=~ /^[0-9]*:[0-9]*/) {
				my $fullref = $wholeline[2];   		#replaced $ref for $fullref v2.6
				my $ref = $fullref;   		#new line v2.6
				$ref =~ s/^([0-9]*:[0-9]*):.*/$1/;   		#new line v2.6
				#print "$ref\n";
				# if the loci exist in the list, keep them
				#if (grep {/^$ref/} keys %goodloci) { 
				if (exists $goodloci{$ref}) {
					push (@keptlisted, $line);
					push (@indexes, $index);
					push (@namelist, $ref);
				}
				$sepsymbol = ":";
			}
			else { die "\n\nERROR!\nLoci format does not match any of the types accepted, check help information: vcf_loci_subseter help\n\n"; }
			$index++;
		}
		$savedloci = scalar @keptlisted;
		print "$savedloci loci matching the list at $locifile were found.\n";
		
		my $selectedfile = "candidate$savedloci" . "loci_$filename";
		
		@subsetlisted = (@headers, @keptlisted);
		#save
		my $listpath = "$subdirpath" . "/" . "$selectedfile";
		print "Saving as $selectedfile... ";
		open my $LIST, '>', $listpath or die "\nUnable to create or save \"$listpath\": $!\n";
		foreach (@subsetlisted) {print $LIST "$_\n";} # Print each entry in our array to the file
		close $LIST; 
		
		my $noext = $filename;
		$noext =~ s/^(.*?).vcf/$1/;
		$noext =~ s/^(.*?).VCF/$1/;
		my $coincident = "$subdirpath" . "/" . "$noext" . "_$savedloci" . "locinames";
		open my $NAMES, '>', $coincident or die "\nUnable to create or save \"$coincident\": $!\n";
		foreach (@namelist) {print $NAMES "$_\n";} # Print each entry in our array to the file
		close $NAMES; 
		
		print "done!\n\n";
	}


	if ($random != 0) {
		my @randomloci = ();
		print "Now generating random SNP positions... ";
		# extract as many loci as loci there are in the list
		if ($random == 1) { $random = $savedloci; }
		my $metadata = scalar @headers;
		my $totalrows = scalar @wholefile;
		my $range = $totalrows - $metadata;
		my $addthis = $metadata -1;
		#save the chosen loci indexes as a hash
		my %used = map {$_ => 42} @indexes;
		my %rndindex = ();
		my $null = 0;
		my $size = 0;
		
		while ($size < $random) {
			#generate random
			my $rndm = int(rand($range)) + $addthis;
			
			if (exists $used{$rndm}) { $null++; }
			elsif (exists $rndindex{$rndm}) { $null++; }
			else { $rndindex{$rndm} = 42; }
			$size = scalar keys %rndindex;
			#print "\nrandom = $rndm, size = $size / $random";
		}
		print "$size random positions for loci generated.\n";
		
		my @chosen = keys %rndindex;
		
		#save chosen lines in the file
		foreach (@chosen) { push (@randomloci, $wholefile[$_]); }
		
		#save random loci names
		my @rndlist = ();
		foreach my $tosave (@randomloci) {
			my @saveline = split ("\t", $tosave);
			my $scaffold = $saveline[0];
			$scaffold =~ s/\D//g;
			my $ref = $saveline[2];
			#print "\n$ref\n";
			my @refnum = split(':', $ref);
			#my $new_name = "$scaffold" . "_" . "$refnum[0]";
			
			if ($sepsymbol eq "0") {
			if ($ref=~ /[0-9]*:[0-9]*/) { $sepsymbol = ":"; } else { my $tempsymb = $ref; $tempsymb =~ s/[0-9]*(\D)*[0-9]*/$1/; $sepsymbol = $tempsymb; }
			
			}
			
			my $new_name = "$scaffold$sepsymbol$refnum[0]";
			push (@rndlist, $new_name);
		}
		
		
		my $randomfile = "random$size" . "loci_$filename";
		my $randomname = $filename;
		$randomname =~ s/^(.*?).vcf/$1/;
		$randomname =~ s/^(.*?).VCF/$1/;
		my $randomlist = "random$size" . "names_$randomname";
		
		
		my @subsetrandom = (@headers, @randomloci);
		#save
		my $rndpath = "$subdirpath" . "/" . "$randomfile";
		my $rndlist = "$subdirpath" . "/" . "$randomlist";
		print "Saving as $randomfile... ";
		open my $AZAR, '>', $rndpath or die "\nUnable to create or save \"$rndpath\": $!\n";
		open my $RNDM, '>', $rndlist or die "\nUnable to create or save \"$rndlist\": $!\n";
		foreach (@subsetrandom) {print $AZAR "$_\n";} # Print each entry in our array to the file
		foreach (@rndlist) {print $RNDM "$_\n";} # Print each entry in our array to the file
		close $AZAR; 
		close $RNDM; 
		
	}
	print "done!\n\n";
}
elsif ($singlefile eq "no") {

	#variables saved refer to directories, not files
	my $inputdir = $filename;
	my $inputpath = $filepath;
	my $parentdir = $subdirpath;
	
	
	#read files
	print "Reading files at $inputdir  -->  ";
	opendir(DIR, $inputpath);						#open the directory 
	my @infiles = readdir(DIR);					#extract filenames
	closedir(DIR);
	
	#filter files
	my @vcffiles = ();
	foreach my $infile (@infiles) {
		next if ($infile =~ /^\.$/);				#don't use any hidden file
		next if ($infile =~ /^\.\.$/);			
		if ($infile =~ /.*?\.vcf$/) { push (@vcffiles, $infile); }		#save the bayescan files
		elsif ($infile =~ /.*?\.VCF$/) { push (@vcffiles, $infile); }		#save the bayescan files
	}
	my $filenumber = scalar @vcffiles;
	
	if ($filenumber == 0) { die "\n\n$filenumber files found matching \'.vcf'\ or \'.VCF'\, program will exit.\nCheck the help information if needed: $version -help\n\n"; }
	else { print "$filenumber files found. "; }
	
	my $outpath = "no default";
	$inputname =~ s/(.*?)\/$/$1/;  		#clean "/" from the end of the name
	$outdir =~ s/(.*?)\/$/$1/;  		#clean "/" from the end of the name
	my @passedir = split ('/', $outdir);
	my $dirlength = scalar @passedir;
	if ($dirlength > 1) {
		$outpath = $outdir;
		#print "\nlong path\n";
		#print "\n$outdir\n";
	}
	elsif ($dirlength <= 1) {
		$outpath = "$subdirpath/$inputname/$outdir";
		#print "\nshortpath\n";
		#print "\n\"$outpath\"\n";
	}
	else { die "\n\nWrong shit here\n\n"; }
	unless(-e $outpath or mkdir $outpath) {die "Unable to create the directory \"$outpath\"\nMay be you don't have the rights: $!\n"; }
	print "Output files will be saved at\n$outpath\n\n";

	
	
	my $filedone = 0;
	foreach my $onefile (@vcffiles) {
		$filename = "$onefile"; 
	
		# SUBSETING
		if ($filedone == 0) { print "Processing $onefile... \n" }
		$subdirpath = "$inputpath"; 
		$filepath = "$subdirpath/$filename"; 
		
		open my $VCFFILE, '<', $filepath or die "\nUnable to find or open $filepath: $!\n";



		# READ VCF FILE
		my @wholefile = ();

		while (<$VCFFILE>) {
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			my $row = $_;  		#save line
			$row =~ s/\s+$//;  		#clean white tails in lines
			push (@wholefile, $row);
		}
			
		close $VCFFILE;

		
		# EXTRACT LOCI FROM LIST


		my @headers = ();
		my @indexes = ();
		my @namelist = ();
		my $savedloci = 0;
		my $printloci = " ";
		# FIRST LETS GET THE LOCI LIST
		if ($locilist ne "no default") {
			my @subsetlisted = ();
			my @keptlisted = ();
			foreach (@wholefile) {
				my $index=0;
				my $line = $_;
				my @wholeline= split('\t', $line);  		#split columns as different elements of an array
				#save the headers
				if ($wholeline[0]=~ /^#.*?/) {
					push (@headers, $line);  		#save the first rows with metadata (they all start with "#")
					#print "Saving metadata $wholeline[0]\n";
				}
				# save new loci names listed
				elsif ($sample =~ /^[0-9]*_[0-9]*$/) {
					# generate the "new locinames" from vcf
					my $scaffold = $wholeline[0];
					$scaffold =~ s/\D//g;
					my $ref = $wholeline[2];
					my @refnum = split(':', $ref);
					my $new_name = "$scaffold" . "_" . "$refnum[0]";
					# if they exist in the list keep them
					if (exists $goodloci{$new_name}) { 
						push (@keptlisted, $line);
						push (@indexes, $index);
						push (@namelist, $new_name);
					}
				}
				elsif ($sample=~ /^[0-9]*:[0-9]*/) {
					my $ref = $wholeline[2];
					#if the loci exist in the list, keep them
					if (grep {/^$ref/} keys %goodloci) { 
						push (@keptlisted, $line);
						push (@indexes, $index);
						push (@indexes, $index);
						push (@namelist, $ref);
					}
				}
				$index++;
			}
			$savedloci = scalar @keptlisted;
			if ($filedone == 0) { print "$savedloci loci matching the list at $locifile were found.\n"; }
			
			my $selectedfile = "candidate$savedloci" . "loci_$filename";
			
			@subsetlisted = (@headers, @keptlisted);
			#save
			my $listpath = "$outpath" . "/$selectedfile";
			if ($filedone == 0) { print "Saving as $selectedfile... "; }
			open my $LIST, '>', $listpath or die "\nUnable to create or save \"$listpath\": $!\n";
			foreach (@subsetlisted) {print $LIST "$_\n";} # Print each entry in our array to the file
			$printloci = "File with $savedloci candidate loci saved.";
			close $LIST; 
			
			my $noext = $filename;
			$noext =~ s/^(.*?).vcf/$1/;
			$noext =~ s/^(.*?).VCF/$1/;
			my $coincident = "$outpath" . "/" . "$noext" . "_$savedloci" . "locinames";
			open my $NAMES, '>', $coincident or die "\nUnable to create or save \"$coincident\": $!\n";
			foreach (@namelist) {print $NAMES "$_\n";} # Print each entry in our array to the file
			close $NAMES; 
			
			if ($filedone == 0) { print "done!\n\n"; }
		}
		
		my $randomnum = 0;
		my $printrandom = " ";
		if ($random != 0) {
			my @randomloci = ();
			if ($filedone == 0) { print "Now generating random SNP positions... "; }
			
			# extract as many loci as loci there are in the list
			if ($random == 1) { $randomnum = $savedloci; }
			my $metadata = scalar @headers;
			my $totalrows = scalar @wholefile;
			my $range = $totalrows - $metadata;
			my $addthis = $metadata -1;
			#save te chose loci indexes as a hash
			my %used = map {$_ => 42} @indexes;
			my %rndindex = ();
			my $null = 0;
			my $size = 0;
			
			while ($size < $randomnum) {
				#generate random
				my $rndm = int(rand($range)) + $addthis;
				
				if (exists $used{$rndm}) { $null++;; }
				elsif (exists $rndindex{$rndm}) { $null++; }
				else { $rndindex{$rndm} = 42; }
				$size = scalar keys %rndindex;
				#print "\nrandom = $rndm, size = $size / $random";
			}
			if ($filedone == 0) { print "$size random positions for loci generated.\n"; }
			
			my @chosen = keys %rndindex;
			
			foreach (@chosen) { push (@randomloci, $wholefile[$_]); }
			
					#save random loci names
			my @rndlist = ();
			foreach my $tosave (@randomloci) {
				my @saveline = split ("\t", $tosave);
				my $scaffold = $saveline[0];
				$scaffold =~ s/\D//g;
				my $ref = $saveline[2];
				my @refnum = split(':', $ref);
				my $new_name = "$scaffold" . "_" . "$refnum[0]";
				push (@rndlist, $new_name);
			}
			
			
			
			my $randomfile = "random$size" . "loci_$filename";
			
			my $randomname = $filename;
			$randomname =~ s/^(.*?).vcf/$1/;
			$randomname =~ s/^(.*?).VCF/$1/;
			my $randomlist = "random$size" . "names_$randomname";
			
			my @subsetrandom = (@headers, @randomloci);
			#save
			my $rndpath = "$outpath/" . "$randomfile";
			if ($filedone == 0) { print "Saving as $randomfile... "; }
			open my $AZAR, '>', $rndpath or die "\nUnable to create or save \"$rndpath\": $!\n";
			foreach (@subsetrandom) {print $AZAR "$_\n";} # Print each entry in our array to the file
			close $AZAR;
			#save list of loci
			my $rndlist = "$outpath" . "/" . "$randomlist";
			open my $RNDM, '>', $rndlist or die "\nUnable to create or save \"$rndlist\": $!\n";
			foreach (@rndlist) {print $RNDM "$_\n";} # Print each entry in our array to the file
			close $RNDM; 
			$printrandom = "File with $randomnum random loci saved.";
		}
		if ($filedone == 0) { print " done!\n\n"; }
		$filedone++;
		print "File $filedone of $filenumber done ($filename). $printloci $printrandom \n";
	}
	print "All files done!\n\n";
}
print "$version finished.\n\n";

