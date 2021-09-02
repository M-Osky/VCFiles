#!/usr/bin/perl
use strict ; use warnings;

# vcf_joiner.pl   			# by M'Óscar 
my $version = "vcf_joiner_v1.2.pl";

############################

# Use this script to join different vcf files in one
# For options, usage and other information check the help typing the name of the program version and "help" or "--h" or so...
# vcf_pop_subseter.pl -help




#######################   PARAMETERS   #########################
# All of them can and should be set from the command line, check the help.

my $def = "default";  # better do not chage this.

my $inputdir = $def;  		# directory with the vcf files to join.
#my $inputdir = ".";  		# directory with the vcf files to join.

my $filter = ".vcf"; #ending all the files must have in common
#my $filter = "onepop.vcf"; #ending all the files must have in common

my $infocols = 9; # columns of the vcf file before the first allele.

my $outname = $def;

#################################################################################################################




#Help!
my %arguments = map { $_ => 1 } @ARGV;
if(exists($arguments{"help"}) || exists($arguments{"--help"}) || exists($arguments{"-help"}) || exists($arguments{"-h"}) || exists($arguments{"--h"})) {
	die "\n\n\t   $version   Help Information\n\t-------------------------------------------------\n
	This program will join a group of vcf files from the same directory.
	Mostly just call the program from the directory with the vcf files to join, and it should work.
	\n\tCommand line arguments and defaults:
	--dir                     Name (or path) of the directory with the files. By default the local working diretory.
	--filter                  Something all filenames to read must end in in order to be processed. Default: $filter
	--infocols                Number of locus information columns before the first sample. Default: $infocols
	--out                     Name for the output file. Otherwise will generate one from the name of the vcf files opened.\n
	Command line call example:\n\tvcf_joiner.pl --input /home/refmap/out/perpop/ --outname all_pops_reunited.vcf\n\n\n";
}



################ PASSING ARGUMENTS


use Getopt::Long;

GetOptions( "dir=s" => \$inputdir,    #   --dir
            "filter=s" => \$filter,      #   --filter
            "out=s" => \$outname,      #   --out
            "infocols=i" => \$infocols );   #   --infocols



#read files

use Cwd qw(cwd);
my $localdir = cwd;

my $dirname = "no default";

if ($inputdir ne $def) { $dirname = $inputdir; }  else { $dirname = $localdir; }

print "\n\n$version is reading files from $dirname\n";

opendir(DIR, $dirname);
my @files = readdir(DIR);
closedir(DIR);



# process the files one by one

my $f = 0;
#my @metadata = ();
#my @samplealleles = ();
my @chosenfiles = ();
my @firstrows = ();
#my $lastinfo = $infoloci -1;
#my $firstsamp = $infoloci;
my $reffile = "\"no file name was read\"";
my $headers = "none";
my %ref_listloci = ();
my @deleted_first = ("List of loci deleted per file", "--------------------------", "");
my @deleted_loci = ();
my %locilist = ();
my $snpnum = 0;
my $shit = 0;
my $common = 0;


foreach my $file (@files) {
	next unless ($file =~ m/^.*$filter$/);
	push (@chosenfiles, $file);
	# create output name if needed
	if ($f == 0 && $outname eq $def) { my $rawname = $file; $rawname =~ s/^.+?_(.*?)\.vcf/$1/; $outname = "$rawname" . "_joined.vcf"; }
	#save the info common in all the files.
	if ($f == 0) {
		$reffile = $file;
		#open
		open my $VCFFILE, '<', $file or die "\nUnable to open $file: $!\n";
		my $k = 0;
		my $head = 0;
		while (<$VCFFILE>) {
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			
			my $line = $_;  		#save line
			$line =~ s/\s+$//;  		#clean white tails in lines
			if ($line=~ /^##.*?/) {
				push (@firstrows, $line); 	#save the first rows with metadata (they all start with "##")
			}
			elsif ($head == 0) {
				$headers = $line;
				$head = 1;
			}
			else {
				my @wholeline= split('\t', $line);  		#split columns as different elements of an array
				my $snpid = "$wholeline[0]" . "_$wholeline[1]" . "_$wholeline[2]";  #save the loci info
				$ref_listloci{$snpid} = $line;
				$k++;
			}
		}
		close $VCFFILE;
		$snpnum = $k;
		print "$k SNPs from $file saved as reference.\n\n";
	}
	elsif ($f > 0) {
		#open
		open my $VCFFILE, '<', $file or die "\nUnable to open $file: $!\n";
		my $head = 0;
		my $k = 0;
		my $deleted=0;
		my $kept = 0;
		my %new_locilist = ();
		
		while (<$VCFFILE>) {
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			
			next if /^##.*?/;  		#skip first lines with comments 
			
			my $line = $_;  		#save line
			$line =~ s/\s+$//;  		#clean white tails in lines
			if ($head == 0) {
				my @wholeline= split('\t', $line);  		#split columns as different elements of an array
				my $colnum = scalar @wholeline;
				my $lastindex = $colnum - 1;
				my @alleles = @wholeline[$infocols..$lastindex];
				my $newline = join("\t", @alleles);
				$headers = "$headers" . "\t" . "$newline";
				$head = 1;
			}
			else {
				my @wholeline= split('\t', $line);  		#split columns as different elements of an array
				my $newlocus = "$wholeline[0]" . "_$wholeline[1]" . "_$wholeline[2]";  #save the loci info
				$new_locilist{$newlocus} = 42;  #save all loci from the new file
				
				#check if the new loci exist in the old loci list
				if (exists $ref_listloci{$newlocus}) {
					my $colnum = scalar @wholeline;
					my $lastindex = $colnum - 1;
					my @alleles = @wholeline[$infocols..$lastindex];
					my $newline = join("\t", @alleles);
					$ref_listloci{$newlocus} = "$ref_listloci{$newlocus}" . "\t" . "$newline";
					$kept++;
					$k++;
				}
				else {
					my $rejected = "$file    $wholeline[0] $wholeline[1] $wholeline[2]";
					push (@deleted_loci, $rejected);
					$deleted++;
					$k++;
				}
			}
		}
		close $VCFFILE;
		
		$common=0;
		my @savedsnps = keys %ref_listloci;
		
		#now check how many loci from the reference hash (with all the loci processed till now) are not in the new opened file
		foreach my $snp (@savedsnps) {
			if (exists $new_locilist{$snp}) { 
				$common++;
			}
			else {
				delete $ref_listloci{$snp};
				my @wholeline= split('_', $snp);  		
				my $rejected = "$reffile    $wholeline[0] $wholeline[1] $wholeline[2]";
				push (@deleted_first, $rejected);
				$shit++;
			}
		}
		if ($f ==1 ) { print "$file done: From $k SNPs, $kept were saved, $deleted not found in the previous file.\n"; }
		else { print "$file done: From $k SNPs, $kept were saved, $deleted not found in previous files.\n"; }
		
	}
	$f++;
}

print "$reffile done: From $snpnum SNPs, $common were saved, $shit not found in the rest of files.\n";


#save a log file with all the deleted loci
my @alldeleted = (@deleted_first, @deleted_loci);
my $numdeleted = scalar @alldeleted;
my $programname = $version;
$programname =~ s/^(.*)\..+?$/$1/;
my $logfile = "$programname" . ".log";
open my $LOG, '>', $logfile or die "\nUnable to create or save \"$logfile\": $!\n";

if ($numdeleted > 0) { foreach my $row (@alldeleted) { print $LOG "$row\n"; } } # print each entry to file 
else { print $LOG "\n None deleted!\n"; }

close $LOG;

#now save the information;
push (@firstrows, $headers);
my @list_loci = keys %ref_listloci;

# sort back the SNPs
#first sort by position (second column)
my @sorted_list1 = sort { ($a =~ /^[a-zA-Z0-9]*_([0-9]*).*$/)[0] <=> ($b =~ /^[a-zA-Z0-9]*_([0-9]*).*$/)[0] } @list_loci;

#then sort by scaffold
my @sorted_list2 = sort { ($a =~ /^[a-zA-Z]*([0-9]*).*$/)[0] <=> ($b =~ /^[a-zA-Z]*([0-9]*).*$/)[0] } @sorted_list1; ##this one works fine, but just for scaffolds to use with another by possition

my $finalnum = scalar @sorted_list2;
print "\nAll $f files done! Saving $finalnum SNPs as $outname\n";


open my $SAVE, '>', $outname or die "\nUnable to create or save \"$outname\": $!\n";
foreach my $row (@firstrows) { print $SAVE "$row\n"; } # print each entry to file
foreach my $row (@sorted_list2) { print $SAVE "$ref_listloci{$row}\n"; } # print each entry to file
close $SAVE; 

print "Done!\n\n";


