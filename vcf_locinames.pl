#!/usr/bin/perl
use strict ; use warnings;

# names_by_position.pl   			# by M'Óscar 
my $version = "vcf_locinames 1.0";

#extract a list of locinames from a vcf file


##############changelog########################
####   07/12/2020       v 1.0         #########
### added option for other lociname formats ###
###############################################


### ARGUMENTS ALL CAN BE INPUT FROM COMMAND LINE ####
#####################################################
#####################################################


my $vcffile = "no default";		  #vcf file
#my $vcffile = "no default";		  #vcf file
my $dir = "no default"; 	#working dir
#my $dir = "no default"; 	#working dir
my $scaffold = "id";		#to use default loci names
#my $scaffold = "locus";		#to use scaffold_locus position
#my $scaffold = "snp";		#to use scaffold_SNP position

################################################
################################################


#Help!
my %arguments = map { $_ => 1 } @ARGV;
if(exists($arguments{"help"}) || exists($arguments{"--help"}) || exists($arguments{"-help"}) || exists($arguments{"-h"}) || exists($arguments{"--h"})) {
	die "\n\n\t   $version   Help Information\n\t-------------------------------------------------\n
	Use this program will extract the column with locus IDs from vcf files, withouth the final +/- symbols
	Just that, that's it.\n
	--vcf                   vcf file. If none parsed will read all vcf file in the directory.
	--dir                   [if no --vcf] set directory to read the vcf files from. Default: working directory
	--name                  Three modes: \'--name id\' Use SNP ID (third column). This is the dafault option.
	                        \'--name locus\' Use \"scaffold number\"_\"locus position in scaffold\" (2nd column).
	                        \'--name snp\' Use \"scaffold number\"_\"SNP position in locus\" (first part of SNP ID).\n\n"
}


################ PASSING ARGUMENTS


use Getopt::Long;

GetOptions( "dir=s" => \$dir,        #   --names
            "name=s" => \$scaffold,       #   --scaffold
            "vcf=s" => \$vcffile );       #   --output



############### DIRECTORY PATH

use Cwd qw(cwd);
my $localdir = cwd;



if ($dir eq "no default") { $dir = $localdir; }



if ($vcffile eq "no default") {
	opendir(DIR, $dir);						#open the directory with the files
	my @files = readdir(DIR);					#extract filenames
	closedir(DIR);

	# parse all files
	my @infiles = ();
	foreach my $file (@files) {
		next if ($file =~ /^\.$/);				#don't use any hidden file
		next if ($file =~ /^\.\.$/);			
		if ($file =~ /.*\.vcf$/) { push (@infiles, $file); }		#save only the right files
	}

	my $filenum = scalar @infiles;

	print "$version found $filenum \"vcf\" files found\n\n";
	
	print "processing loci IDs:\n";
	foreach my $infile (@infiles) {
		print "$infile... ";
		
		my @locilist = ();
		
		open my $VCFFILE, '<', $infile or die "\nUnable to find or open $infile: $!\n";
		
		while (<$VCFFILE>) {
			chomp;	#clean "end of line" symbols
			next if /^$/;  		#skip if blank
			next if /^\s*$/;  		#skip if only empty spaces
			next if /^#.*/;  		#skip commented lines
			
			my $line = $_;  		#save line
			$line =~ s/\s+$//;  		#clean white tails in lines
			
			my @wholeline= split('\t', $line);  		#split columns as different elements of an array
			
			my $fixedid = "nodefault";
			
			if ($scaffold eq "id") {
				my $rawid = $wholeline[2];
				$fixedid = $rawid;
				$fixedid =~ s/^([0-9]*:[0-9]*):.$/$1/;  		#delete the last part of the locus ID
			} elsif ($scaffold eq "locus") {
				my $scaffoldstring = $wholeline[0];
				my $scaffoldnum = $scaffoldstring;
				$scaffoldnum =~ s/^Scaffold([0-9]*)$/$1/;  		#delete the word scaffold
				my $posloc = $wholeline[1];
				$fixedid = "$scaffoldnum" . "_" . "$posloc";
			} elsif ($scaffold eq "snp") {
				my $scaffoldstring = $wholeline[0];
				my $scaffoldnum = $scaffoldstring;
				$scaffoldnum =~ s/^Scaffold([0-9]*)$/$1/;  		#delete the word scaffold
				my $rawid = $wholeline[2];
				my $cleanid = $rawid;
				$cleanid =~ s/^([0-9]*):.*$/$1/;  		#delete the second part of ID
				$fixedid = "$scaffoldnum" . "_" . "$cleanid";
			} else { die "\n\n\n\tERROR!\n\tArgument scaffold must be one of this options: id, snp, locus; check help information" }
			
			
			
			
			
			
			
			
			push (@locilist, $fixedid);
		
		}
		
		close $VCFFILE;
		
		#name the output file
		my $locifile = $infile;
		$locifile =~ s/(.*)\.vcf$/locilist_$1/;
		$locifile = "$locifile" . "_$scaffold";
		
		#save list of loci
		open my $LIST, '>', $locifile or die "\nUnable to create or save \"$locifile\": $!\n";
		foreach (@locilist) {print $LIST "$_\n";} # Print each entry in our array to the file
		close $LIST; 
		
		my $locinum = scalar @locilist;
		
		print "$locinum loci names saved\n";
	}
}
else {
	
	print "Reading $vcffile\n";
	
	my @locilist = ();
	
	open my $VCFFILE, '<', $vcffile or die "\nUnable to find or open $vcffile: $!\n";
	
	while (<$VCFFILE>) {
		chomp;	#clean "end of line" symbols
		next if /^$/;  		#skip if blank
		next if /^\s*$/;  		#skip if only empty spaces
		next if /^#.*/;  		#skip commented lines
		
		my $line = $_;  		#save line
		$line =~ s/\s+$//;  		#clean white tails in lines
		
		my @wholeline= split('\t', $line);  		#split columns as different elements of an array
		
		my $fixedid = "nodefault";
		
		if ($scaffold eq "id") {
			my $rawid = $wholeline[2];
			$fixedid = $rawid;
			$fixedid =~ s/^([0-9]*:[0-9]*):.$/$1/;  		#delete the last part of the locus ID
		} elsif ($scaffold eq "locus") {
			my $scaffoldstring = $wholeline[0];
			my $scaffoldnum = $scaffoldstring;
			$scaffoldnum =~ s/^Scaffold([0-9]*)$/$1/;  		#delete the word scaffold
			my $posloc = $wholeline[1];
			$fixedid = "$scaffoldnum" . "_" . "$posloc";
		} elsif ($scaffold eq "snp") {
			my $scaffoldstring = $wholeline[0];
			my $scaffoldnum = $scaffoldstring;
			$scaffoldnum =~ s/^Scaffold([0-9]*)$/$1/;  		#delete the word scaffold
			my $rawid = $wholeline[2];
			my $cleanid = $rawid;
			$cleanid =~ s/^([0-9]*):.*$/$1/;  		#delete the second part of ID
			$fixedid = "$scaffoldnum" . "_" . "$cleanid";
		} else { die "\n\n\n\tERROR!\n\tArgument scaffold must be one of this options: id, snp, locus; check help information" }
		
		push (@locilist, $fixedid);
	
	}
	
	close $VCFFILE;
	
	#name the output file
	my $locifile = $vcffile;
	$locifile =~ s/(.*)\.vcf$/locilist_$1/;
	$locifile = "$locifile" . "_$scaffold";
	#save list of loci
	open my $LIST, '>', $locifile or die "\nUnable to create or save \"$locifile\": $!\n";
	foreach (@locilist) {print $LIST "$_\n";} # Print each entry in our array to the file
	close $LIST; 
	
	my $locinum = scalar @locilist;
	
	 print "$locinum loci names saved\n";
	
	
}





print "\n$version is done!\n\n";



