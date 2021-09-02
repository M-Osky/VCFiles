# VCFiles
Variant Call Format files are difficult to understand at first, I already had to dive into them in order to learn how to edit them, so you can use my scripts. Some of these functions already exist in vcftools, but it seems it stopped being updated.
so you can use these scripts to do various tasks with them

All scripts have their own HELP information, just call them with the usual flags: -h --help 

# vcf_pop_subseter
outputs a new vcf file with only samples from the specified populations

# vcf_loci_subseter
outputs a new vcf file with only the specified loci or with a number n of random loci

# vcf_locinames
extracts the locus IDs from a vcf file

# vcf_random_resampler
outputs a new vcf file with a number n of samples per population, useful if there are differences in your sample size per population, samples per population are pulled randomly.

# coverage_vcfilter
this will analyse mean coverage/depth (DP) per locus. 
Contrary to vcftools you can choose how to handle missing depth information: genotype with missing DP can ither be ignored when computing the average or otherwise they can be considered as DP=0. Caution: vcftools seem to use both this criteria, for the table generated it ignores missing, but when filtering out loci according to coverage it seems to consider missing as zero.

# vcf_joiner
joins multiple vcf files in one file. Careful if they have different loci.

