UCE BioInformatics Pipeline for Phylogenetic Studies

by Carl H. Oliveros

Based on:

https://github.com/faircloth-lab/phyluce/blob/master/docs/pre-processing-qc.rst

https://github.com/faircloth-lab/phyluce/blob/master/docs/pre-processing-assembly.rst

https://github.com/faircloth-lab/phyluce/blob/master/docs/uce-processing.rst

Read these pages in conjunction with this document.

Outline
-------
1 - Run illumiprocessor.py
2 - Run Velvet to Assemble Contigs for each species
3 - Match Probes to Themselves
4 - Extract UCE contigs 
5 - Inspect data using sqlite
6 - Assemble dataset
7 - Sequence alignment
8 - Format data for phylogenetic analysis


STEP 1 - Run illumiprocessor.py
-------------------------------

Create a configuration file for illumiprocessor.py.  This file will tell 
illumiprocessor.py which index from the sequence reads corresponds to which 
individual (combos section), the sequences corresponding to each index,
what adapter sequences was used for the Illumina run (adapters section) 
so it can remove fragments of those (with the appropriate indices inserted 
in the asterisks), the input file names (params section), and what name to
use in the output files (remap section). USE ALL SMALL LETTERS FOR NAMES.

Open a text editor and create the mapping file as follows:

[adapters]
N7:CTGTCTCTTATACACATCTCCGAGCCCACGAGAC*ATCTCGTATGCCGTCTTCTGCTTG
N5:CTGTCTCTTATACACATCTGACGCTGCCGACGA*GTGTAGATCTCGGTGGTCGCCGTATCATT

[indexes]
N701:TAAGGCGA
N702:CGTACTAG
N703:AGGCAGAA
N704:TCCTGAGC
N705:GGACTCCT
N706:TAGGCATG
N707:CTCTCTAC
N708:CAGAGAGG
N709:GCTACGCT
N710:CGAGGCTG
N711:AAGAGGCA
N712:GTAGAGGA
N501:TAGATCGC
N502:CTCTCTAT
N503:TATCCTCT
N504:AGAGTAGA
N505:GTAAGGAG
N506:ACTGCATA
N507:AAGGAGTA
N508:CTAAGCCT

[params]
separate reads:True
read1:{name}_L001_R1_001.fastq.gz
read2:{name}_L001_R2_001.fastq.gz

[combos]
acryllium-vulturinum:N503,N707

[remap]
acryllium_vulturinum_CTCTCTAC-TATCCTCT:acryllium-vulturinum


Save the file on the /public/uce/work folder as pre-process.conf.

The illumiprocessor python script trims adapter contamination and performs 
quality checking of reads. Create a folder in your work folder and copy all raw
read files in that folder.

In the command line, call the program and insert arguments for the input (the 
folder containing raw reads), the output 
directory, your config file, and the number of processors to be used. 

Navigate to your working directory.

XX mkdir cleaned-reads

version 2
XX illumiprocessor --input L001/ --output cleaned-reads/ --config illumiprocessor.conf --cores 12

version 1
XX illumiprocessor.py /public/uce/work/L006/ /public/uce/work/cleaned-reads/ /public/uce/work/preprocess.conf --remap --complex --cores 12

If CBOT was not used in the Illumina sequencing run, you will have to run 
illumiprocessor.py twice (assuming a rapid run of two lanes).  You will 
have to make two versions of the config file differing only in the file name 
extensions given in the params settings.  Store the output of one run in 
cleaned-reads, and the other in cleaned-reads2.  The outputs of both 
illumiprocessor runs should be concatenated.

$ cd /public/uce/work
$ mkdir cleaned-reads2
$ illumiprocessor.py /public/uce/work/L007/ /public/uce/work/cleaned-reads/ /public/uce/work/preprocess1.conf --remap --complex --cores 12
$ illumiprocessor.py /public/uce/work/L008/ /public/uce/work/cleaned-reads2/ /public/uce/work/preprocess2.conf --remap --complex --cores 12
$ for iter in * ; do cat /public/uce/work/cleaned-reads2/$iter/interleaved-adapter-quality-trimmed/$iter-READ1and2-interleaved.fastq.gz >> /public/uce/work/cleaned-reads/$iter/interleaved-adapter-quality-trimmed/$iter-READ1and2-interleaved.fastq.gz; cat /public/uce/work/cleaned-reads2/$iter/interleaved-adapter-quality-trimmed/$iter-READ-singleton.fastq.gz >> /public/uce/work/cleaned-reads/$iter/interleaved-adapter-quality-trimmed/$iter-READ-singleton.fastq.gz; done
$ for iter in * ; do cat /public/uce/work/cleaned-reads2/$iter/split-adapter-quality-trimmed/$iter-READ1.fastq.gz >> /public/uce/work/cleaned-reads/$iter/split-adapter-quality-trimmed/$iter-READ1.fastq.gz; cat /public/uce/work/cleaned-reads2/$iter/split-adapter-quality-trimmed/$iter-READ2.fastq.gz >> /public/uce/work/cleaned-reads/$iter/split-adapter-quality-trimmed/$iter-READ2.fastq.gz; cat /public/uce/work/cleaned-reads2/$iter/split-adapter-quality-trimmed/$iter-READ-singleton.fastq.gz >> /public/uce/work/cleaned-reads/$iter/split-adapter-quality-trimmed/$iter-READ-singleton.fastq.gz; done

XX  $ for iter in * ; do cat /public/uce/aug/cleaned-reads8/$iter/interleaved-adapter-quality-trimmed/$iter-READ-singleton.fastq.gz >> /public/uce/aug/cleaned-reads/$iter/interleaved-adapter-quality-trimmed/$iter-READ-singleton.fastq.gz; cat /public/uce/aug/cleaned-reads8/$iter/split-adapter-quality-trimmed/$iter-READ1.fastq.gz >> /public/uce/aug/cleaned-reads/$iter/split-adapter-quality-trimmed/$iter-READ1.fastq.gz; cat /public/uce/aug/cleaned-reads8/$iter/split-adapter-quality-trimmed/$iter-READ2.fastq.gz >> /public/uce/aug/cleaned-reads/$iter/split-adapter-quality-trimmed/$iter-READ2.fastq.gz; done


STEP 2a - Run Velvet to Assemble Contigs for each species
---------------------------------------------------------

For detailed instructions on this step, visit:
https://github.com/faircloth-lab/phyluce/blob/master/docs/pre-processing-assembly.rst

Move into the cleaned-reads folder.
	
$ cd /public/uce/work/cleaned-reads/

To run Velvet Optimiser (a wrapper script that determines the best settings for 
Velvet given your dataset), use the following commands. First, print the names 
of the cleaned files you will be using.

$ VelvetOptimiser.pl --s 69 --e 75 --k=n50 --c=tbp -t 6 -a -f "-fastq.gz -shortPaired $(echo $( ls */split-adapter-quality-trimmed/*-READ?.fastq.gz ) ) -short $(echo $( ls */split-adapter-quality-trimmed/*-READ-singleton.fastq.gz ) )"

Once you have an idea of the kmer range that will work with your data, run 
assemblo.py, which will use velvet to assemble contigs for each species.

$ cd /public/uce/work/
$ assemblo.py ./cleaned-reads 71 83 6
[problematic = run with python and full path name of script]

In the command above specify the kmer range and the number of processors you 
want to use.  The resulting contigs (symlinks) will be found in the contigs 
directory in the cleaned reads directory.


STEP 2b - Run Trinity to Assemble Contigs for each species
----------------------------------------------------------

You will run the python script assemblo_trinity.py to assemble contigs for
each species.  You need to create a run config file, that tells the script
where the read data are and what you want the final data named. Below is an 
example config file that was named "trinity-assemblies.conf":

[samples]
aegithina-lafresnayei-23213:/public/uce/work/cleaned-reads/aegithina-lafresnayei-23213
aleadryas-rufinucha-16569:/public/uce/work/cleaned-reads/aleadryas-rufinucha-16569
artamus-cinereus-6183:/public/uce/work/cleaned-reads/artamus-cinereus-6183

The name you want things called on output is on the left side. Within each of 
the paths on the right side, there is a folder called "split-adapter-quality-
trimmed" (for PE data) or "interleaved-adapter-quality-trimmed" where the reads
should be located.  

Invoke assemblo_trinity.py giving the name of the config file, the name of
an output path, the name of the subfolder where reads are located, and the 
number of cores you want the program to run.

$ mkdir trinity-assemblies

XX python ~/phyluce/bin/assembly/assemblo_trinity.py --conf trinity-assemblies.conf --output trinity-assemblies --subfolder 'split-adapter-quality-trimmed' --cores 12 [--clean]

This should assemble everything in it's own folder in "trinity-assemblies".  
There will also be a "contigs" folder in that directory that contains 
symlinks to the contigs that are assembled, named appropriately for the 
"match_contigs_to_probes.py" script.  The --clean option is optional and is
used to clean up unnecessary trinity files that take up a lot of disk space.


STEP 3 - Match Probes to Themselves
-----------------------------------

This step can be skipped if previously done.  Just copy the probe and dupe files
to your working directory.

The probes file (LSU-Custom-Array-Jan-2013.fasta) needs to be present in the 
/public/uce/work directory.

Run easy_lastz:

XX  python ~/phyluce/bin/share/easy_lastz.py --target 5472_probes.fasta --query 5472_probes.fasta --output 5472_probes_to_self.fasta --identity 85


STEP 4 - Extract UCE contigs 
----------------------------

Match assembled contigs to self-matched probes and extract the contigs matching 
UCE probes.  

[OLD VERSION]
$ match_contigs_to_probes.py /public/uce/work/cleaned-reads/contigs /public/uce/work/LSU-Custom-Array-Jan-2013.fasta /public/uce/work/lastz --regex "_p[1-9]+$" --repl "" --dupefile /public/uce/work/LSU-Custom-Array-Jan-2013-to-self.fasta

XX  python ~/phyluce/bin/assembly/match_contigs_to_probes.py --contigs trinity-assemblies/contigs/ --probes 5k-probes.fasta  --dupefile 5k-probes-to-self.fasta --output ./lastz/ --regex "(chr\w+)(?:_probe\d+)"


STEP 5 - Inspect data using sqlite
----------------------------------

See Brant's github page for more details.  It's all there!

https://github.com/faircloth-lab/phyluce/blob/master/docs/uce-processing.rst


STEP 6 - Assemble dataset
-------------------------

Create a dataset configuration in a file, e.g. datasets.conf, and save it
in the working directory.  It should look something like:

[dataset1]
genus_species1
genus_species2
genus_species3
genus_species4
genus_species5

[dataset2]
genus_species4
genus_species5
genus_species6


A. Run a query on the database:

XX  python ~/phyluce/bin/assembly/get_match_counts.py --locus-db lastz/probe.matches.sqlite --taxon-list-config datasets.conf --taxon-group 'trogons' --output trogons.conf --log-path trogons_log

This will produce a config file listing taxon names from dataset1 as well as a 
list of loci for which all taxa in the list has data.  If an incomplete data
matrix is desired, add the flag --incomplete-matrix as below:

XX  python ~/phyluce/bin/assembly/get_match_counts.py --locus-db lastz/probe.matches.sqlite --taxon-list-config datasets.conf --taxon-group 'trogons' --output trogons.inc.conf --log-path trogons_inc_log --incomplete-matrix

B. Generate a fasta file for your dataset.

XX  python ~/phyluce/bin/assembly/get_fastas_from_match_counts.py --contigs trinity-assemblies/contigs/ --locus-db lastz/probe.matches.sqlite --match-count-output trogons.conf --output trogons.fasta --log-path trogons_log

If working with an incomplete data matrix, add the flag --incomplete-matrix and
the path to an incomplete.nostrict file that will hold missing locus information 
as below:

XX  python ~/phyluce/bin/assembly/get_fastas_from_match_counts.py --contigs trinity-assemblies/contigs/ --locus-db lastz/probe.matches.sqlite --match-count-output trogons.inc.conf --output trogons.inc.fasta --incomplete-matrix trogons.inc.nostrict --log-path trogons_inc_log


STEP 7 - Calculate coverage for Trinity-aligned datasets
--------------------------------------------------------

For this step you will need:
- cleaned reads folder that contains a folder for each sample and raw reads 1
and 2 of each sample in their respective folder
- trinity assemblies folder that contains a folder for each sample and the 
contigs.fasta file
- a configuration file that maps each sample to the location of its raw reads
files location
- a database of contigs-uce mapping (probe.matches.sqlite)
- the output of get_match_counts.py run for an incomplete matrix (so that
coverage will be calculated for all enriched loci).

The following commands will help you copy the needed files from your cleaned-
reads folder and trinity-assemblies folder to the cluster.

$ mkdir trogons_clean
$ dest=/media/Dicrurus/uce/trogons_clean
$ [cd to cleaned reads folder]
$ for i in apalharpactes_mackloti_b49104 apaloderma_aequatoriale_8461 euptilotis_neoxenus_prs2606 harpactes_ardens_26958 harpactes_erythrocephalus_9970 harpactes_oreskios_23185 pharomachrus_antisianus_b22870 priotelus_roseigaster_6363 priotelus_temnurus_5565 trogon_personatus_gfb2125 trogon_violaceus_rop258 ceyx_argentata_19269 otus_elegans_10975; do mkdir $dest/$i ; mkdir $dest/$i/split-adapter-quality-trimmed; cp $i/split-adapter-quality-trimmed/*READ?.fastq.gz $dest/$i/split-adapter-quality-trimmed; done;

$ mkdir trogons_trinity
$ dest=/media/Dicrurus/uce/trogons_trinity
$ [cd to trinity assemblies folder]
$ for i in apalharpactes_mackloti_b49104 apaloderma_aequatoriale_8461 euptilotis_neoxenus_prs2606 harpactes_ardens_26958 harpactes_erythrocephalus_9970 harpactes_oreskios_23185 pharomachrus_antisianus_b22870 priotelus_roseigaster_6363 priotelus_temnurus_5565 trogon_personatus_gfb2125 trogon_violaceus_rop258 ceyx_argentata_19269 otus_elegans_10975; do mkdir $dest/$i ; cp -L $i/contigs.fasta $dest/$i; done;

Use the following command in your PBS script

A. Calculate coverage for all contigs

XX python /scratch/oliveros/phyluce/bin/assembly/get_trinity_coverage.py --assemblies /scratch/oliveros/trogons/trogons_trinity/ --assemblo-config /scratch/oliveros/trogons/trinity.conf --cores 4 --subfolder 'split-adapter-quality-trimmed' --log-path /scratch/oliveros/trogons/trogons_covlog/

B. Calculate coverage for UCE contigs.  This module is ok to run on a single processor interactive session:
XX python /scratch/oliveros/phyluce/bin/assembly/get_trinity_coverage_for_uce_loci.py --assemblies /scratch/oliveros/trogons/trogons_trinity/ --match-count-output /scratch/oliveros/trogons/trogons.inc.conf --locus-db /scratch/oliveros/trogons/probe.matches.sqlite --output /scratch/oliveros/trogons/trogons_coverage --log-path /scratch/oliveros/trogons/trogons_covlog/

 ** *.jar paths hardcoded in /phyluce/bwa.py
 ** add /scratch/oliveros/phyluce/:/scratch/oliveros/biopython-1.64 to PYTHONPATH

Summarizing output:
cat get_trinity_coverage.log | grep 'Processing\|mean coverage' | sed -r 's/.+Processing (\w+).+/\1/' | sed -r 's/.+INFO -\s+([0-9]+) contigs, mean coverage = ([0-9]+\.[0-9]), mean length = ([0-9]+\.[0-9])/\1\t\2\t\3/' | sed '$!N;s/\n/\t/' 

cat get_trinity_coverage_for_uce_loci.log | grep 'Processing\|mean trimmed coverage' | sed -r 's/.+Processing (\w+).+/\1/' | sed -r 's/.+INFO -\s+([0-9]+) contigs, mean trimmed length = ([0-9]+\.[0-9]), mean trimmed coverage = ([0-9]+\.[0-9])x, on-target bases \(uce contigs\) = ([0-9]+\.[0-9])%, unique reads aligned \(all contigs\) = ([0-9]+\.[0-9])%/\1\t\2\t\3\t\4\t\5/' | sed '$!N;s/\n/\t/'

run time (refine this)
cat get_trinity_coverage.log | grep 'Processing\|Completed' |sed -r 's/(.+) (.+),.+/\1\t\2/'
 
STEP 8 - Sequence alignment
---------------------------

A. Align the dataset.  The following script will create aligned nexus files by
locus.

XX python ~/phyluce/bin/align/seqcap_align_2.py --fasta trogons.fasta --output trogons_aligned --output-format fasta --taxa 71 --aligner mafft --cores 12 --log-path trogons_log

If working with an incomplete matrix:

1.  Perform alignment as above but add the flag --incomplete-matrix:
XX python ~/phyluce/bin/align/seqcap_align_2.py --fasta trogons.inc.fasta --output trogons_inc_aligned --taxa 71 --aligner mafft --cores 12 --log-path trogons_inc_log --incomplete-matrix

2.  You can set a minimum number of taxa for a locus to be included in your 
dataset:
XX python ~/phyluce/bin/align/get_only_loci_with_min_taxa.py --alignments trogons_inc_aligned --taxa 71 --output trogons_inc_min_75percent --percent 0.75 --cores 12 --log-path trogons_inc_log

3.  Insert missing data designators for taxa missing from the alignment of a 
given locus:

XX python ~/phyluce/bin/align/add_missing_data_designators.py --alignments trogons_inc_min_75percent --output trogons_inc_min_75percent_with_missing --output-format fasta --match-count-output trogons.inc.conf --incomplete-matrix trogons.inc.nostrict --cores 12 --log-path trogons_inc_log

B. Trim alignments using GBlocks

XX python ~/phyluce/bin/align/get_gblocks_trimmed_alignments_from_untrimmed.py --alignments trogons_aligned --output trogons_gbtrimmed --b2 0.65 --cores 12 --log-path trogons_log 

If working with an incomplete matrix:

XX python ~/phyluce/bin/align/get_gblocks_trimmed_alignments_from_untrimmed.py --alignments trogons_inc_min_75percent_with_missing --output trogons_inc_min_75percent_gbtrimmed --b2 0.65 --cores 12 --log-path trogons_inc_log 

C. You can get summary stats from your aligned dataset.

XX python ~/phyluce/bin/align/get_align_summary_data.py --alignments trogons_gbtrimmed --input-format nexus --cores 12 --log-path trogons_log

If working with an incomplete matrix:
XX python ~/phyluce/bin/align/get_align_summary_data.py --alignments trogons_inc_min_75percent_gbtrimmed --input-format nexus --cores 12 --log-path trogons_inc_log


D. Screen alignments for bases that are not in the set of IUPAC base codes.

XX  python ~/phyluce/bin/align/screen_alignments_for_problems.py --alignments trogons_gbtrimmed --input-format nexus --output trogons_screened  --cores 12 --log-path trogons_log

****change output to .nex  (changed code so you don't have to do this)

If working with an incomplete matrix:
XX  python ~/phyluce/bin/align/screen_alignments_for_problems.py --alignments trogons_inc_min_75percent_gbtrimmed --input-format nexus --output trogons_inc_min_75percent_screened  --cores 12 --log-path trogons_inc_log

E. After inspecting individual alignments, remove loci names from the taxon names in
the nexus files.

XX python ~/phyluce/bin/align/remove_locus_name_from_nexus_lines.py --taxa 71 --alignment trogons_screened --output trogons_renamed --cores 12 --log-path trogons_log

OR (IF YOU NEED TO CONVERT TO SHORTER NAMES):

XX python ~/phyluce/bin/align/rename_taxa_from_nexus_lines.py --taxa 71 --alignment trogons_screened --output trogons_renamed --cores 12 --log-path trogons_log

If working with an incomplete matrix:
XX python ~/phyluce/bin/align/remove_locus_name_from_nexus_lines.py --taxa 71 --alignment trogons_inc_min_75percent_screened --output trogons_inc_min_75percent_renamed --cores 12 --log-path trogons_inc_log

The dataset is ready for phylogenetic analysis.


STEP 9 - Formatting data for phylogenetic analysis
--------------------------------------------------

For PHYLIP, PhyML, and CloudForest (also MrBayes, since we will need CloudForest's
models info) convert the dataset into strict phylip format.  To keep taxon names 
consistent between analyses, use the strict phylip format for RAxML as well.

XX python ~/phyluce/bin/align/convert_one_align_to_another.py --alignment trogons_renamed/ --output trogons_phylip/ --input-format nexus --output-format phylip --cores 12 --shorten-names --log-path trogons_log

If you do not specify --shorten-names, program will take first 10 characters of name.  Consider using rename_taxa_from_nexus_lines.py

(INCOMPLETE MARTRIX)
XX python ~/phyluce/bin/align/convert_one_align_to_another.py --alignment trogons_inc_min_75percent_renamed/ --output trogons_inc_min_75percent_phylip/ --input-format nexus --output-format phylip --cores 12 --shorten-names --log-path trogons_inc_log


Convert the dataset back to nexus if you wish to have the same short taxon name format across your analyses.
XX python ~/phyluce/bin/align/convert_one_align_to_another.py --alignments trogons_phylip/ --output trogons_nexus/ --input-format phylip --output-format nexus --cores 12 --log-path trogons_log

(INCOMPLETE MARTRIX)
XX python ~/phyluce/bin/align/convert_one_align_to_another.py --alignments trogons_inc_min_75percent_phylip/ --output trogons_inc_min_75percent_nexus/ --input-format phylip --output-format nexus --cores 12 --log-path trogons_inc_log


Format the dataset for RAxML analysis.

XX python ~/phyluce/bin/align/format_nexus_files_for_raxml.py --alignments trogons_renamed/ --output trogons_raxml --log-path trogons_log

(INCOMPLETE MARTRIX)
XX python ~/phyluce/bin/align/format_nexus_files_for_raxml.py --alignments trogons_inc_min_75percent_renamed/ --output trogons_inc_min_75percent_raxml --log-path trogons_inc_log


Estimate gene trees and best fitting substitution models using CloudForest.

$ mkdir trogons_cloudforest
XX python ~/CloudForest/cloudforest/cloudforest_mpi.py trogons_phylip/ trogons_cloudforest/ genetrees /usr/bin/phyml --parallelism multiprocessing --cores 12 2>&1 | tee trogons_cloudforest/genetrees.out

python /tools/cluster/6.2/cloudforest/0.1/scripts/cloudforest_mpi.py trogons_phylip/ trogons_cloudforest/ genetrees /tools/cluster/6.2/cloudforest/0.1/bin/phyml --parallelism multiprocessing --cores 12
(INCOMPLETE MARTRIX)
mkdir trogons_inc_min_75percent_cloudforest
XX python ~/CloudForest/cloudforest/cloudforest_mpi.py trogons_inc_min_75percent_phylip/ trogons_inc_min_75percent_cloudforest/ genetrees /usr/bin/phyml --parallelism multiprocessing --cores 12 2>&1 | tee trogons_inc_min_75percent_cloudforest/genetrees.out

(on the Cluster)
#PBS -N trogons.cf.gt
#PBS -l nodes=1:ppn=20:avx,mem=16000m,walltime=168:00:00
#PBS -M oliveros@ku.edu
#PBS -m abe
#PBS -d /scratch/oliveros/trogons
#PBS -j oe
#PBS -o /scratch/oliveros/trogons/trogons.cf.gt.out
unbuffer python /scratch/oliveros/CloudForest/cloudforest/cloudforest_mpi2.py /scratch/oliveros/trogons/trogons_phylip/ /scratch/oliveros/trogons/trogons_cloudforest/ genetrees /tools/cluster/6.2/cloudforest/0.1/bin/phyml --parallelism multiprocessing --cores 20 > /scratch/oliveros/trogons/trogons.cf.gt.out

(on the Cluster, INCOMPLETE MATRIX)
#PBS -N trogons.cf.gt
#PBS -l nodes=1:ppn=20:avx,mem=16000m,walltime=168:00:00
#PBS -M oliveros@ku.edu
#PBS -m abe
#PBS -d /scratch/oliveros/trogons
#PBS -j oe
#PBS -o /scratch/oliveros/trogons/trogons.cf.gt.out
unbuffer python /scratch/oliveros/CloudForest/cloudforest/cloudforest_mpi2.py /scratch/oliveros/trogons/trogons_inc_min_75percent_phylip/ /scratch/oliveros/trogons/trogons_inc_min_75percent_cloudforest/ genetrees /tools/cluster/6.2/cloudforest/0.1/bin/phyml --parallelism multiprocessing --cores 20 > /scratch/oliveros/trogons/trogons.cf.gt.out


Strip the models from CloudForest output.

XX python ~/phyluce/bin/genetrees/split_models_from_genetrees.py --genetrees trogons_cloudforest/genetrees.tre --output trogons.models.txt

(INCOMPLETE MARTRIX)
XX python ~/phyluce/bin/genetrees/split_models_from_genetrees.py --genetrees trogons_inc_min_75percent_cloudforest/genetrees.tre --output trogons.inc.min.75percent.models.txt

Create a nexus file for MrBayes.

** Can't use dataset1_nexus because output of convert_one_align_to_another.py
is *.nexus whereas format_nexus_files_for_mrbayes.py expects *.nex.  Rename?

XX python ~/phyluce/bin/align/format_nexus_files_for_mrbayes.py --alignments trogons_renamed/ --models trogons.models.txt --output trogons.mrbayes.nex --unlink

(INCOMPLETE MARTRIX)
XX python ~/phyluce/bin/align/format_nexus_files_for_mrbayes.py --alignments trogons_inc_min_75percent_renamed/ --models trogons.inc.min.75percent.models.txt --output trogons.inc.min.75percent.mrbayes.nex --unlink

For a fully partitioned nexus file:
XX python ~/phyluce/bin/align/format_nexus_files_for_mrbayes.py --alignments trogons_nexus/ --models trogons-models.txt --output trogons.fully.mrbayes.nex --unlink --fully-partition

For bootstrapping with Cloudforest:
XX python ~/CloudForest-master/cloudforest/cloudforest_mpi2.py trogons_phylip/ trogons_cloudforest_bootstrap/ bootstraps /usr/bin/phyml --genetrees trogons_cloudforest/genetrees.tre --bootreps 500 --parallelism multiprocessing --cores 12 > /dev/null

Create bootstrapping jobs:
# numbers<-c(040,060,080,140,160,180,240,260,280,340,360,380,440,460,480)

interval<-1
numbers<-seq(from=0,to=499,by=interval)  #indicate start, end, and interval here
dir<-"/scratch/oliveros/trogons/"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=2000m,walltime=168:00:00
#PBS -M oliveros@ku.edu
#PBS -m n
#PBS -r n
#PBS -o /dev/null"

for(i in numbers)
{
	runfile<-paste("trogons.cf.",sprintf("%03d",i),sep="")

	pbsn<-paste("#PBS -N",runfile)
	pbsd<-paste("#PBS -d",dir)
	pbse<-paste("#PBS -e ",dir,runfile,".err",sep="")
	output<-paste(dir,"trogons_cloudforest_bootstrap",i,sep="")
	dir.create(output)
	command<-paste("python /scratch/oliveros/CloudForest/cloudforest/cloudforest_mpi2.py /scratch/oliveros/trogons/trogons_phylip/",output, "bootstraps /tools/cluster/6.2/cloudforest/0.1/bin/phyml --genetrees /scratch/oliveros/trogons/trogons_cloudforest/genetrees.tre --bootreps", interval, "--parallelism single")
	a<-paste(pbsn,pbs,pbsd,pbse,command,sep="\n")
	write.table(a,runfile, row.names=F,col.names=F,quote=F)
}

# To submit jobs to cluster:
for i in trogons.cf.???; do qsub $i; done

EXABAYES
Formatting Exabayes files from mrbayes files.

head -3 trogons.inc.min.75percent.mrbayes.nex | grep 'dimensions' | sed -r 's/\tdimensions ntax=([0-9]+) nchar=([0-9]+);/\1 \2/' > trogons.inc.min.75percent.phylip
[head -(numtaxa + 5) tail (numtaxa)]
head -18 trogons.inc.min.75percent.mrbayes.nex | tail -13 >> trogons.inc.min.75percent.phylip
cat trogons.inc.min.75percent.mrbayes.nex | grep 'charset' | sed -r 's/\tcharset (.+);/DNA, \1/' > trogons.inc.min.75percent.part


Create Exabayes single chain jobs:

numbers<-seq(from=1,to=1,by=1)  #indicate start, end, and interval here
dir<-"/scratch/oliveros/trogons/trogons_inc_min_75percent_exabayes/"
phylipfile<-paste(dir,"trogons.inc.min.75percent.phylip",sep="")
partitionsfile<-paste(dir,"trogons.inc.min.75percent.part",sep="")
configfile<-paste(dir,"trogons.inc.min.75percent.config.nex",sep="")

pbs<-"#PBS -l nodes=1:ppn=20:avx,mem=30000m,walltime=168:00:00
#PBS -M oliveros@ku.edu
#PBS -m abe
#PBS -r n
#PBS -o /dev/null"
pbsd<-paste("#PBS -d",dir)

for(i in numbers)
{
	runfile<-paste("trogons.exb",sprintf("%01d",i),sep="")
	pbsn<-paste("#PBS -N",runfile)
	pbse<-paste("#PBS -e ",dir,runfile,".err",sep="")
	output<-paste(dir,runfile,".out",sep="")
	name<-paste("trogons.c",i,sep="")
	randseed<-floor(runif(1)*799736+1111)
	command<-paste("unbuffer mpirun -np 20 exabayes_avx -f",phylipfile,"-q",partitionsfile,"-n",name,"-R 1 -C 1 -s",randseed,"-c",configfile,"> ",output)
	a<-paste(pbsn,pbs,pbsd,pbse,command,sep="\n")
	write.table(a,runfile, row.names=F,col.names=F,quote=F)
}

# To submit jobs to cluster:
for i in trogons.exb.???; do qsub $i; done


XX python ~/phyluce/bin/align/get_informative_sites1.py --input trogons_renamed/ --output-prefix trogons --input-format nexus --cores 12 --log-path trogons_log/


UTILITIES:
For cleaning trinity assemblies:

for iter in trinity-assemblies-her* ; do rm -r $iter/*/chrysalis ; rm $iter/*/both.fa* ; rm $iter/*/bowtie.* ; rm $iter/*/inchworm.* ; rm $iter/*/iworm_scaffolds* ; rm $iter/*/jellyfish.*; rm $iter/*/scaffolding*; rm $iter/*/target.*; done

For cleaning trinity coverage files:
for iter in trinity-assemblies; do rm $iter/*/bwa-index-file.log; rm $iter/*_*/*_*; rm $iter/*/contigs.dict; rm $iter/*/contigs.fasta.*; done


sudo apt-get update && sudo apt-get upgrade -y

sudo shutdown -h now  OR  sudo shutdown -r now
 
SCREEN
======
Screen can be started by typing

screen

in a terminal. Press Enter after reading the introductory text.

Virtual terminals in Screen can be manipulated by pressing the Ctrl+A key combination, and subsequently pressing a key to execute one of the commands given below:

    c creates a new virtual console

    n switches to the next available virtual console

    p switches back to the previous virtual console

    " lists all available virtual consoles and their assigned numbers
    hitting a number key brings the corresponding virtual console to the foreground

    Esc lets you scroll back and forth in your terminal output

    d detatches the current screen sessions and brings you back to the normal terminal 

When a Screen session is detached, the processes that were running inside it aren't stopped. You can re-attach a detached session by typing

screen -r

in a terminal 
 
