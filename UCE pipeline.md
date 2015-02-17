# UCE Pipeline for Phylogenetic Studies

### by Carl H. Oliveros

Based on:

https://github.com/faircloth-lab/phyluce/blob/master/docs/uce-processing.rst

Read these pages in conjunction with this document.

## Outline

1 - Run illumiprocessor.py

2 - Assemble reads into contigs

3 - Extract UCE contigs 

4 - Inspect data using sqlite

5 - Assemble dataset

6 - Calculate coverage for Trinity-aligned datasets

7 - Sequence alignment

8 - Formatting data for phylogenetic analysis

## STEP 1 - Run illumiprocessor.py

The illumiprocessor python script trims adapter contamination and performs 
quality checking of reads. Create a folder in your work folder and copy all raw
read files in that folder.

Create a configuration file for illumiprocessor.py.  USE ALL SMALL LETTERS FOR NAMES.

For illumiprocessor Version 2:

```
[adapters]
i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT

[tag sequences]
i7-0708:TGCAAGAC
i5-03A:AACACCAC
i5-03B:TGAGCTGT
i5-03D:TGACAACC
i5-03H:CAGTGCTT
i5-03E:TGTTCCGT

[tag map]
genus1_species1_15512_GTCTTGCA-AACACCAC:i7-0708,i5-03A
genus2_species2_15700_GTCTTGCA-TGAGCTGT:i7-0708,i5-03B
genus3_species3_26718_GTCTTGCA-TGACAACC:i7-0708,i5-03D
genus4_species4_15475_GTCTTGCA-CAGTGCTT:i7-0708,i5-03H
genus5_species5_8703_GTCTTGCA-TGTTCCGT:i7-0708,i5-03E

[names]
genus1_species1_15512_GTCTTGCA-AACACCAC:genus1_species1_15512
genus2_species2_15700_GTCTTGCA-TGAGCTGT:genus2_species2_15700
genus3_species3_GTCTTGCA-TGACAACC:genus3_species3_26718
genus4_species4_15475_GTCTTGCA-CAGTGCTT:genus4_species4_15475
genus5_species5_8703_GTCTTGCA-TGTTCCGT:genus5_species5_8703

```

For illumiprocessor version 1:

```
[adapters]
N7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
N5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT

[indexes]
N7_07_01:GAGAAGGT
N7_07_02:ACTGGTGT
N7_07_03:CTTCCTTC
N7_07_04:TCACTCGA
N5_03_A:AACACCAC
N5_03_B:TGAGCTGT
N5_03_C:CACAGGAA
N5_03_D:TGACAACC
N5_03_E:TGTTCCGT
N5_03_F:CCTAGAGA
N5_03_G:GCATAACG
N5_03_H:CAGTGCTT

[params]
separate reads:True
read1:{name}_L001_R1_001.fastq.gz
read2:{name}_L001_R2_001.fastq.gz

[combos]
genus1_species1_1234:N7_07_01,N5_03_A

[remap]
genus1_species1_1234_ACCTTCTC-AACACCAC:genus1_species1_1234

```

Save the file on the working directory as preprocess.conf.

Navigate to your working directory and create a directory that will contain the output (only for version 1).

```
mkdir cleaned_reads
```

Call the program and insert arguments for the input (the folder containing raw reads),
the output directory, your config file, and the number of processors to be used. 

version 2
```
illumiprocessor --input raw_reads --output cleaned_reads --config preprocess.conf --trimmomatic /anaconda/jar/trimmomatic.jar --log-path logs --cores 12
```

version 1
```
illumiprocessor.py raw_reads cleaned_reads preprocess.conf --remap --complex --cores 12
```


## STEP 2 - Assemble reads into contigs

In this step, you could use Trinity (preferred) or Velvet.

### 2a - Run Trinity to Assemble Contigs for each species

Run the python script assemblo_trinity.py to assemble contigs for
each species.  You need to create a run config file, that tells the script
where the read data are and what you want the final data named. Below is an 
example config file that was named "trinity-assemblies.conf":

```
[samples]
genus1_species1_1234:/path-to/cleaned-reads/genus1_species1_1234
genus2_species2_5678:/path-to/cleaned-reads/genus2_species2_5678
genus3_species3_9012:/path-to/cleaned-reads/genus3_species3_9012
```

The name you want things called on output is on the left side. Within each of 
the paths on the right side, there is a folder called "split-adapter-quality-
trimmed" (for PE data) or "interleaved-adapter-quality-trimmed" where the reads
should be located.  

Invoke assemblo_trinity.py giving the name of the config file, the name of
an output path, the name of the subfolder where reads are located, and the 
number of cores you want the program to run.

```
mkdir trinity_assemblies

python /public/uce/phyluce/bin/assembly/assemblo_trinity.py --conf trinity.conf --output trinity_assemblies --subfolder 'split-adapter-quality-trimmed' --cores 2

```

This should assemble everything in it's own folder in "trinity-assemblies".  
There will also be a "contigs" folder in that directory that contains 
symlinks to the contigs that are assembled, named appropriately for the 
`"match_contigs_to_probes.py"` script.  The --clean option is optional and is
used to clean up unnecessary trinity files that take up a lot of disk space.


### 2b - Run Velvet to Assemble Contigs for each species

For detailed instructions on this step, visit:
https://github.com/faircloth-lab/phyluce/blob/master/docs/pre-processing-assembly.rst

Move into the cleaned-reads folder.

```	
cd /public/uce/work/cleaned-reads/
```

To run Velvet Optimiser (a wrapper script that determines the best settings for 
Velvet given your dataset), use the following commands. First, print the names 
of the cleaned files you will be using.

```
VelvetOptimiser.pl --s 69 --e 75 --k=n50 --c=tbp -t 6 -a -f "-fastq.gz -shortPaired $(echo $( ls `*/split-adapter-quality-trimmed/*`-READ?.fastq.gz ) ) -short $(echo $( ls `*/split-adapter-quality-trimmed/*`-READ-singleton.fastq.gz ) )"
```

Once you have an idea of the kmer range that will work with your data, run 
assemblo.py, which will use velvet to assemble contigs for each species.

```
cd /public/uce/work/
assemblo.py ./cleaned-reads 71 83 6
```

In the command above specify the kmer range and the number of processors you 
want to use.  The resulting contigs (symlinks) will be found in the contigs 
directory in the cleaned reads directory.


## STEP 3 - Extract UCE contigs 

Download the probe set

```
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta
```

Match assembled contigs to self-matched probes and extract the contigs matching 
UCE probes.  

```
NEW VERSION
python /public/uce/phyluce/bin/assembly/match_contigs_to_probes.py --contigs /path/to/assembly/contigs/ --probes uce-5k-probes.fasta --output lastz --log-path logs

OLD VERSIONS
match_contigs_to_probes.py /public/uce/work/cleaned-reads/contigs /public/uce/work/LSU-Custom-Array-Jan-2013.fasta /public/uce/work/lastz --regex "_p[1-9]+$" --repl "" --dupefile /public/uce/work/LSU-Custom-Array-Jan-2013-to-self.fasta

python /public/uce/phyluce/bin/assembly/match_contigs_to_probes.py --contigs trinity_assemblies/contigs/ --probes 5k-probes.fasta  --dupefile 5k-probes-to-self.fasta --output ./lastz/ --regex "(chr\w+)(?:_probe\d+)"
```

To convert log into tab delimited text file:

```
cat path-to/match_contigs_to_probes.log |grep 'dupe probe matches' | sed -r 's/.+ - INFO - //' | sed -r 's/: /\t/' | sed -r 's/ \(/\t/' |sed -r 's/\) uniques of /\t/' |sed -r 's/ contigs, /\t/' |sed -r 's/ dupe probe matches, /\t/' | sed -r 's/ UCE loci removed for matching multiple contigs, /\t/' |sed -r 's/ contigs removed for matching multiple UCE loci//' > contig.summary.txt
```


## STEP 4 - Inspect data using sqlite

See Brant's github page for more details.  It's all there!

https://github.com/faircloth-lab/phyluce/blob/master/docs/uce-processing.rst


## STEP 5 - Assemble dataset

Create a dataset configuration in a file, e.g. datasets.conf, and save it
in the working directory.  It should look something like:

```
[dataset1]
genus1_species1
genus2_species2
genus3_species3
genus4_species4
genus5_species5

[dataset2]
genus4_species4
genus5_species5
genus6_species6
```

Use short (max 8 chars) and unique dataset names.  These dataset names will be used for directory names and cluster job names later.

### A. Run a query on the database:

```
python /public/uce/phyluce/bin/assembly/get_match_counts.py --locus-db lastz/probe.matches.sqlite --taxon-list-config datasets.conf --taxon-group 'dataset1' --output dataset1.conf --log-path dataset1_log
```

This will produce a config file listing taxon names from dataset1 as well as a 
list of loci for which all taxa in the list has data.  If an incomplete data
matrix is desired, add the flag --incomplete-matrix as below:

```
python /public/uce/phyluce/bin/assembly/get_match_counts.py --locus-db lastz/probe.matches.sqlite --taxon-list-config datasets.conf --taxon-group 'dataset1' --output dataset1.inc.conf --log-path dataset1_inc_log --incomplete-matrix
```

### B. Generate a fasta file for your dataset.

```
python /public/uce/phyluce/bin/assembly/get_fastas_from_match_counts.py --contigs trinity_assemblies/contigs/ --locus-db lastz/probe.matches.sqlite --match-count-output dataset1.conf --output dataset1.fasta --log-path dataset1_log
```

If working with an incomplete data matrix, add the flag --incomplete-matrix and
the path to an incomplete.nostrict file that will hold missing locus information 
as below:

```
python /public/uce/phyluce/bin/assembly/get_fastas_from_match_counts.py --contigs trinity_assemblies/contigs/ --locus-db lastz/probe.matches.sqlite --match-count-output dataset1.inc.conf --output dataset1.inc.fasta --incomplete-matrix dataset1.inc.nostrict --log-path dataset1_inc_log
```

## STEP 6 - Calculate coverage for Trinity-aligned datasets

For this step you will need:
- cleaned reads folder that contains a folder for each sample and raw reads 1
and 2 of each sample in their respective folder
- trinity assemblies folder that contains a folder for each sample and the 
contigs.fasta file
- a configuration file that maps each sample to the location of its raw reads
files location
- a database of contigs-uce mapping (probe.matches.sqlite)
- the output of `get_match_counts.py` run for an incomplete matrix (so that
coverage will be calculated for all enriched loci).

The following commands will help you copy the needed files from your cleaned-
reads folder and trinity-assemblies folder to the cluster.

```
mkdir dataset1_clean

dest=/media/Dicrurus/uce/dataset1_clean

[cd to cleaned reads folder]

for i in `ls`; do mkdir $dest/$i ; mkdir $dest/$i/split-adapter-quality-trimmed; cp $i/split-adapter-quality-trimmed/*READ?.fastq.gz $dest/$i/split-adapter-quality-trimmed; done;
```

```
mkdir dataset1_trinity

dest=/media/Dicrurus/uce/dataset1_trinity

[cd to trinity assemblies folder]

for i in `ls`; do mkdir $dest/$i ; cp -L $i/contigs.fasta $dest/$i; done;
```

Use the following commands in your PBS script:

To calculate coverage for all contigs:

```
python /scratch/oliveros/phyluce/bin/assembly/get_trinity_coverage.py --assemblies /scratch/username/dataset1/dataset1_trinity/ --assemblo-config /scratch/username/dataset1/trinity.conf --cores 4 --subfolder 'split-adapter-quality-trimmed' --log-path /scratch/username/dataset1/dataset1_covlog/
```

To calculate coverage for UCE contigs.  This module is ok to run on a single processor interactive session:

```
python /scratch/oliveros/phyluce/bin/assembly/get_trinity_coverage_for_uce_loci.py --assemblies /scratch/username/dataset1/dataset1_trinity/ --match-count-output /scratch/username/dataset1/dataset1.inc.conf --locus-db /scratch/username/dataset1/probe.matches.sqlite --output /scratch/username/dataset1/dataset1_coverage --log-path /scratch/username/dataset1/dataset1_covlog/
```

`*.jar` paths hardcoded in /phyluce/bwa.py

Add /scratch/oliveros/phyluce/:/scratch/oliveros/biopython-1.64 to PYTHONPATH

Summarizing output:

```
cat get_trinity_coverage.log | grep 'Processing\|mean coverage' | sed -r 's/.+Processing (\w+).+/\1/' | sed -r 's/.+INFO -\s+([0-9]+) contigs, mean coverage = ([0-9]+\.[0-9]), mean length = ([0-9]+\.[0-9])/\1\t\2\t\3/' | sed '$!N;s/\n/\t/' 

cat get_trinity_coverage_for_uce_loci.log | grep 'Processing\|mean trimmed coverage' | sed -r 's/.+Processing (\w+).+/\1/' | sed -r 's/.+INFO -\s+([0-9]+) contigs, mean trimmed length = ([0-9]+\.[0-9]), mean trimmed coverage = ([0-9]+\.[0-9])x, on-target bases \(uce contigs\) = ([0-9]+\.[0-9])%, unique reads aligned \(all contigs\) = ([0-9]+\.[0-9])%/\1\t\2\t\3\t\4\t\5/' | sed '$!N;s/\n/\t/'
```

Collecting run times (refine this)

```
cat get_trinity_coverage.log | grep 'Processing\|Completed' |sed -r 's/(.+) (.+),.+/\1\t\2/'
```

## STEP 7 - Sequence alignment

### A. Align the dataset.  

The following script will create aligned nexus files by locus.

```
python /public/uce/phyluce/bin/align/seqcap_align_2.py --fasta dataset1.fasta --output dataset1_aligned --output-format fasta --taxa 71 --aligner mafft --cores 12 --log-path dataset1_log
```

If working with an incomplete matrix, perform the three steps:

1.  Perform alignment as above but add the flag --incomplete-matrix:

```
python /public/uce/phyluce/bin/align/seqcap_align_2.py --fasta dataset1.inc.fasta --output dataset1_inc_aligned --taxa 71 --aligner mafft --cores 12 --log-path dataset1_inc_log --incomplete-matrix
```

2.  You can set a minimum number of taxa for a locus to be included in your 
dataset:

```
python /public/uce/phyluce/bin/align/get_only_loci_with_min_taxa.py --alignments dataset1_inc_aligned --taxa 71 --output dataset1_inc_min_75percent --percent 0.75 --cores 12 --log-path dataset1_inc_log
```

3.  Insert missing data designators for taxa missing from the alignment of a 
given locus:

```
python /public/uce/phyluce/bin/align/add_missing_data_designators.py --alignments dataset1_inc_min_75percent --output dataset1_inc_min_75percent_with_missing --output-format fasta --match-count-output dataset1.inc.conf --incomplete-matrix dataset1.inc.nostrict --cores 12 --log-path dataset1_inc_log
```

### B. Trim alignments using GBlocks

```
python /public/uce/phyluce/bin/align/get_gblocks_trimmed_alignments_from_untrimmed.py --alignments dataset1_aligned --output dataset1_gbtrimmed --b2 0.65 --cores 12 --log-path dataset1_log 
```

If working with an incomplete matrix:

```
python /public/uce/phyluce/bin/align/get_gblocks_trimmed_alignments_from_untrimmed.py --alignments dataset1_inc_min_75percent_with_missing --output dataset1_inc_min_75percent_gbtrimmed --b2 0.65 --cores 12 --log-path dataset1_inc_log 
```

### C. Remove loci names from the taxon names

After inspecting individual alignments, remove loci names from the taxon names in
the nexus files.

```
python /public/uce/phyluce/bin/align/remove_locus_name_from_nexus_lines.py --taxa 71 --alignment dataset1_gbtrimmed --output dataset1_renamed --cores 12 --log-path dataset1_log
```

OR (IF YOU NEED TO CONVERT TO SHORTER NAMES):

```
python /public/uce/phyluce/bin/align/rename_taxa_from_nexus_lines.py --taxa 71 --alignment dataset1_gbtrimmed --output dataset1_renamed --cores 12 --log-path dataset1_log
```

If working with an incomplete matrix:

```
python /public/uce/phyluce/bin/align/remove_locus_name_from_nexus_lines.py --taxa 71 --alignment dataset1_inc_min_75percent_gbtrimmed --output dataset1_inc_min_75percent_renamed --cores 12 --log-path dataset1_inc_log
```

### D. Summary stats

You can get summary stats from your aligned dataset.

```
python /public/uce/phyluce/bin/align/get_align_summary_data.py --alignments dataset1_renamed --input-format nexus --cores 12 --log-path dataset1_log
```

If working with an incomplete matrix:

```
python /public/uce/phyluce/bin/align/get_align_summary_data.py --alignments dataset1_inc_min_75percent_renamed --input-format nexus --cores 12 --log-path dataset1_inc_log
```

## STEP 8 - Formatting data for phylogenetic analysis

### A. RAxML

To assemble a non-partitioned dataset for RAxML:

```
python /public/uce/phyluce/bin/align/format_nexus_files_for_raxml.py --alignments dataset1_renamed/ --output dataset1_raxml --log-path dataset1_log
```

If working with an incomplete matrix:

```
python /public/uce/phyluce/bin/align/format_nexus_files_for_raxml.py --alignments dataset1_inc_min_75percent_renamed/ --output dataset1_inc_min_75percent_raxml --log-path dataset1_inc_log
```

### B. Cloudforest

If you wish to perform concatenated analysis (e.g. MrBayes) on an unpartitioned dataset, skip section B and proceed to the MrBayes section.  However, if you wish to perform concatenated analysis on a partitioned dataset (partitioned by model), then perform step B1 before going to the MrBayes section.  I have found that running Cloudforest on the cluster is way faster.

#### B1. Estimating gene trees

You need to perform the following steps if:
- you wish to perform concatenated analysis on a partitioned dataset
- you wish to estimate genetrees from bootstrapped data using Cloudforest 

First, convert the dataset into strict phylip format.

```
python /public/uce/phyluce/bin/align/convert_one_align_to_another.py --alignment dataset1_renamed/ --output dataset1_phylip/ --input-format nexus --output-format phylip --cores 12 --shorten-names --log-path dataset1_log
```

Note: If you do not specify --shorten-names, program will take first 10 characters of name.  Consider using `rename_taxa_from_nexus_lines.py`

If working with an incomplete matrix:

```
python /public/uce/phyluce/bin/align/convert_one_align_to_another.py --alignment dataset1_inc_min_75percent_renamed/ --output dataset1_inc_min_75percent_phylip/ --input-format nexus --output-format phylip --cores 12 --shorten-names --log-path dataset1_inc_log
```

Next, estimate gene trees and best fitting substitution models using CloudForest.

```
mkdir dataset1_cloudforest
python ~/CloudForest/cloudforest/cloudforest_mpi.py dataset1_phylip/ dataset1_cloudforest/ genetrees /usr/bin/phyml --parallelism multiprocessing --cores 12 2>&1 | tee dataset1_cloudforest/genetrees.out

```

If working with an incomplete matrix:

```
mkdir dataset1_inc_min_75percent_cloudforest
python ~/CloudForest/cloudforest/cloudforest_mpi.py dataset1_inc_min_75percent_phylip/ dataset1_inc_min_75percent_cloudforest/ genetrees /usr/bin/phyml --parallelism multiprocessing --cores 12 2>&1 | tee dataset1_inc_min_75percent_cloudforest/genetrees.out
```

For doing the above on the cluster, copy your phylip folder on to a working directory such as /scratch/username/dataset1, create the output directory (`/scratch/username/dataset1/dataset1_cloudforest or /scratch/username/dataset1_inc_min_75percent_cloudforest`), then use one of the following PBS scripts:

For a complete matrix:

```
#PBS -N dataset1.cf.gt
#PBS -l nodes=1:ppn=20:avx,mem=16000m,walltime=168:00:00
#PBS -M username@ku.edu
#PBS -m abe
#PBS -d /scratch/username/dataset1
#PBS -j oe
#PBS -o /dev/null
unbuffer python /scratch/oliveros/CloudForest/cloudforest/cloudforest_mpi2.py /scratch/username/dataset1/dataset1_phylip/ /scratch/username/dataset1/dataset1_cloudforest/ genetrees /tools/cluster/6.2/cloudforest/0.1/bin/phyml --parallelism multiprocessing --cores 20 > /scratch/username/dataset1/dataset1.cf.gt.out
```

For an incomplete matrix:

```
#PBS -N dataset1.cf.gt
#PBS -l nodes=1:ppn=20:avx,mem=16000m,walltime=168:00:00
#PBS -M username@ku.edu
#PBS -m abe
#PBS -d /scratch/username/dataset1
#PBS -j oe
#PBS -o /dev/null
unbuffer python /scratch/oliveros/CloudForest/cloudforest/cloudforest_mpi2.py /scratch/username/dataset1/dataset1_inc_min_75percent_phylip/ /scratch/username/dataset1/dataset1_inc_min_75percent_cloudforest/ genetrees /tools/cluster/6.2/cloudforest/0.1/bin/phyml --parallelism multiprocessing --cores 20 > /scratch/username/dataset1/dataset1.cf.gt.out
```

#### B2. Cloudforest bootstrapping

You need to perform this step if you wish to use Cloudforest to estimate genetrees from bootstrapped data.

For bootstrapping with Cloudforest:

```
python ~/CloudForest-master/cloudforest/cloudforest_mpi2.py dataset1_phylip/ dataset1_cloudforest_bootstrap/ bootstraps /usr/bin/phyml --genetrees dataset1_cloudforest/genetrees.tre --bootreps 500 --parallelism multiprocessing --cores 12 > /dev/null
```

You can run Cloudforest bootstrapping on the cluster.  I do this by generating bootstrapping job scripts using R in your working directory (/scratch/username/dataset1) and launching all these jobs.  Note that the cloudforest folder containing the genetrees.tre file should be in the working directory.

Create bootstrapping jobs in R:

```
interval<-1
numbers<-seq(from=0,to=499,by=interval)  #indicate start, end, and interval here
dir<-"/scratch/username/dataset1/"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=2000m,walltime=168:00:00
#PBS -M username@ku.edu
#PBS -m n
#PBS -r n
#PBS -o /dev/null"

for(i in numbers)
{
	runfile<-paste("dataset1.cf.",sprintf("%03d",i),sep="")

	pbsn<-paste("#PBS -N",runfile)
	pbsd<-paste("#PBS -d",dir)
	pbse<-paste("#PBS -e ",dir,runfile,".err",sep="")
	output<-paste(dir,"dataset1_cloudforest_bootstrap",i,sep="")
	dir.create(output)
	command<-paste("python /scratch/oliveros/CloudForest/cloudforest/cloudforest_mpi2.py /scratch/username/dataset1/dataset1_phylip/",output, "bootstraps /tools/cluster/6.2/cloudforest/0.1/bin/phyml --genetrees /scratch/username/dataset1/dataset1_cloudforest/genetrees.tre --bootreps", interval, "--parallelism single")
	a<-paste(pbsn,pbs,pbsd,pbse,command,sep="\n")
	write.table(a,runfile, row.names=F,col.names=F,quote=F)
}
```

To submit jobs to cluster:

```
for i in dataset1.cf.???; do qsub $i; done
```

Create a directory for performing species tree analyses and save all bootstrapping results into one file.

```
mkdir /scratch/username/dataset1/dataset1_speciestree

cat /scratch/username/dataset1/dataset1_cloudforest_bootstrap*/bootrep001 > /scratch/username/dataset1/dataset1_speciestree/500-bootreps.tre
```


### C. MrBayes

Perform the two steps below to create a MrBayes file with the data partitioned according to the type of best fitting substitution model.

First, strip the models from CloudForest output.

```
python /public/uce/phyluce/bin/genetrees/split_models_from_genetrees.py --genetrees dataset1_cloudforest/genetrees.tre --output dataset1.models.txt
```

If working with an incomplete matrix:

```
python /public/uce/phyluce/bin/genetrees/split_models_from_genetrees.py --genetrees dataset1_inc_min_75percent_cloudforest/genetrees.tre --output dataset1.inc.min.75percent.models.txt
```

Second, create a nexus file for MrBayes.

```
python /public/uce/phyluce/bin/align/format_nexus_files_for_mrbayes.py --alignments dataset1_renamed/ --models dataset1.models.txt --output dataset1.mrbayes.nex --unlink
```

If working with an incomplete matrix:

```
python /public/uce/phyluce/bin/align/format_nexus_files_for_mrbayes.py --alignments dataset1_inc_min_75percent_renamed/ --models dataset1.inc.min.75percent.models.txt --output dataset1.inc.min.75percent.mrbayes.nex --unlink
```

### D. ExaBayes

Formatting Exabayes files from mrbayes files.

```
head -3 dataset1.inc.min.75percent.mrbayes.nex | grep 'dimensions' | sed -r 's/\tdimensions ntax=([0-9]+) nchar=([0-9]+);/\1 \2/' > dataset1.inc.min.75percent.phylip
[head -(numtaxa + 5) tail (numtaxa)]
head -18 dataset1.inc.min.75percent.mrbayes.nex | tail -13 >> dataset1.inc.min.75percent.phylip
cat dataset1.inc.min.75percent.mrbayes.nex | grep 'charset' | sed -r 's/\tcharset (.+);/DNA, \1/' > dataset1.inc.min.75percent.part
```

Create Exabayes single chain jobs:

```
numbers<-seq(from=1,to=1,by=1)  #indicate start, end, and interval here
dir<-"/scratch/username/dataset1/dataset1_inc_min_75percent_exabayes/"
phylipfile<-paste(dir,"dataset1.inc.min.75percent.phylip",sep="")
partitionsfile<-paste(dir,"dataset1.inc.min.75percent.part",sep="")
configfile<-paste(dir,"dataset1.inc.min.75percent.config.nex",sep="")

pbs<-"#PBS -l nodes=1:ppn=20:avx,mem=30000m,walltime=168:00:00
#PBS -M username@ku.edu
#PBS -m abe
#PBS -r n
#PBS -o /dev/null"
pbsd<-paste("#PBS -d",dir)

for(i in numbers)
{
	runfile<-paste("dataset1.exb",sprintf("%01d",i),sep="")
	pbsn<-paste("#PBS -N",runfile)
	pbse<-paste("#PBS -e ",dir,runfile,".err",sep="")
	output<-paste(dir,runfile,".out",sep="")
	name<-paste("dataset1.c",i,sep="")
	randseed<-floor(runif(1)*799736+1111)
	command<-paste("unbuffer mpirun -np 20 exabayes_avx -f",phylipfile,"-q",partitionsfile,"-n",name,"-R 1 -C 1 -s",randseed,"-c",configfile,"> ",output)
	a<-paste(pbsn,pbs,pbsd,pbse,command,sep="\n")
	write.table(a,runfile, row.names=F,col.names=F,quote=F)
}
```

To submit jobs to cluster:

```
for i in dataset1.exb.???; do qsub $i; done
```

## UTILITIES

Getting informative sites

```
python /public/uce/phyluce/bin/align/get_informative_sites1.py --input dataset1_renamed/ --output-prefix dataset1 --input-format nexus --cores 12 --log-path dataset1_log/
```

For cleaning trinity assemblies:

```
for iter in trinity-assemblies-her* ; do rm -r $iter/*/chrysalis ; rm $iter/*/both.fa* ; rm $iter/*/bowtie.* ; rm $iter/*/inchworm.* ; rm $iter/*/iworm_scaffolds* ; rm $iter/*/jellyfish.*; rm $iter/*/scaffolding*; rm $iter/*/target.*; done

For cleaning trinity coverage files:
for iter in trinity-assemblies; do rm $iter/*/bwa-index-file.log; rm $iter/*_*/*_*; rm $iter/*/contigs.dict; rm $iter/*/contigs.fasta.*; done
```

Updating and shutting down machine remotely

```
sudo apt-get update && sudo apt-get upgrade -y

sudo shutdown -h now  OR  sudo shutdown -r now
```
