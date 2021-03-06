This document explains how you can extract UCE data from genomes from the Jarvis et al. 2014 Science paper.  The python scripts are in the working branch of https://github.com/faircloth-lab/phyluce.  Thanks to Brant Faircloth (http://faircloth-lab.org/) for the advice!

### Download genomes

You can download bird genomes from http://gigadb.org/dataset/101000.  Download the assembled scaffolds in fasta format for the species of interest.

### Convert to 2bit format

You will need to unzip the downloaded files then convert them from fasta into 2bit format.  faToTwoBit download is available at http://hgdownload.cse.ucsc.edu/admin/exe/.

```
gunzip *.gz
faToTwoBit Acanthisitta_chloris.fa AcaChl.2bit
```

You can batch convert files by looping through a text file (say file.txt) containing fasta names and 2bit names: 

```
Acanthisitta_chloris.fa AcaChl.2bit
Buceros_rhinoceros.fa BucRhi.2bit
Colius_striatus.fa ColStr.2bit
Leptosomus_discolor.fa LepDis.2bit
Manacus_vitellinus.fa ManVit.2bit
Merops_nubicus.fa MerNub.2bit
Nestor_notabilis.fa NesNot.2bit
Picoides_pubescens.fa PicPub.2bit
```

To loop through this text file with faToTwoBit:

```
while read p; do faToTwoBit $p; done < file.txt
```

### In silico alignment

Align probe sequences to all genomes.  You will need to install bx-python and add the bx-python library to your PYTHONPATH.  Specify an sqlite database file, a lastz output directory, the probe file, the list of genome scaffolds, and the path to the location of the genomes in 2bit format.  You will need to create the lastz output directory before running the python script.  You can also change the coverage and identity values to make the matching more or less stringent.

```
phyluce_probe_run_multiple_lastzs_sqlite \
    --db jarvis.sqlite \
    --output jarvis_lastz \
    --probefile uce-5k-probes.fasta \
    --scaffoldlist AcaChl BucRhi ColStr LepDis ManVit MerNub NesNot PicPub \
    --genome-base-path avian-genomes/ \
    --no-dir \
    --coverage 83 \
    --identity 92.5 \
    --cores 12
```

Set up a conf file (say, genomes.conf) that contains a mapping of the scaffold names and paths to the genome files. 

```
[scaffolds]
AcaChl:genomes/AcaChl.2bit
BucRhi:genomes/BucRhi.2bit
ColStr:genomes/ColStr.2bit
LepDis:genomes/LepDis.2bit
ManVit:genomes/ManVit.2bit
MerNub:genomes/MerNub.2bit
NesNot:genomes/NesNot.2bit
PicPub:genomes/PicPub.2bit
```

Slice out fastas from each respective genome with the script below.  Adjust the size of flanking regions to include using the flank parameter.

```
python ~/phyluce/bin/probes/phyluce_probe_slice_sequence_from_genomes \
    --conf genomes.conf \
    --lastz jarvis_lastz \
    --output jarvis_1000_flank_fasta \
    --name-pattern "uce-5k-probes.fasta_v_{}.lastz.clean" \
    --flank 1000
```

The fasta files produced can now be used for `match_contigs_to_probes.py`.
