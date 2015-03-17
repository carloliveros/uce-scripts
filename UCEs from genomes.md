This document explains how you can extract UCE data from genomes from the Jarvis et al. 2014 Science paper.

### Download genomes

You can download bird genomes from http://gigadb.org/dataset/101000.  Download the assembled scaffolds for the relevant species.

### Convert to 2bit format

You will need to unzip the downloaded files then convert them from fasta into 2bit format.

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

Align probe sequences to all genomes.  Specify an sqlite database file, a lastz output file, the probe file, the list of genome scaffolds, and the path to the location of the genomes in 2bit format.

```
python /public/uce/phyluce/bin/align/run_multiple_lastzs_sqlite.py jarvis.sqlite jarvis_lastz uce-5k-probes.fasta --scaffoldlist AcaChl BucRhi ColStr LepDis ManVit MerNub NesNot PicPub --genome-base-path genomes/ --cores 12 --coverage 67 --identity 80
```

Set up a conf file (say, genomes.conf) that contains a mapping of the scaffold names and genome files. 

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

Slice out fastas from each respective genome

```
python /public/uce/phyluce/bin/share/slice_sequence_from_genomes2.py genomes.conf jarvis_lastz jarvis_1000_flank_fasta --flank=1000 --name-pattern "uce-5k-probes.fasta_v_{}.lastz.clean"
```