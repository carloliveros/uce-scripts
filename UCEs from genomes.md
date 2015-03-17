Convert downloaded scaffolds in fasta format to 2bit format

faToTwoBit Acanthisitta_chloris.fa AcaChl.2bit

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

while read p; do faToTwoBit $p; done < file.txt

## In silico alignment

- align probe sequences to all genomes
- run the hymenoptera probeset against genome assemblies (NOTE change in default coverage and %
  identity. --genome-base-path is the path to the location of the genomes, stored
  in 2bit (faToTwoBit from Kent source) format:

    python /public/uce/phyluce/bin/align/run_multiple_lastzs_sqlite.py jarvis.sqlite jarvis_lastz uce-5k-probes.fasta --scaffoldlist AcaChl BucRhi ColStr LepDis ManVit MerNub NesNot PicPub --cores 12 --genome-base-path genomes/ --coverage 67 --identity 80

- slice out fastas from each respective genome, after setting up the conf file `genome-sequence-location.conf`::

    python /public/uce/phyluce/bin/share/slice_sequence_from_genomes2.py genomes.conf jarvis_lastz jarvis_1000_flank_fasta --flank=1000 --name-pattern "uce-5k-probes.fasta_v_{}.lastz.clean"

    # conf file above contains:

[scaffolds]
AcaChl:genomes/AcaChl.2bit
BucRhi:genomes/BucRhi.2bit
ColStr:genomes/ColStr.2bit
LepDis:genomes/LepDis.2bit
ManVit:genomes/ManVit.2bit
MerNub:genomes/MerNub.2bit
NesNot:genomes/NesNot.2bit
PicPub:genomes/PicPub.2bit
