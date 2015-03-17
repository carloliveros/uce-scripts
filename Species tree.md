# Summary species tree analysis

### by Carl H. Oliveros

Perform the following steps on the KU cluster in a working directory such as /scratch/username/dataset1/dataset1_speciestree.

## 1. Prepare bootstrap replicates

This step assumes that you have gene trees estimated from 500 bootstrap replicates of the data saved in a file called 500-bootreps.tre.  If your dataset has L loci and N bootstrap replicates, the file should contain:

```
genetree_from_locus1_of_bootstrap1
genetree_from_locus2_of_bootstrap1
...
genetree_from_locusL_of_bootstrap1
genetree_from_locus1_of_bootstrap2
genetree_from_locus2_of_bootstrap2
...
genetree_from_locusL_of_bootstrap2
...
...
...
genetree_from_locus1_of_bootstrapN
genetree_from_locus2_of_bootstrapN
...
genetree_from_locusL_of_bootstrapN
```

Split and rename bootstrap replicates.

```
split --lines=1933 --suffix-length=3 -d 500-bootreps.tre boot
```

The --lines argument contains the number of loci in your dataset.  This command will split this big file (500-bootreps.tre) into boot000, boot001, ... bootN-1, each with a set of L genetrees.

## 2. Clean up phylip trees and root them

Create cleaning and rooting R scripts.  The following R script will clean up genetrees from Cloudforest (and save them as boot???.phy) and root them with the specified outgroup (and save them as boot???.rooted). 

```
numboot<-500  # number of bootstrap replicates
increment<-20   # increment
outgroup<-"GenSpe"  #outgroup taxon

start<-0
end<-increment - 1

while (end < numboot)
{
	fname<-paste("root",sprintf("%03d.R",start),sep="")
	phybase<-'library("phybase")'
	forloop<-paste("for(i in ",start,":",end,")",sep="")
	script<-paste('{
	\tprint(i)
	\tgenetreefname<-paste("boot",sprintf("%03d",i),sep="")
	\tgenetreefnamephy<-paste("boot",sprintf("%03d",i),".phy",sep="")
	\t# redundant reading, writing, reading, necessary to get rid of unnecessary characters
	\ttemp<-read.tree(file=genetreefname)
	\twrite.tree(temp,file=genetreefnamephy)

	\t# rooting trees
	\toutfile<-paste("boot",sprintf("%03d",i),".rooted",sep="")
	\ta<-read.tree.string(genetreefnamephy,format="phylip")$tree
	\tngenetree<-length(a)
	\td<-rep("",ngenetree)
	\tfor(k in 1:ngenetree)
	\t{
		\t\tspname<-species.name(a[k])
		\t\tnspecies<-length(spname)
		\t\toutgroup<-which(spname=="',outgroup,'")
		\t\tb<-read.tree.nodes(a[k],spname)$node
		\t\td[k]<-write.subtree(2*nspecies-1,root.tree(b,outgroup),spname,2*nspecies-1)
	\t}
	\twrite.tree.string(d,format="phylip",file=outfile)
}',sep="")
	a<-paste(phybase,forloop,script,sep="\n")
	write.table(a,fname, row.names=F,col.names=F,quote=F)
	start<-start + increment
	end<-end + increment
}
```

The following R script creates cleaning/rooting job scripts for the cluster.

```
numbers<-seq(from=0,to=480,by=20)  #indicate start, end, and interval here
dir<-"/scratch/username/dataset1/dataset1_speciestree"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=5000m,walltime=48:00:00
#PBS -M username@ku.edu
#PBS -r n
#PBS -m n
#PBS -j oe"

for(i in numbers)
{
	runfile<-paste("dataset1.root.",sprintf("%03d",i),sep="")

	pbsn<-paste("#PBS -N",runfile)
	pbsd<-paste("#PBS -d",dir)
	pbso<-paste("#PBS -o ",dir,"/",runfile,".out",sep="")
	rootfile<-paste("root",sprintf("%03d",i),".R",sep="")
	command<-paste("R --vanilla <",rootfile)
	a<-paste(pbsn,pbs,pbsd,pbso,command,sep="\n")
	write.table(a,runfile, row.names=F,col.names=F,quote=F)
}
```

To submit jobs to cluster:

```
for i in dataset1.root.???; do qsub $i; done
```

Wait for all cleaning/rooting jobs to finish before running any of the analyses below.  It's ok to get them set up while waiting for the jobs to finish.

## 3. STAR, STEAC, NJst

The following R script creates R scripts to infer the STAR, STEAC and NJst trees for each set of gene trees.  Specify the outgroup name in the script.

```
numboot<-500  # number of bootstrap replicates
increment<-1   # increment
outgroup<-"GenSpe"  # outgroup taxon
wd<-"/scratch/username/dataset1/dataset1_speciestree"  # working directory

start<-0
end<-increment - 1

while (end < numboot)
{
	fname<-paste("starsteacnjst",sprintf("%03d",start),".R",sep="")
	phybase<-'library("phybase")'
	workdir<-paste('setwd("',wd,'")',sep="")
	readtree<-paste('mytrees<-read.tree(paste("boot", sprintf("%03d",',start,'), ".rooted", sep=""))',sep="")
	importnjst<-'source("/scratch/oliveros/NJst.R")'
	details<-'taxaname<-mytrees[[1]]$tip.label
speciesname<-taxaname
ntaxa<-length(taxaname)
ngene<-length(mytrees)'
	outg1<-paste('outgrouptaxon<-"',outgroup,'"',sep="")
	outg2<-'print(paste("Outgroup taxon", outgrouptaxon))'
	forloop<-paste("for(i in ",start,":",end,")",sep="")
	script = '{
	\tgenetreefname<-paste("boot", sprintf("%03d",i), ".rooted", sep="")
	\ttreestring<-read.tree.string(genetreefname,format="phylip")
	\ttrees<-treestring$tree
	
	\tspecies.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
	\tdiag(species.structure)<-1
	\tprint(paste("Estimating STAR tree for", genetreefname))
	\tstar<-star.sptree(trees, speciesname, taxaname, species.structure, outgroup=outgrouptaxon,method="nj")
	\tstarfname<-paste("boot",sprintf("%03d",i),".star.tre", sep="")
	\twrite.table(star,starfname,row.names=F,col.names=F,quote=F,append=FALSE)

	\tspecies.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
	\tdiag(species.structure)<-1
	\tprint(paste("Estimating STEAC tree for", genetreefname))
	\tsteac<-steac.sptree(trees, speciesname, taxaname, species.structure, outgroup=outgrouptaxon,method="nj")
	\tsteacfname<-paste("boot",sprintf("%03d",i),".steac.tre", sep="")
	\twrite.table(steac,steacfname,row.names=F,col.names=F,quote=F,append=FALSE)
	
	\tgenetreefname<-paste("boot", sprintf("%03d",i), ".phy", sep="")
	\ttreestring<-read.tree.string(genetreefname,format="phylip")
	\ttrees<-treestring$tree
	
	\tspecies.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
	\tdiag(species.structure)<-1
	\tprint(paste("Estimating NJst tree for", genetreefname))
	\tnjsttree<-NJst(trees, speciesname, taxaname, species.structure)
	\tnjstfname<-paste("boot",sprintf("%03d",i),".njst.tre", sep="")
	\twrite.table(njsttree,njstfname,row.names=F,col.names=F,quote=F,append=FALSE)

}'
	a<-paste(phybase,workdir,readtree,importnjst,details,outg1,outg2,forloop,script,sep="\n")
	write.table(a,fname, row.names=F,col.names=F,quote=F)
	start<-start + increment
	end<-end + increment
}
```

The following R script creates job scripts to run the STAR/STEAC/NJst R scripts on the cluster.

```
numbers<-seq(from=0,to=499,by=1)  #indicate start, end, and interval here

dir<-"/scratch/username/dataset1/dataset1_speciestree"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=5000m,walltime=48:00:00
#PBS -M username@ku.edu
#PBS -r n
#PBS -m n
#PBS -j oe"

for(i in numbers)
{
	runfile<-paste("dataset1.ssn.",sprintf("%03d",i),sep="")
	pbsn<-paste("#PBS -N",runfile)
	pbsd<-paste("#PBS -d",dir)
	pbso<-paste("#PBS -o ",dir,"/",runfile,".out",sep="")
	ssnfile<-paste("starsteacnjst",sprintf("%03d",i),".R",sep="")
	command<-paste("R --vanilla <",ssnfile)
	a<-paste(pbsn,pbs,pbsd,pbso,command,sep="\n")
	write.table(a,runfile, row.names=F,col.names=F,quote=F)
}
```

To submit jobs to cluster:

```
for i in dataset1.stst.???; do qsub $i; done
```

## 4. ASTRAL

The following R script creates ASTRAL job scripts for the cluster.  Make sure you have the correct astral command invocation in the astralcom variable.

```
interval<-50  #indicate interval here
numbers<-seq(from=0,to=450,by=interval)  #indicate start, end here
astralcom<-"unbuffer java -jar /scratch/oliveros/Astral/astral.4.4.0.jar"
dir<-"/scratch/username/dataset1/dataset1_speciestree"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=25000m,walltime=96:00:00
#PBS -M username@ku.edu
#PBS -m n
#PBS -r n
#PBS -j oe
#PBS -o /dev/null"

for(i in numbers)
{
	runfile<-paste("dataset1.ast.",sprintf("%03d",i),sep="")
	pbsn<-paste("#PBS -N",runfile)
	pbsd<-paste("#PBS -d",dir)
	comlist<-""
	end<-i+interval-1
	for(j in i:end)
	{
		phylipfile<-paste("boot",sprintf("%03d",j),".phy",sep="")
		astralfile<-paste("boot",sprintf("%03d",j),".phy.astral.tre",sep="")
		outfile<-paste("boot",sprintf("%03d",j),".phy.astral.out",sep="")
		command<-paste(astralcom,"-i",phylipfile,"-o",astralfile,">",outfile)
		comlist<-paste(comlist,command,sep="\n")
	}
	a<-paste(pbsn,pbs,pbsd,comlist,sep="\n")
	write.table(a,runfile, row.names=F,col.names=F,quote=F)	
}
```

To submit jobs to cluster:

```
for i in dataset1.ast.???; do qsub $i; done
```

Collecting ASTRAL run times

This assumes Astral output files are named *.astral.out.
Before saving this to a file, it's good idea to pipe it to wc to check that you
have that correct number of run times.

```
cat `ls *.astral.out` | grep "Optimal tree" | sed -e 's/Optimal tree inferred in //' -e 's/ secs//' > astral.runtimes.txt
```

## 5 MP-EST

Create MPEST Control Files

Here's how you can create a species-allele table from Astral output from the command line:

```
cat boot000.phy.astral.out |grep "Taxa" |sed 's/Taxa: //'|sed -r 's/(\w+)/\1 1 \1\n/g' |sed 's/, //' | tr '[]' '\n'
```

Paste the appropriate species-allele table in the R script below.  Enter the appropriate number 
of taxa and number of genes.  The following R script 

```
nsim<-499  # number of replicates - 1
ntaxa<-5  # number of taxa
ngenes<-979  # number of genes

#species-allele table below
c<-"GenAbc 1 GenAbc
GenDef 1 GenDef
GenGhi 1 GenGhi
GenJkl 1 GenJkl
GenSpe 1 GenSpe"

for(i in 0:nsim)
{
	file<-paste("control",sprintf("%03d",i),sep="")   #filename of control file
	treefile<-paste("boot",sprintf("%03d",i),".rooted",sep="")  #filename of input rooted tree
	b<-floor(runif(1)*799736+1111)  # creates random seed number
	a<-paste(treefile,"0",b,paste(ngenes, ntaxa),c ,"0",sep="\n")  # contents of control file including num genetrees and num species
	write.table(a, file,row.names=F,col.names=F,quote=F)
}
```

The following R script creates MP-EST job script files for the cluster.

```
interval<-1  #indicate interval here
numbers<-seq(from=0,to=499,by=interval)  #indicate start, end here
#numbers<-c(9,28,89,109,129,189,248,288,308,348,369,468,489)
dir<-"/scratch/username/dataset1/dataset1_speciestree"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=5000m,walltime=96:00:00
#PBS -M username@ku.edu
#PBS -r n
#PBS -m n
#PBS -j oe"

for(i in numbers)
{
	runfile<-paste("dataset1.mpest.",sprintf("%03d",i),sep="")

	pbsn<-paste("#PBS -N",runfile)
	pbsd<-paste("#PBS -d",dir)
	pbso<-paste("#PBS -o ",dir,"/",runfile,".out",sep="")
	comlist<-""
	end<-i+interval-1
	for(j in i:end)
	{
		contfile<-paste("control",sprintf("%03d",j),sep="")
		command<-paste("mpest",contfile)
		comlist<-paste(comlist,command,sep="\n")
	}
	a<-paste(pbsn,pbs,pbsd,pbso,comlist,sep="\n")
	write.table(a,runfile, row.names=F,col.names=F,quote=F)
}
```

Submit MP-EST jobs to cluster:

```
for iter in dataset1.mpest.*; do qsub $iter; done
```

Collect all MP-EST trees

```
cat boot???.rooted.tre | grep "tree mpest" > summary
```

Create nexus file with all trees

```
cat zhead summary ztail > mpest.500.tre
```

Note:  Make zhead and ztail from nexus headers and tail of each output file.

Collect MP-EST run times

This assumes you have saved MPEST output in files named *rooted.tre
Before saving this to a file, good idea to pipe it to wc to check that you
have that correct number of run times.

```
cat $(echo $(ls *rooted.tre)) | grep "Analysis completed" | sed -e 's/\[Analysis completed in //' -e 's/ seconds\]//' > mpest.runtimes.txt
```

## 6. Running STAR, STEAC, ASTRAL, MP-EST on original data

The following R script performs Steps 2 and 3 above on the original data.  Run the script in the same directory as your Cloudforest output directory (e.g., ./dataset1_cloudforest/) where the genetrees.tre file is saved.

```
outgrouptaxon<-"GenSpe"  # indicate outgroup name here

library("phybase")
# clean up phylip file
temp<-read.tree(file="genetrees.tre")  
write.tree(temp,file="genetrees.phy")
	
# rooting trees 
a<-read.tree.string("genetrees.phy",format="phylip")$tree
ngenetree<-length(a)
d<-rep("",ngenetree)
for(k in 1:ngenetree) 
{
	spname<-species.name(a[k])
	nspecies<-length(spname)
	outgroup<-which(spname==outgrouptaxon)  
	b<-read.tree.nodes(a[k],spname)$node
	d[k]<-write.subtree(2*nspecies-1,root.tree(b,outgroup),spname,2*nspecies-1)
}
write.tree.string(d,format="phylip",file="genetrees.rooted")	

#STAR, STEAC

mytrees<-read.tree("genetrees.rooted")
taxaname<-mytrees[[1]]$tip.label
speciesname<-taxaname
ntaxa<-length(taxaname)
ngene<-length(mytrees)

print(paste("Outgroup taxon", outgrouptaxon))

treestringphy<-read.tree.string("genetrees.phy",format="phylip")
treesphy<-treestringphy$tree
species.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
diag(species.structure)<-1
print("Estimating STAR tree")
star<-star.sptree(treesphy, speciesname, taxaname, species.structure, outgroup=outgrouptaxon,method="nj")
write.table(star,"genetrees.star.tre",row.names=F,col.names=F,quote=F,append=TRUE)

treestringrooted<-read.tree.string("genetrees.rooted",format="phylip")
treesrooted<-treestringrooted$tree
species.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
diag(species.structure)<-1
print("Estimating STEAC tree")
steac<-steac.sptree(treesrooted, speciesname, taxaname, species.structure, outgroup=outgrouptaxon,method="nj")
write.table(steac,"genetrees.steac.tre",row.names=F,col.names=F,quote=F,append=TRUE)
```

ASTRAL and MPEST need to be run from the command line.

```
java -jar /scratch/oliveros/Astral/astral.4.4.0.jar -i genetrees.phy -o genetrees.astral.tre 2>&1 | tee genetrees.phy.astral.out

mpest control
```

## Miscellaneous stuff

Installing R package from source

```
install.packages("file_name_and_path", repos = NULL, type="source")
```

Where `file_name_and_path` would represent the full path and file name of the package. 

Monitoring progress of STAR, STEAC, ASTRAL, and MP-EST runs

```
ls boot???.star.tre | wc -l
ls boot???.steac.tre | wc -l
ls boot???.phy.astral.tre | wc -l
cat boot???.rooted.tre | grep "tree mpest" | wc -l
```

Summarizing trees with Dendropy

```
sumtrees.py -o dataset1.star.con.tre boot???.star.tre
sumtrees.py -o dataset1.steac.con.tre boot???.steac.tre
sumtrees.py -o dataset1.astral.con.tre boot???.phy.astral.tre
sumtrees.py -o dataset1.mpest.con.tre mpest.500.tre
```

Add the -f 0.70 option to sumtrees.py to collapse nodes below 70% support.

Getting frequencies of splits with Dendropy from the Python interpreter

```
import dendropy
trees=dendropy.TreeList()
split_leaves = ['ZosMad', 'ZosEry', 'ZosGri']
f = trees.frequency_of_split(labels=split_leaves)
print('Frequency of split %s: %s' % (split_leaves, f))

      Frequency of split ['ZosMad', 'ZosEry', 'ZosGri']: 0.342153846154
      
split_leaves1 = ['ZosMad', 'ZosEry', 'ZosGri', 'SpeMel']
f1 = trees.frequency_of_split(labels=split_leaves1)
print('Frequency of split %s: %s' % (split_leaves1, f1))

      Frequency of split ['ZosMad', 'ZosEry', 'ZosGri', 'SpeMel']: 0.161230769231
```
