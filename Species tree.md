# Species tree analysis pipeline

INSTALLING R PACKAGE FROM SOURCE

```
install.packages("file_name_and_path", repos = NULL, type="source")
```

Where `file_name_and_path` would represent the full path and file name of the package. 


## 1. Prepare bootstrap replicates

Split output from cloudforest

split --lines=1933 --suffix-length=3 -d 500-bootreps.tre boot

## 2. Clean up phylip trees and root them

Create cleaning and rooting R scripts

```
numboot<-500  # number of bootstrap replicates
increment<-20   # increment
outgroup<-"StaNig"  #outgroup taxon

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

Create rooting run files

```
numbers<-seq(from=0,to=480,by=20)  #indicate start, end, and interval here
dir<-"/scratch/oliveros/zosterops/zosterops_speciestree"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=5000m,walltime=48:00:00
#PBS -M oliveros@ku.edu
#PBS -r n
#PBS -m n
#PBS -j oe"

for(i in numbers)
{
	runfile<-paste("zosterops.root.",sprintf("%03d",i),sep="")

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
for i in zosterops.root.???; do qsub $i; done
```

## 3. STAR and STEAC

Create STAR/STEAC R scripts

```
numboot<-500  # number of bootstrap replicates
increment<-1   # increment
outgroup<-"StaNig"  # outgroup taxon
wd<-"/scratch/oliveros/zosterops/zosterops_speciestree"  # working directory

start<-0
end<-increment - 1

while (end < numboot)
{
	fname<-paste("starsteac",sprintf("%03d",start),".R",sep="")
	phybase<-'library("phybase")'
	workdir<-paste('setwd("',wd,'")',sep="")
	readtree<-paste('mytrees<-read.tree(paste("boot", sprintf("%03d",',start,'), ".rooted", sep=""))',sep="")
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
}'
	a<-paste(phybase,workdir,readtree,details,outg1,outg2,forloop,script,sep="\n")
	write.table(a,fname, row.names=F,col.names=F,quote=F)
	start<-start + increment
	end<-end + increment
}
```

Create STAR/STEAC run files

```
numbers<-seq(from=0,to=499,by=1)  #indicate start, end, and interval here

dir<-"/scratch/oliveros/zosterops/zosterops_speciestree"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=5000m,walltime=48:00:00
#PBS -M oliveros@ku.edu
#PBS -r n
#PBS -m n
#PBS -j oe"

for(i in numbers)
{
	runfile<-paste("zosterops.stst.",sprintf("%03d",i),sep="")
	pbsn<-paste("#PBS -N",runfile)
	pbsd<-paste("#PBS -d",dir)
	pbso<-paste("#PBS -o ",dir,"/",runfile,".out",sep="")
	starsteacfile<-paste("starsteac",sprintf("%03d",i),".R",sep="")
	command<-paste("R --vanilla <",starsteacfile)
	a<-paste(pbsn,pbs,pbsd,pbso,command,sep="\n")
	write.table(a,runfile, row.names=F,col.names=F,quote=F)
}
```

To submit jobs to cluster:

```
for i in zosterops.stst.???; do qsub $i; done
```

## 4. ASTRAL

Create ASTRAL run files

```
interval<-50  #indicate interval here
numbers<-seq(from=0,to=450,by=interval)  #indicate start, end here
astralcom<-"unbuffer java -jar /scratch/oliveros/Astral/astral.4.4.0.jar"
dir<-"/scratch/oliveros/zosterops/zosterops_speciestree"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=25000m,walltime=96:00:00
#PBS -M oliveros@ku.edu
#PBS -m n
#PBS -r n
#PBS -j oe
#PBS -o /dev/null"

for(i in numbers)
{
	runfile<-paste("zosterops.ast.",sprintf("%03d",i),sep="")
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
for i in zosterops.ast.???; do qsub $i; done
```

Collecting ASTRAL run times

This assumes you have saved Astral output in files named *.astral.out.
Before saving this to a file, it's good idea to pipe it to wc to check that you
have that correct number of run times.

```
cat $(echo $(ls *.astral.out)) | grep "Optimal tree" | sed -e 's/Optimal tree inferred in //' -e 's/ secs//' > astral.runtimes.txt
```

## 5 MP-EST

Create MPEST Control Files

Here's how you can create a species-allele table from Astral output:

```
cat boot000.phy.astral.out |grep "Taxa" |sed 's/Taxa: //'|sed -r 's/(\w+)/\1 1 \1\n/g' |sed 's/, //' | tr '[]' '\n'
```

Paste the appropriate species-allele table in the R script below.

```
nsim<-499  # number of replicates - 1
ntaxa<-50  # number of taxa
ngenes<-979  # number of genes

#species-allele table below
c<-"ZosSpl 1 ZosSpl
ZosRen 1 ZosRen
ZosPap 1 ZosPap
ZosMur 1 ZosMur
ZosFus 1 ZosFus
ZosGri 1 ZosGri
ZosVel 1 ZosVel
ZosKul 1 ZosKul
ZosLti 1 ZosLti
ZosTet 1 ZosTet
WooLac 1 WooLac
ZosNov 1 ZosNov
ZosMet 1 ZosMet
ZosJap 1 ZosJap
ZosSem 1 ZosSem
ZosNig 1 ZosNig
ZosCon 1 ZosCon
RukOle 1 RukOle
ZosFin 1 ZosFin
ZosAtr 1 ZosAtr
ChlEmi 1 ChlEmi
LopSup 1 LopSup
StaNig 1 StaNig
LopGoo 1 LopGoo
ZosMon 1 ZosMon
ZosMey 1 ZosMey
ZosAtf 1 ZosAtf
YuhEve 1 YuhEve
ZosVir 1 ZosVir
ZosMad 1 ZosMad
ZosSen 1 ZosSen
SpeMel 1 SpeMel
ZosPal 1 ZosPal
ZosEry 1 ZosEry
ZosEve 1 ZosEve
ZosStr 1 ZosStr
ZosChl 1 ZosChl
ZosFla 1 ZosFla
ZosUgi 1 ZosUgi
ZosSan 1 ZosSan
ZosRnl 1 ZosRnl
ZosLut 1 ZosLut
ZosCit 1 ZosCit
ZosLue 1 ZosLue
ZosHyp 1 ZosHyp
ZosExp 1 ZosExp
WooSup 1 WooSup
ZosLat 1 ZosLat
ZosLte 1 ZosLte
ZosLtr 1 ZosLtr"

for(i in 0:nsim)
{
	file<-paste("control",sprintf("%03d",i),sep="")   #filename of control file
	treefile<-paste("boot",sprintf("%03d",i),".rooted",sep="")  #filename of input rooted tree
	b<-floor(runif(1)*799736+1111)  # creates random seed number
	a<-paste(treefile,"0",b,paste(ngenes, ntaxa),c ,"0",sep="\n")  # contents of control file including num genetrees and num species
	write.table(a, file,row.names=F,col.names=F,quote=F)
}
```

Create MP-EST run files

```
interval<-1  #indicate interval here
numbers<-seq(from=0,to=499,by=interval)  #indicate start, end here
#numbers<-c(9,28,89,109,129,189,248,288,308,348,369,468,489)
dir<-"/scratch/oliveros/zosterops/zosterops_speciestree"
pbs<-"#PBS -l nodes=1:ppn=1:avx,mem=5000m,walltime=96:00:00
#PBS -M oliveros@ku.edu
#PBS -r n
#PBS -m n
#PBS -j oe"

for(i in numbers)
{
	runfile<-paste("zosterops.mpest.",sprintf("%03d",i),sep="")

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

Submit MP-EST jobs to cluster

```
for iter in zosterops.mpest.*; do qsub $iter; done
```

Monitor MP-EST progress

```
cat *.rooted.tre | grep " tree mpest" | wc
```

Collect all MP-EST trees

```
cat *.rooted.tre | grep " tree mpest" > summary
```

Create nexus file with all trees

```
cat zhead summary ztail > mpest.all.tre
```

Note:  Make zhead and ztail from nexus headers and tail of each output file.

Collect MP-EST run times

This assumes you have saved MPEST output in files named *rooted.tre
Before saving this to a file, good idea to pipe it to wc to check that you
have that correct number of run times.

```
cat $(echo $(ls *rooted.tre)) | grep "Analysis completed" | sed -e 's/\[Analysis completed in //' -e 's/ seconds\]//' > mpest.runtimes.txt
```

DENDROPY
========
sumtrees.py -o upgma.tre treefiles
sumtrees.py -f 0.95 -o star.con.tre treefiles



## 6. Running STAR, STEAC, ASTRAL, MP-EST on original data

```
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
	outgroup<-which(spname=="CyrAnn")  #indicate outgroup name here
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

outgrouptaxon<-"CyrAnn"
print(paste("Outgroup taxon", outgrouptaxon))

treestringphy<-read.tree.string("genetrees.phy",format="phylip")
treesphy<-treestringphy$tree
	
species.structure<-matrix(0,ncol=ntaxa,nrow=ntaxa)
diag(species.structure)<-1
print("Estimating STAR tree")
star<-star.sptree(treesphy, speciesname, taxaname, species.structure, outgroup=outgrouptaxon,method="nj")
starfname<-paste("star",sprintf("%03d",start),".tre", sep="")
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
java -jar /public/uce/Astral/astral.4.4.0.jar -i genetrees.phy -o genetrees.astral.tre 2>&1 | tee genetrees.phy.astral.out

mpest control
```

## Miscellaneous stuff

Getting frequencies of splits with Dendropy

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
