# **_brevirostris_ Group Reticulation from UCE Data**

Ahoy.

This pipeline is the workflow used to estimate reticulation in the _Plestiodon brevirostris_ group using 3,282 Ultra-Conserved Elements (UCEs) from [Bryson Jr. et al. 2017](https://doi.org/10.1111/jbi.12989).

To do so, we will be mapping the UCEs to a chromosome-level reference genome and then concatenating UCEs that are within 200,000bps of eachother to approximate 'cogenic' UCEs. Then, we will take the 'psuedo-genes' and continue the analysis with those larger than 1,000bps. With this concatenation scheme we are creating more informative partitions that attempt to be biologically relevant.

## **Data Preperation**

In order to map UCEs on concatenated plestiodonUCE.txt file in the supplementary data of Bryson Jr. et al 2017 needs to be split into each of the 59 individuals
```
counter=0; cat plestiodonUCE.txt | while read LINE; do ((counter++)); echo $LINE > "$counter.txt"; done
```

Then they are manipulated into .fasta format
```
#add carrot in front of name
#!/bin/bash
for i in {2..60}
do
	sed -i '' '1s/^/>/' ${i}.txt
done
```
```
#.txt to .fasta
#!/bin/bash
for i in {2..60}
do
	awk '{ print }' ${i}.txt > ${i}.fas
done
```
OK, so I needed to add a line break after the >sampleID part but I just did it manually. I'm sure you could use sed but I was tired of sed

Now, we split each of the species fastas into UCEs using the partition.txt file from the supplementary material
```
#This script was written by Nathan Whelan.  

library(Biostrings)
specimen <- as.vector(c(2:60))

split <- function(x) {
setwd(paste0("/Users/jonhoff/Documents/AMNH/Thesis_Research/Phylo/Brevirostris_Project/all_ind_fas/", x)); ##Change as needed

##The file name should be a supermatrix in phylip format. Must modify variable names as needed
sequenceData<-unmasked(readAAMultipleAlignment(paste0(x,".fas"),format="fasta"))

##Gene list should be in PartitionFinder format (i.e. GENE_NAME = startPostition-stopPosition;)
##Semicolon at end of each line in gene list is not necessary
table<-read.table("/Users/jonhoff/Documents/AMNH/Thesis_Research/Phylo/Brevirostris_Project/3282partitions_PF-format.txt")
list<-data.frame(do.call('rbind', strsplit(as.character(table$V3),split="-",fixed=TRUE)))
list$X2<-gsub(";","",list$X2)
startVar<-as.data.frame(list$X1)
endVar<-as.data.frame(list$X2)
numberRows<-as.integer(nrow(startVar))

#This for loop will go through and split a supermatrix as specified by gene list
#A single alignment in fasta format will be created for each gene/partition
##See bash script in repository for removing taxa that was not sampled for any given gene output by this for loop
for(n in 1:numberRows){
  subset <- AAStringSet(sequenceData, start=as.integer(as.vector(startVar[n,1])), end=as.integer(as.vector(endVar[n,1])))
  LINE1=toString(n)
  NAME=paste(LINE1,"_UCE.fas",sep="",collapse=NULL)
  writeXStringSet(x=subset, filepath = NAME, format="fasta") #This will give a unique name
  print(n)
  print(subset)
}

}

lapply(specimen, split)
```

Here's a step that annoyed me: removing spaces after the names for each UCE for each species
```
#!/bin/bash
for i in {2..60}
do
 cd all_ind_fas/${i}
	for x in {1..3282}
	do
		sed "s/[[:space:]]*$//g" ${x}_UCE.fas > ${x}_ns.fas
	done
 cd ..
 cd ..
done
```

Last step of the preparation is to rename the UCE fastas to inlcude the UCE number, rather than just the specimen name
```
library(ape)
library(phytools)
library(phylotools)

setwd("/Users/jonhoff/Documents/AMNH/Thesis_Research/Phylo/Brevirostris_Project/all_ind_fas")

#vector list of individuals
critter <- as.vector(c("xxxxxxx","mx41_plb","mx43_pll","mx44_plb","mx46_pll","mx49_plb","mx50_pll","mxh10_plb","mxh11_plb","mxh12_plb","mxh13_plb","mxh14_plb","mxh15_plb","mxh16_plb","mxh17_plc","mxh187_plb","mxh188_pll","mxh189_plpv","mxh18_plb","mxh192_plpv","mxh193_plb","mxh194_plb","mxh19_plb","mxh1_plb","mxh20_plb","mxh21_plpv","mxh22_plpv","mxh25_plpvia","mxh26_plpvia","mxh27_plo","mxh28_plo","mxh29_pls","mxh32_plb","mxh33_plb","mxh34_plb","mxh35_pll","mxh36_pll","mxh37_pll","mxh38_pll","mxh3_pld","mxh40_pld","mxh42_plb","mxh44_out","mxh45_plb","mxh46_plb","mxh4_plpvia","mxh50_plpv","mxh52_plcol","mxh53_pld","mxh54_pll","mxh55_pll","mxh56_pll","mxh57_pll","mxh58_pll","mxh59_pll","mxh5_plc","mxh6_pld","mxh7_plb","mxh8_plb","mxh9_plb"))

#function to rename 
rename <- function(number) {
for (i in 1:3282) {
  df <- data.frame(og_name=c(paste0(critter[number])), new_name=c(paste0(critter[number],"_", i)))
  rename.fasta(paste0(number,"/",i,"_UCE.fas"), df, outfile=paste0(number,"/",i, "_final.fas"))
}
}

#vector to tell function 2-60
dog <- as.vector(c(2:60))

lapply(dog, rename)
```

Now we're ready to blast!

## **Blasting UCEs against closest chromosome-level reference genome _Tiliqua scincoides_**

First we need to install blast with Homebrew (we're using MacOs)
```
brew install blast
```

Then, we need to create a custom database to blast against using a downloaded _Tiliqua scincoides_ genome
```
makeblastdb -in GCA_035046505.1_rTilSci1.hap2_genomic.fna -dbtype nucl -out tiliqua
```

Ok, time to blast! Blasting will tell us where these UCEs are located on the genome. For the blasting and some of the subsequent analyses, I used multiple specimens to compare the outputs in order to decide if the final partitions for one specimen will be representative of all specimens.
```
#!/bin/bash
for i in {1..3282}
do 	
	blastn -db tiliqua -query final_UCE_fastas/${i}_final.fas -out blast_output_tiliqua/UCE_${i}.txt -strand plus -outfmt 6
done
```

Now we sort the blast results to isolate the best hit for each UCE
```
#!/bin/bash
for i in {1..3282}
do 
	sort -k1,1 -k12,12nr -k11,11n  blast_output_tiliqua/UCE_${i}.txt | sort -u -k1,1 --merge > blast_output_tiliqua/besthit/bUCE_${i}.txt
done
```

Finally, we concatenate all of the blast ouputs from all the UCEs back into one file
```
setwd("~/Documents/AMNH/Thesis_Research/Phylo/Brevirostris_Project/blast_output_tiliqua/besthit")


list.files()
uces <- list.files()

datalist <- lapply(uces, function(x)read.table(x, header=F))
datafr <- do.call("rbind", datalist)

write.csv(datafr, file = "all_UCE_blast.csv", row.names = FALSE)
```
## **Creating psuedo-genes**



## **Infer Gene Trees with IQ-TREE**

The first step in this project is to take the aligned UCE data and use [IQ-TREE](https://github.com/iqtree/iqtree2) to infer unrooted gene trees for each of the 3,282 UCEs
```

```

## **Infer Species Tree with ASTRALIII**

Now that we have gene trees from IQtree2, we can use [ASTRALIII](https://github.com/smirarab/ASTRAL) to infer an unrooted species tree. ASTRALIII is a is a tool for estimating an unrooted species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under the multi-species coalescent model, allowing for appropriate inference of Incomplete Lineage Sorting (ILS) and introgression. 
```

```

## **Visualize Gene Tree Discordance**


## **Estimate Introgression with QuIBL**
