# Epigenomics Final Task 

By: Thijs Dormans

Date: 19-02-2024

# Introduction

This practical is part of the course in epigenomics. This practical builds on the handout from the epigonimics_uvic github of Beatrice Borsari (https://github.com/bborsari/epigenomics_uvic/wiki). Several steps from this github (section 2 &3) have been followed whereby EN‐TEx ChIP‐seq data has been obtained and downstream analyses have been performed. For this practical we focus specifically on sections 4 & 5 which has the following  objectives:

**Section 4: EN‐TEx ATAC‐seq data: downstream analyses:** 

In this section we obtain  EN-TEx ATAC-seq data from donor ENCDO451RUA on the ENCODE portal (https://www.encodeproject.org). Next we run an intersection analysis using BEDTools to find the number of ATAC-seq peaks that intersect promotor regions and the number of ATAC-seq peaks that fall outside of gene coordinates.

**Section 5: Distal regulatory activity:**

In this section we will utilize the ATAC-seq peaks that have been obtained in section 4. We will explore regions that are  flanked by H3K27ac and H3K4me1. Additionally, we will perform numerous tasks that allow additional insights into regulatory elements.

# Setup
First we run the docker image, mount it and enter the container:
```{r setup, echo=FALSE }
knitr::opts_chunk$set(eval = FALSE)
library(rmdformats)
```

```{bash}
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```
# Tasks of hands-on section 4

## 4.1 Set up folders

We create numerous folders to store our outouts.
```{bash}
cd ./ATAC-seq
mkdir analyses
mkdir -p data/bigBed.files data/bigWig.files
mkdir analyses/peaks.analysis 
mkdir data/bed.files
mkdir annotation
```

## 4.2 Retrieve required bigBed files, check MD5 and obtain Gencode annotation files.

After navigating the ENCODE database we obtain the metadata for our selected files. 
```{bash}
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=ATAC-seq&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&type=Experiment"
```

Let's take a look at the file. We link head and awk to obtain information about the columns.
```{bash}
head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'
```

The information we need can be found in column 1 (file accession), for tissue column 11 (Biosample_term_name) and column 23 (Experiment_target)

Let's start parsing the metadata file:

```{bash}
# filter bigBed_Narrowpeak, pseudoreplicated_peaks and accession GRChc38
grep -F "bigBed_narrowPeak" metadata.tsv|\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |
# Now we extract the desired columns to the peaks file
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
# We perform primary sort on tissue column and secondary sort on file accession and remove duplicates (-u) and store the output in a file
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt
```

After checking the file, we utlize the obtained IDs to get the files of interest.
```{bash}
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

Let's check out the integrity of the files by checking MD5 hash using *md5sum*. First we extract the original MD5 hash from the metadata. Next generate MD5 hashes for our obtained files. Lastly, we check that outputs are similar. 
```{bash}
for file_type in bigBed; do

  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
```

We received no input back, indicating that there are no differences. If you want to visually verify, you can do accordiingly:
```{bash}
less data/bigBed.files/md5sum.txt
```
```{latex}
ENCFF762IFP     f6a97407b6ba4697108e74451fb3eaf4        f6a97407b6ba4697108e74451fb3eaf4
ENCFF287UHP     46f2ae76779da5be7de09b63d5c2ceb9        46f2ae76779da5be7de09b63d5c2ceb9
```

As you can see the MD5 hashes are identical.

Now that we have our bigBed files, we need our annotation files. As these files are identical to the procedure followed during class for the ChIP-seq hands-on, we will copy these to be efficient!

```{bash}
cp ../ChIP-seq/annotation/*v24*.bed ./annotation
```

## 4.3 Intersection analysis using BEDTOOLS

First let's convert our bigBed files to bed format.

```{bash}
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

We will now perform the intersection analyses for both the tissues. We extract the first two columns and set their values as filename and as tissue variable respectively in order to calculate the intersection. We store this in a txt file.

#### Number of peaks that intersect promotor region
```{bash}
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -u |\
  sort -u > analyses/peaks.analysis/ATACpeaks.within.promoters."$tissue".txt
done
```

The number of peaks at the promotor region for both tissues can be found  using wc -l.

```{bash}
wc -l analyses/peaks.analysis/*.txt
```
```{latex}
47871 analyses/peaks.analysis/ATACpeaks.within.promoters.sigmoid_colon.txt
44749 analyses/peaks.analysis/ATACpeaks.within.promoters.stomach.txt
92620 total
```

We find 47871 ATAC-seq peaks in promotor regions for our sigmoid colon tissue and 44749 for the stomach tissue.

#### Number of peaks that fall outside gene coordinates

Since we need this file for the next assignment we will create a BED file straight away rather than making a txt file. We set -v to focus on areas that do NOT overlap.

```{bash}
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
	bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.gene.body.bed -v |\
	sort -u > analyses/peaks.analysis/ATACpeaks.outside.gene."$tissue".bed
done
```

Let's check the number of peaks outside of the gene coordinates.
```{bash}
wc -l analyses/peaks.analysis/*.bed
```
```{latex}
37035 analyses/peaks.analysis/ATACpeaks.outside.gene.sigmoid_colon.bed
34537 analyses/peaks.analysis/ATACpeaks.outside.gene.stomach.bed
71572 total
```

We find 37035 ATAC-seq peaks in outside of gene coordinates for sigmoid colon tissue and 34537 for stomach tissue.

The above generated .bed files containing peaks outside of the gene coordinates will serve as the basis for the following steps.

# Tasks of hands-on section 5

## Task 5.1

First we create the required folders.
```{bash}
mkdir ../regulatory_elements
cd ../regulatory_elements
mkdir data
mkdir analyses
mkdir analyses/peaks.analysis
mkdir data/bed.files
mkdir data/bigBed.files
```

## Task 5.2

Next we need to obtain the required files. 

For H3K4me3:
```{bash}
grep -F H3K4me1 ../ChIP-seq/metadata.tsv|\
grep -F "bigBed_narrowPeak"|\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/H3K4me1.bigBed.peaks.ids.txt

cut -f1 analyses/H3K4me1.bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

For H3K27ac: 
```{bash}
grep -F H3K27ac ../ChIP-seq/metadata.tsv|\
grep -F "bigBed_narrowPeak"|\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/H3K27ac.bigBed.peaks.ids.txt

cut -f1 analyses/H3K27ac.bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

Just like before we need to perform the hash check.
```{bash}
for file_type in bigBed; do

  ../bin/selectRows.sh <(cut -f1 analyses/*"$file_type".peaks.ids.txt) ../ChIP-seq/metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
```

We again get no input as a return which is a good sign. If you want, perform a visual check using:
```{bash}
less data/bigBed.files/md5sum.txt
```
```{latex}
ENCFF844XRN     de679228721fb4055aa1f657c77c21a6        de679228721fb4055aa1f657c77c21a6
ENCFF724ZOF     c87fefbf41de3d291fa1d340a26627f5        c87fefbf41de3d291fa1d340a26627f5
ENCFF977LBD     be29636550527e36c4755ea036531e75        be29636550527e36c4755ea036531e75
ENCFF872UHN     2207b7b3378df7776e7ecdc2aa1a5de0        2207b7b3378df7776e7ecdc2aa1a5de0
```
Again they are identical.

Next we convert our files from bigBed to bed. 
H3K4me1:
```{bash}
cut -f1 analyses/H3K4me1.bigBed.peaks.ids.txt|\
while read filename; do
        bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```
H3K27ac:
```{bash}
cut -f1 analyses/H3K27ac.bigBed.peaks.ids.txt|\
while read filename; do
        bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

As distal regulatory regions are located outside of protein coding areas so we need to perform intersecting steps. We will perform the following procedures for both tissues:

* A:Find intersect between H3K4me1 and ATAC-peaks outside of genes and create bed file
* B:Find intersect between H3K4me1 and ATAC-peaks outside of genes and create bed file
* C:Find intersect between H3K27ac two above generated bed files.

TEST!!

----------------------------------------------------------------------
A: Find intersect between H3K4me1 and ATAC-peaks outside of genes and create .bed file.


```{bash}
cut -f-2 analyses/H3K4me1.bigBed.peaks.ids.txt |
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b ../ATAC-seq/analyses/peaks.analysis/ATACpeaks.outside.gene."$tissue".bed -u > analyses/peaks.analysis/"$tissue".H3K4me1.outside.bed
done
```


B:Find intersect between H3K27ac and ATAC-peaks outside of genes and create .bed file.
```{bash}
cut -f-2 analyses/H3K27ac.bigBed.peaks.ids.txt |
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b ../ATAC-seq/analyses/peaks.analysis/ATACpeaks.outside.gene."$tissue".bed -u > analyses/peaks.analysis/"$tissue".H3K27ac.outside.bed
done
```

C:Find intersect between H3K27ac two above generated bed files and create final .bed file.
```{bash}
for tissue stomach sigmoid_colon; do
  cut -f-2 analyses/H3K27ac.bigBed.peaks.ids.txt |
  while read filename tissue; do 
    bedtools intersect -a analyses/peaks.analysis/"$tissue".H3K27ac.outside.bed -b analyses/peaks.analysis/"$tissue".H3K4me1.outside.bed -u > analyses/peaks.analysis/overlap.outside.histone."$tissue".bed
done
```

Obtain the counts.

```{bash}
wc -l analyses/peaks.analysis/*histone*.bed

wc -l analyses/peaks.analysis/*hostone*.bed
```
```{latex}
  7367 analyses/peaks.analysis/overlap.outside.histone.sigmoid_colon.bed
  4342 analyses/peaks.analysis/overlap.outside.histone.stomach.bed
 11709 total
```


## Task 5.3

Task description: Focus on regulatory elements that are located on chromosome 1, and generate a file *regulatory.elements.starts.tsv* that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.

Let's follow the reqiured steps for each tissue. But first we need to know which columns contain this info. Lets take a look at the stomach file to explore.

```{bash}
head -1 analyses/peaks.analysis/overlap.outside.histone.stomach.bed | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'
```

```{latex}
chr1    1
1067466 2
1070082 3
Peak_8834       4
813     5
.       6
15.55870        7
81.38033        8
78.60110        9
1035    10
```


The chromosome information is located in column 1, the start coordinate in column 2 and the name of the peak is in column 4.
Let's go ahead and obtain the infomation we need fo rboth tissues.

```{bash}
for tissue in sigmoid_colon stomach; do
  awk 'BEGIN{FS=OFS="\t"} $1=="chr1" {print $4,$2}' analyses/peaks.analysis/overlap.outside.histone."$tissue".bed > analyses/peaks.analysis/"$tissue".regulatory_elements.starts.tsv
done
```

```{bash}
wc -l analyses/peaks.analysis/*.regulatory_elements.starts.tsv
```
```{latex}
795 analyses/peaks.analysis/sigmoid_colon.regulatory_elements.starts.tsv
528 analyses/peaks.analysis/stomach.regulatory_elements.starts.tsv
1323 total
```

## Task 5.4

Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here, prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column

```{bash}
awk 'BEGIN{FS=OFS="\t"} $1=="chr1" {if ($6=="+"){start=$2} else {start=$3}; print $4, start}' ../ATAC-seq/annotation/gencode.v24.protein.coding.gene.body.bed > analyses/gene.starts.tsv
```

If we want we we can check out the file:
```{bash}
head -5 analyses/gene.starts.tsv
```
```{latex}
ENSG00000186092.4       69090
ENSG00000279928.1       182392
ENSG00000279457.3       200322
ENSG00000278566.1       451678
ENSG00000273547.1       686654
```


## Task 5.5

The scripted that was provided was altered as follows:

```{latex}
#!/usr/bin/env python


#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser


#*****************
# OPTION PARSING *
#*****************

parser = OptionParser()
parser.add_option("-i", "--input", dest="input")
parser.add_option("-s", "--start", dest="start")
options, args = parser.parse_args()

open_input = open(options.input)
enhancer_start = int(options.start)


#********
# BEGIN *
#********

x=1000000 # set maximum distance to 1 Mb
selectedGene="" # initialize the gene as empty
selectedGeneStart=0 # initialize the start coordinate of the gene as empty

for line in open_input.readlines(): # for each line in the input file
	gene, y = line.strip().split('\t') # split the line into two columns based on a tab 
	position = int(y) # define a variable called position that correspond to the integer of the start of the gene
	diff_pos_start = abs(position-enhancer_start)# compute the absolute value of the difference between position and enhancer_start

	if diff_pos_start < x: # if this absolute value is lower than x
		x = diff_pos_start# this value will now be your current x
		selectedGene = gene# save gene as selectedGene
		selectedGeneStart = position # save position as selectedGeneStart

print "\t".join([selectedGene, str(selectedGeneStart), str(x)])
```

Now lets test the script and see if we get the desired output (ENSG00000187642.9	982093 2093)
```{bash}
python ../bin/get.distance.py --input analyses/gene.starts.tsv --start 980000
```
```{latex}
ENSG00000187642.9       982093  2093
```  

It works!

## Task 5.6

For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above.

```{bash}
for tissue in stomach sigmoid_colon; do
  cat analyses/peaks.analysis/"$tissue".regulatory_elements.starts.tsv |
  while read element start; do
    python ../bin/get.distance.py --input analyses/gene.starts.tsv --start "$start"
  done > analyses/peaks.analysis/"$tissue".gene_distances.tsv
done
```

## Task 5.7

Now using R we can caluclate the mean and median of the stored distances. Let's write a small function
```{bash}
R
compute_mean_median <- function(file) {
  distances <- read.delim(file, header = FALSE)
  cat("Mean Distance:", mean(distances[,3]), "\n")
  cat("Median Distance:", median(distances[,3]), "\n")
}
```

Next we calculate mean and median using this function.
For stomach:

```{bash}
compute_mean_median("analyses/peaks.analysis/stomach.gene_distances.tsv")
```
```{latex}
Mean Distance: 47878.12
Median Distance: 27245
```

For sigmoid colon:
```{bash}
compute_mean_median("analyses/peaks.analysis/sigmoid_colon.gene_distances.tsv")
```
```{latex}
Mean Distance: 76026.47
Median Distance: 36609
```
Summary:

**Stomach**

* Mean Distance: 47878.12
* Median Distance: 27245

**Sigmoid colon**

* Mean Distance: 76026.47
* Median Distance: 36609