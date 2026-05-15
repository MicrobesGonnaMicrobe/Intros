# Introduction 

This is an introduction to methods used in some smaller projects, such as a master thesis.

- For theoretical background on 'omics, I recommend having a look at this webinar series Microbial 'Omics for Beginners from the Meren lab:
https://youtube.com/playlist?list=PL7133RHfhW-MwCLz-c2DZxAmtoHipqBcL&si=lvaZ9sGuyEsVXmVP
- Explanations of couple of 'omics concepts: https://anvio.org/vocabulary/

## Content
  - [Intro to working with command line / Basic Unix (bioinformatics)](#intro-to-working-with-command-line-basic-unix-bioinformatics)
  - [Setting up the Windows subsystem for Linux (WSL)](#setting-up-the-windows-subsystem-for-linux-wsl)
  - [Downloading programs and setting up conda environments](#downloading-programs-and-setting-up-conda-environments)
  - [Useful pieces of code](#useful-pieces-of-code)
  - [Downloading genomes on NCBI + check GTDB](#downloading-genomes-on-ncbi-check-gtdb)
  - [Anvio](#anvio)
    - [Importing genomes to Anvio](#importing-genomes-to-anvio)
    - [Annotations](#annotations)
    - [Extracting 16S rRNA genes](#extracting-16s-rrna-genes)
  - [Phylogenomics](#phylogenomics)
    - [Make a list of single copy genes to use](#make-a-list-of-single-copy-genes-to-use)
    - [Get amino acid sequences](#get-amino-acid-sequences)
    - [Align individual sequences with mafft](#align-individual-sequences-with-mafft)
    - [Trimming](#trimming)
    - [Concatenate](#concatenate)
    - [Build tree](#build-tree)
  - [Genome similarity](#genome-similarity)
    - [Average nucleotide identity (ANI)](#average-nucleotide-identity-ani)
    - [Average aminoacid identity (AAI)](#average-aminoacid-identity-aai)
  - [Gene annotation](#gene-annotation)
  - [Pangenomics](#pangenomics)
  - [Plotting in R](#plotting-in-r)

## Intro to working with command line / Basic Unix (bioinformatics)
Learning the basics of the command line is valuable so you can use the most commonly used bioinformatics tools. 
To start and get a feeling of it, have a look at this video series of Basic Unix lessons for bioinformatics: https://youtube.com/playlist?list=PL4K-daRUS2A-knoIZVV5SJkEkZhqfdHVB&si=nV1uEyOjEV5asOSH
The written part of the lessons can be found here: https://astrobiomike.github.io/unix/unix-intro

In case you want to play around with command line, here is a very nice UNIX tutorial with example files: https://ndombrowski.github.io/Unix_tutorial/

It might feel a bit overwhelming at first, but don’t get discouraged.

## Setting up the Windows subsystem for Linux (WSL)
Most people first dipping their toes into bioinformatics analyses will have a Windows computer.
To perform the following analyses, you will need to install the Windows subsystem for Linux: https://apps.microsoft.com/detail/9nz3klhxdjp5?hl=en-us&gl=US

For easier management, I would also recommend installing Windows Terminal and use this to work in the Ubuntu tab: https://apps.microsoft.com/detail/9n0dx20hk701?hl=en-us&gl=US

## Downloading programs and setting up conda environments
To understand what is conda and how we use it to manage programs, have a look at this conda tutorial: https://youtu.be/ukgogSsW7oI?si=qI2OWrUQxIDwJrD3

1. In WSL, install Miniconda, following these instructions: [https://www.anaconda.com/docs/getting-started/miniconda/install/overview](https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install)
Close and re-open your terminal window for the installation to fully take effect.

2. Create a conda environment and activate it

3. Install conda packages (like the anvio program, see below).

## Useful pieces of code
Ask the intrenet for code advice if you get stuck.
Here is a list of some useful commands that might come in handy: https://github.com/MicrobesGonnaMicrobe/tools

## Downloading genomes on NCBI + check GTDB

Don't forget to gather metadata information!

## Anvio
For some analyses we will be using anvi'o, a platform for microbial 'omics. 
- Instructions to install here: https://anvio.org/install/
Video on what it can be used for: https://merenlab.org/2016/11/08/pangenomics-v2/

### Importing genomes to Anvio

#### Re-formatting input FASTA
```
mkdir anvio_reformatfasta
```

For loop:
```
for i in *.fa; do anvi-script-reformat-fasta --simplify-names --seq-type NT -r "report_simplify_names"${i}.txt  -o anvio_reformatfasta/${i} ${i}; done
```

#### Creating an anvi’o contigs database
With reformated fasta.

```
cd anvio_reformatfasta
mkdir contigs_database
for i in *.fa; do anvi-gen-contigs-database -f ${i} -o contigs_database/${i%.fa}.db --num-threads 6; done
```

#### How to make external_selected_genomes.txt
```
echo -e "name\tcontigs_db_path" > external_selected_genomes.txt
for i in *.db; do echo -e "${i%.db}\t$PWD/$i"; done >> external_selected_genomes.txt
```

### Annotations

#### anvi-run-hmms
The results are stored in the contigs database.
```
for i in *.db; do anvi-run-hmms -c ${i} --num-threads 20; done 
```

#### anvi-run-kegg-kofams
Essentially, this program uses the KEGG database to annotate functions and metabolic pathways in a contigs-db.
```
for i in *.db; do anvi-run-kegg-kofams -c ${i} --num-threads 20; done
```

### Extracting 16S rRNA genes
```
anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt --hmm-source Ribosomal_RNA_16S --list-available-gene-names
```

## Phylogenomics
Have a look at this video explaining phylogenomics: https://www.youtube.com/live/hfHu8Lnwgzs?si=6uKKPtg2f7JXcL08

Overview of the workflow:
- Anvio
- MAFFT
- AliView
- trimal
- https://github.com/nylander/catfasta2phyml
- IQ-TREE
- ITOL: https://itol.embl.de/help.cgi#annot

### Make a list of single copy genes to use
List available HMM sources in the contigs database
```bash
anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt --list-hmm-sources
```
### Get amino acid sequences
```bash
anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt --hmm-source Bacteria_71 --list-available-gene-names
```

Decide which set of markers to use
- only select the single copy marker genes (check with a matrix if they are in duplicates)
- select the ones that are mostly present in the genomes
- more good quality marker genes, the better

Not always needed, only when deciding which markers to use: To check presence of marker genes in all genomes (matrix):
```bash
anvi-script-gen-hmm-hits-matrix-across-genomes --external-genomes external_selected_genomes.txt --hmm-source Bacteria_71 -o Zeta_Bacteria71_markers_matrix.txt
```

Choose a subset of those genes and get them in one fasta file (without alignment and concatenation)
```bash
anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt -o Zetaproteobacteria_selectedmarkers_Bacteria71.fa --hmm-source Bacteria_71 --gene-names Zetaproteobacteria_selectedmarkers_Bacteria71.txt --return-best-hit --get-aa-sequences
```
For --gene-names, you can either provide a list of selected markers (for example --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3) or a text file including names of selected markers (--gene-names Zetaproteobacteria_selectedmarkers_Bacteria71.txt).


Not always needed, only when deciding which markers to use: Make separate trees for all separate proteins (to check if monophyletic/good marker genes)
```bash
anvi-gen-phylogenomic-tree -f Zetas_ribosomal_L15.fa -o Zetas_ribosomal_L15_tree.txt
```

Split multifasta markers from anvio into single marker files
- The file with all markers is divided based on the header before the "___" symbol into one file per marker:

```bash
seqkit split -i --id-regexp "^(\\S+)\___\s?" Zeta_ribosomal_markers20_proteins.fa
```

### Align individual sequences with mafft
* `MAFFT L-INS-i v7.397`
```bash
mkdir individual_mafft
for i in *.fa; do mafft-linsi $i | awk 'BEGIN{FS=":|[|]"}{if(/^>/){print ">"$2}else{print $0}}' > individual_mafft/${i%.fa}_mafft.fa; done
```

- Inspect the alignment manually
* `AliView v1.26`

### Trimming
* `trimAl v1.4.rev15`

Before using trimal, remove spaces in >fasta headers, so that trimal does not remove the taxonomic classification reported in the header
```bash
sed -i 's/ /_/' *mafft.fa
sed -i 's/:/_/' *mafft.fa

mkdir trimal
for i in *mafft.fa; do trimal -in $i -gt 0.5 -cons 60 |cut -f 1 -d ' ' > trimal/${i%.fa}_trimal.fa; done

-gt 0.5 -cons 60
#### Removes all positions in the alignment with gaps in 50% or more of the sequences, unless this leaves less than 60%. In such case, print the 60% best (with less gaps) positions.
```

### Concatenate
* `catfasta2.phyml v07.04.20`: (https://github.com/nylander/catfasta2phyml)
```bash
catfasta2phyml -v -c -f *mafft_trimal.fa > Zetaproteobacteria_concat_mafft_trimal.fa
```

### Build tree
* `IQ-TREE v2.0.3`: https://github.com/Cibiv/IQ-TREE

Choose the best substitution model (Best-fit model)
- "By default, substitution models are not included in these tests. If we want to test them we have to add them. Generally, it is recommended to include them in the test and the following selection would be quite comprehensive for testing models."

```bash
iqtree -s Zetaproteobacteria_concat_mafft_trimal.fa -m MFP -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60,LG+C10+R+F,LG+C20+R+F,LG+C30+R+F,LG+C40+R+F,LG+C50+R+F,LG+C60+R+F -v -nt 4
```

#### Bootstrapping
Various methods allow to assess the confidence in branching patterns or branch supports.
Run another tree with a selected best-fit model (in the case below it was LG+F+R7) and bootstrapping.

##### Ultrafast bootstrapping and SH-like approximate likelihood ratio test (SH-aLRT)
```bash
iqtree -s Zetaproteobacteria_concat_mafft_trimal.fa -m LG+F+R7 -v -bb 1000 -alrt 1000 -nt 4 -pre Zetaproteobacteria_bb1000alrt1000_LG_F_R7
```
Numbers in parentheses are SH-aLRT support (%) / ultrafast bootstrap support (%).
From IQ-TREE page: "One would typically start to rely on the clade if its SH-aLRT >= 80% and UFboot >= 95%."

##### The standard non-parametric bootstrapping with 1000 replicates
This can take very long time, so only do this when everything else is ready and this much use of computational resources is needed (always run on the server).
```bash
iqtree -s Zetaproteobacteria_concat_mafft_trimal.fa -m LG+F+R7 -b 1000 -nt 4 -pre Zetaproteobacteria_b1000_LG_F_R7
```

#### Tree visualisation and annotation
Visualise the tree by importing to iTOL and use templates for tree annotation: https://itol.embl.de/help.cgi#annot

## Genome similarity

### Average nucleotide identity (ANI)
```
anvi-compute-genome-similarity -e external_selected_genomes.txt -o DPANN_HighMediumQC_ANI --program pyANI --method ANIb -T 25
```
In the resulting graphic, define thresholds (for example yellow for under 95% ANI, red for over 95% ANII (species threshold)).

### Average aminoacid identity (AAI)

One program option: ezAAI
Follow the tutorial on this page: https://endixk.github.io/ezaai/

In the resulting graphic, define thresholds (for example yellow for under 65% AAI, orange for over 65% AAI (genus threshold)).

## Gene annotation
A quick introduction to annotation in this video (watch the first 27 min) - https://www.youtube.com/watch?v=OaIgO7VKpS8 
Introduction to annotation with KEGG (watch the first 6 min) - https://www.youtube.com/watch?v=1v4UzjE7K2g

Annotation of several genomes using several databases can be computationally demanding. You will most likely receive already annotated genomes, following this annotation protocol: https://ndombrowski.github.io/Annotation_workflow/

It is great to screen through the annotation files manually to get a feel for what kind of gene annotations you find in your genomes.
Not all annotations might be accurate. Be aware of length of the predicted gene, position on contig, and genes around it (downstream and upstream).

To visualise presence / absence of specific genes, you can visualise them in a heatmap
- R visualisation tutorial: https://ndombrowski.github.io/R_tutorial/#working-with-annotation-data

## Pangenomics
There is an anvi'o workflow for making pangenomes (usually focusing on one genus) - https://merenlab.org/2016/11/08/pangenomics-v2/

## Plotting in R

Here are some examples of code for plots:
- Relative abundance (Relative_abundances_plot_MAGs.r)

Background on stacked barcharts: https://www.youtube.com/watch?v=NVym44SdcaE and https://www.youtube.com/watch?v=w4X3o6MQjVA

- Heatmap for ANI / AAI

Background on heatmaps: https://www.youtube.com/watch?v=-4sGamyqGAQ

- Metabolism presence/absence
