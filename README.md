# Intros

## Content

##Anvio: Importing genomes
- https://anvio.org/

###Re-formatting input FASTA
```
mkdir anvio_reformatfasta
```

For loop:
```
for i in *.fa; do anvi-script-reformat-fasta --simplify-names --seq-type NT -r "report_simplify_names"${i}.txt  -o anvio_reformatfasta/${i} ${i}; done
```

###Creating an anvi’o contigs database
With reformated fasta.

```
cd anvio_reformatfasta
mkdir contigs_database
for i in *.fa; do anvi-gen-contigs-database -f ${i} -o contigs_database/${i%.fa}.db --num-threads 6; done
```

