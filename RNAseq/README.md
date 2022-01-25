

# Quality control
using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.8


# trimming:
We first trimmed the 10 bp in the 5' and the last 2 bp in the 3'. The length of females samples are longer. we used [fastx-toolkits](http://hannonlab.cshl.edu/fastx_toolkit/index.html) 0.0.13

### for males samples
```
fastx_trimmer \
-f 11 \
-l 72 \
-i samples/14.fastq \
-o trim1/14.fastq
``` 
### for females samples
```
fastx_trimmer \
-f 11 \
-l 77 \
-i samples/16.fastq \
-o trim1/16.fastq
```
Then, we trimes the Illumina adaptors and low quality read using [trim-galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) 0.6.7. and [cutadapt](https://cutadapt.readthedocs.io/en/stable/#) 3.5
with defualt parameters, and re-run FastQC.
```
trim_galore \
-o trim2/16.fastq \
--fastqc \
trim1/16.fastq 
```



# mapping:
We mapped the trimed squence to the [N.furzeri](https://www.ncbi.nlm.nih.gov/genome/?term=furzeri) genome using [STAR](https://github.com/alexdobin/STAR) 2.7.9a
```
STAR \
--runThreadN 16 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFileNamePrefix STAR_output/16 \
--genomeDir RefGenomes/Nfu_20140520/STAR/ \
--readFilesIn trim2/16.fastq
```

We perform differential gene expression as a function of age, genotype, and the interaction between age and genotype was performed using the [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) 3.32.1 and 
Enriched Gene Ontology (GO) terms associated with transcripts levels were identified using Gene Set Enrichment Analysis (GSEA) implemented in R package [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) 3.18.1 

The analysis done by [Tehila Atlan](https://github.com/tehila-atlan) (Harel-lab)
RNA editing analyzed by Kobi Shapira from [Erez Levanon's lab](https://www.levanonlab.com/)
