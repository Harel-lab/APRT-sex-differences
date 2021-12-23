# Create a conda virtual envrionment

This is the main envrionment for our analysis.

```bash
conda env create -f Killifish-APRT-RNAEditing.yml
```

# Genomic annotations 

We assume you have the following files:
1. `genome.fasta`
2. `transcriptome.fasta`
3. `genomic.gff`
4. `CDS.gff`
5. `repeats.gff`  

In our case, files 1-3 correspond to (1) `GCF_001465895.1_Nfu_20140520_genomic.fna`, (2) `GCF_001465895.1_Nfu_20140520_rna.fna`, and (3) `GCF_001465895.1_Nfu_20140520_genomic.gff` from the [Turquoise killifish RefSeq annotations](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Nothobranchius_furzeri/latest_assembly_versions/GCF_001465895.1_Nfu_20140520/), unzipped.
The `CDS.gff` is easily created by selecting all `CDS` features from the `genomic.gff` file.
The `repeats.gff` is the aligment of repetitive elements to the `genome.fasta`. In our case, the repetitive elements were annotated de-novo and aligned using [EDTA](https://github.com/oushujun/EDTA).
In this manual, these files are assumed to be found in the `Data/Annotations` dir.

# RNA data

## Quality control

Certain filtartion steps are needed to eliminate common biases in A-to-I RNA editing.

### Deduplicating & trimming

We used [PRINSEQ-lite](http://prinseq.sourceforge.net/) 0.20.4 with [perl](https://www.softwarecollections.org/en/scls/rhscl/perl516/) 5.16 to remove exact duplicates, retain reads no shorter than 62 bp, and trim reads longer than 63 bp. Duplicate RNA reads – defined as reads with the same sequence on the same strand or with the reverse-compliment sequence on the opposite strand – can result from PCR cycles conducted before sequencing. Using only reads with similar lengths allowed us to compare the RNA Editing Index of different samples.
The command looked like this:

```bash
perl prinseq-lite.pl \
-fastq $sample.fastq \
-out_bad null \
-derep 14 \
-trim_to_len 63 \
-trim_left 0 \
-trim_right 0 \
-min_len 62 \
-out_good $sample
```

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) 0.39 was used to remove `Nextera Transposase` adapters from sample 24 as follows:

```bash
conda create -n trimmomatic -c bioconda trimmomatic
```

```bash
conda activate trimmomatic
```

```bash
trimmomatic \
SE \
-threads 10 \
-phred33 \
sample24.with-adapters.fastq \
sample24.trimmed.fastq \
ILLUMINACLIP:~/anaconda3/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36
```

```bash
conda deactivate
```

### FastQC & MultiQC

We assessed the quality of the reads with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.8 and [MultiQC](https://multiqc.info/) 1.11. We aimed for a minimal average Phred score of 30 in each postion.

## Run table file

Finally, you should have a `RunTable.csv` file (assumed to be located at `Data/Samples`) which resmbels `SraRunTable.txt` files and consists of the following columns:
1. `Run`
2. `AvgSpotLen`
3. `Tissue`
4. `Group`
5. Additional columns

For example:

| Run      | AvgSpotLen | Group                 | Age(weeks) | Genotype | Diet   | Sex    | Tissue |
| -------- | ---------- | --------------------- | ---------- | -------- | ------ | ------ | ------ |
| sample4  | 63         | 15-WT-Full-Male       | 15         | WT       | Full   | Male   | Liver  |
| sample35 | 63         | 6.5-Het-Fasted-Female | 6.5        | Het      | Fasted | Female | Liver  |

The columns' order doesn't matter.


# Hyper Editing & Cluster Screening

## Hyper Editing

ADAR’s activity includes dense clusters of editing sites, which results in reads that are difficult to align by regular alignment procedures. However, using a method previously published by Porath, 2014 ([see original paper](https://www.nature.com/articles/ncomms5726)), we could correctly align those hyper-edited reads to the genome. 
The outputs of this method are lists of hyper-edited sites and the corresponding ultra-edited (UE) regions they span.
In terms of reproducibility, this is the most challening part of our reserach. While [an old version exits on github](https://github.com/hagitpt/Hyper-editing), it is hard to use. Furthemore, it doesn't deal with stranded data, and thus can only consider 6 `ReferenceBase-to-MutatedBase` mismatches (e.g., `A2G`) instead of all possible 12.
For this research, we supply the needed output in `Data/HyperEditing`.

## Cluster Screening

To get a clear Signal-to-Noise ratio (SNR) as possible, we applied an approach similar to the one described by [Buchumenski, 2021](https://academic.oup.com/nar/article/49/8/4325/6238403), termed here as Cluster Screening.
Using a `distance` ∈ {original=20, 50, 100, 150, …, 500}, we made two demands for each `distance`. First, two different editing sites (e.g., `A2G` and `C2T`) cannot reside within `distance` bp next to each other. Second, each editing site must have a neighbor editing site of the same type, located no further than `distance` bp from him.
We dismissed editing sites failing to satisfy any of these.
In our case, we used UE regions defined by `distance` = 20 as they yielded both the best signal and SNR.

```bash
mkdir Data/HyperEditing/ClusterScreening.SE_0.05_0.6_30_0.6_0.1_0.8_0.2
```

```bash
conda activate Killifish-APRT-RNAEditing
```

```bash
python ClusterScreening/cluster_screening.py \
--in_dir Data/HyperEditing/UEdetect.SE_0.05_0.6_30_0.6_0.1_0.8_0.2/all.ES.bed_files \
--analyse_a2g_table Data/HyperEditing/statistic/analyse.stranded_SE_0.05_0.6_30_0.6_0.1_0.8_0.2.A2G \
--data_table Data/Samples/RunTable.csv \
--out_dir Data/HyperEditing/ClusterScreening.SE_0.05_0.6_30_0.6_0.1_0.8_0.2 \
--genome_file Data/Annotations/genome.fasta \
--repeats_file Data/Annotations/repeats.gff \
--cds_file Data/Annotations/CDS.gff \
--beds_prefix ESuniqS.stranded \
--stranded
```

A report is given in the `ClusterScreeningSummary.ipynb` notebook and its corresponding html version.

# Diff Exp

First, create a direcotry for all Salmon-related outputs.

```bash
mkdir Data/DiffExp
```

## Creating an index

You should have an [argparse.bash](https://github.com/nhoffman/argparse-bash) executable installed at `/private/common/Software/BashLibs/argparse-bash`. If located elsewhere (which is most likely), please go ahead and modify line 5 in the following `create_index_with_decoy.sh` script accordingly.

```
mkdir Data/DiffExp/SalmonIndex
```

```
DiffExp/create_index_with_decoy.sh \
--genome Data/Annotations/genome.fasta \
--transcriptome Data/Annotations/transcriptome.fasta \
--uncompress_fasta_genome cat \
--uncompress_fasta_transcriptome cat \
--name_genome_index N.furzeri \
--salmon_index_dir Data/DiffExp/SalmonIndex
```

We thank [Roni Fultheim](https://github.com/roni-fultheim) for this script.

## Creating a tab-delimited transcript-to-gene map

```
python DiffExp/create_gene_map.py \
--transcriptome Data/Annotations/transcriptome.fasta \
--gff Data/Annotations/genomic.gff \
--id_prefix "rna-" \
--parent_prefix "gene-" \
--out_path Data/DiffExp/GeneMap.tab
```
  
## Running salmon

Finally, run Salmon itself.

```
mkdir Data/DiffExp/Salmon
```

```
python DiffExp/run_salmon.py \
--gene_map Data/DiffExp/GeneMap.tab \
--salmon_index Data/DiffExp/SalmonIndex/N.furzeri \
--in_dir Data/Samples \
--out_dir Data/DiffExp/Salmon
```

# Alignment

```
mkdir Data/BWA
```

## Create index

```
mkdir Data/BWA/Index
```

```
cd Data/BWA/Index
```

```
nohup bwa index \
-p N.furzeri \
../../Annotations/genome.fasta \
> index.out &
```

## Align

```
cd -
```

```
mkdir Data/BWA/Alignment
```

```
python Alignment/run_bwa.py \
--in_dir Data/Samples \
--postfix .fastq.lzma \
--decompress_cmd lzcat \
--out_dir Data/BWA/Alignment \
--index_dir Data/BWA/Index \
--index_name N.furzeri \
--bwa_algorithm mem \
--parity se
```

# RNA Editing Index

Now we quantify all A-to-I editing activity withing the hyper-edited regions we found earlier.
You should have the [RNAEditingIndexer](https://github.com/a2iEditing/RNAEditingIndexer) locally installed so that you could supply its path to the following `he_index.py` script. If you prefer to dockerized version, you'll have to modify that `Index/he_index.py` script, and also the `index.py` script, accordingly (both scripts are found at the `Index` directory).

```
mkdir Data/Index
```

```
python Index/he_index.py \
--he_cluster_dir Data/HyperEditing/ClusterScreening.SE_0.05_0.6_30_0.6_0.1_0.8_0.2/Original \
--root_dir Data/BWA/Alignment \
--stranded_he \
--data_table Data/Samples/RunTable.csv \
--salmon_dir Data/DiffExp/Salmon \
--genome_file Data/Annotations/genome.fasta \
--gff Data/Annotations/genomic.gff \
--id_prefix "gene-" \
--out_dir Data/Index \
--index_path RNAEditingIndex1.1
```

Reports are given in the `HEIndexAnalysis.ipynb` and `PerRegionPerSamplePCA.ipynb` notebooks and their corresponding html versions. 
They include all sorts of visualziations and statistical tests, so if you have a question regarding any of these, you are very wellcome to ask.