---
title: "Assembly and Assessment"
authors:
  - "Viswanathan Satheesh"
  - "Sivanandan Chudalayandi"
date: "2024-11-20"
---

## Genome Assembly and Quality Assessment

This tutorial will guide you through the process of assembling a genome using HiFi reads and assessing its quality. We'll be working with *Arabidopsis thaliana* chromosome 2 data as an example.

## 1. Genome Assembly with `hifiasm`

`hifiasm` is a fast and accurate haplotype-resolved assembler designed specifically for PacBio HiFi reads. It's particularly effective because it can handle the high accuracy of HiFi reads while maintaining computational efficiency.

First, let's create a directory for our assembly:

```bash
mkdir 06_Assembly

# Run hifiasm with the following parameters:
# -o: output prefix
# -t: number of threads (20 in this case)
# -m: minimum number of overlaps to keep a contig (10 in this case)
time hifiasm -o 06_Assembly/chr2_hifi.asm -t 20 -m 10 01_Data/AT_HiFi_chr2.fastq.gz
```

Time:
<pre>
real    3m48.387s
</pre>

The assembly took approximately 3 minutes and 48 seconds to complete.

### Converting Assembly Format from GFA to FASTA

The assembly output is in GFA (Graphical Fragment Assembly) format, which we need to convert to FASTA format for downstream analysis. GFA is a format that represents genome graphs, while FASTA is a more widely used format for sequence data.

```bash
# Extract sequences from GFA and convert to FASTA format
awk '/^S/ { print ">"$2; print $3 }' 06_Assembly/chr2_hifi.asm.bp.p_ctg.gfa > 06_Assembly/chr2_hifi.asm.bp.p_ctg.fa
```

Let's check how many contigs we got from our assembly:
```bash
grep ">" -c 06_Assembly/chr2_hifi.asm.bp.p_ctg.fa
```
Our assembly resulted in 63 contigs.

<pre>
63
</pre>

## 2. Quality Assessment

### BUSCO Analysis with `compleasm`

BUSCO (Benchmarking Universal Single-Copy Orthologs) is a crucial tool for assessing genome assembly completeness. It searches for highly conserved genes that should be present in all species of a given lineage. The presence, absence, or fragmentation of these genes gives us insight into the quality of our assembly.

Key BUSCO metrics:
- **Complete (S)**: Single-copy complete genes
- **Duplicated (D)**: Complete genes present multiple times
- **Fragmented (F)**: Partially assembled genes
- **Missing (M)**: Genes not found in the assembly
- **Total (N)**: Total number of genes in the BUSCO dataset

We'll use `compleasm`, a faster implementation of BUSCO, with the `embryophyta_odb10` database which is specific for plants:

```bash
software=/project/gif_vrsc_workshop/software/

time $software/compleasm_kit/compleasm.py run -t 20 \
  -l eukaryota -L 01_Data/busco_downloads/lineages/embryophyta_odb10/ \
  -a 06_Assembly/chr2_hifi.asm.bp.p_ctg.fa -o 07_Compleasm 
```

Time:
<pre>
real    4m55.897s
</pre>

#### Understanding the BUSCO Results

Our assembly of chromosome 2 shows:
- **21.96% (56)** complete single-copy genes
- **0.00% (0)** duplicated genes
- **1.57% (4)** fragmented genes
- **76.47% (195)** missing genes
- **Total: 255** genes assessed

These numbers are expected to be low since we're only looking at one chromosome rather than the complete genome. For comparison, let's look at the metrics for the complete Arabidopsis genome:

```bash
# Note: Do not run this command - shown for comparison only
time $software/compleasm_kit/compleasm.py run -t 20 \
  -l eukaryota -L 01_Data/busco_downloads/lineages/embryophyta_odb10/ \
  -a references/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -o 07_Compleasm/
```

Time:
<pre>
real    13m26.798s
</pre>

The complete genome shows:
- **90.59% (231)** complete single-copy genes
- **9.41% (24)** duplicated genes
- **0.00% (0)** fragmented genes
- **0.00% (0)** missing genes
- **Total: 255** genes assessed

This comparison helps us understand that our chromosome 2 assembly, while not complete, is reasonable given we're only looking at one chromosome.

## 3. K-mer Based Quality Assessment

K-mer analysis is another powerful method for evaluating genome assembly quality. A k-mer is a substring of length k from a DNA sequence. By analyzing the distribution and frequency of k-mers in both our raw reads and assembly, we can assess assembly completeness and accuracy.

### Step 1: Determining Optimal K-mer Size

First, we need to determine the optimal k-mer length for our analysis. This is important because:
- Too small k-mers may lead to random matches
- Too large k-mers may miss legitimate matches
- The optimal size depends on genome size and error rates

```bash
software=/project/gif_vrsc_workshop/software

# Calculate optimal k-mer size based on genome size
$software/merqury/best_k.sh 19000000
```

Output:
<pre>
genome: 19000000
tolerable collision rate: 0.001
17.0719
</pre>

The script suggests a k-mer size of 17 based on our genome size of 19Mb and a tolerable collision rate of 0.001.

### Step 2: K-mer Analysis with Merqury

Merqury is a tool that evaluates genome assemblies using k-mer frequencies. It helps us understand:
- Assembly completeness
- Base-level accuracy
- Potential misassemblies

Let's run the analysis:

```bash
# Count k-mers in the input reads
time $software/meryl-1.4.1/bin/meryl \
  k=17 count output 08_AT_HiFi.meryl 01_Data/AT_HiFi_chr2.fastq.gz

# Create output directory and run Merqury
mkdir 09_Merqury_Output_HiFi && cd 09_Merqury_Output_HiFi
merqury.sh ../08_AT_HiFi.meryl \
  ../06_Assembly/chr2_hifi.asm.bp.p_ctg.fa merqury_out
```

Time taken:
<pre>
real    0m22.237s
</pre>

### Step 3: Interpreting Merqury Results

#### Quality Value (QV) Statistics
Let's examine the assembly quality values:

```bash
cat merqury_out.qv
```

<pre>
chr2_hifi.asm.bp.p_ctg  40      28858272        70.8866 8.15344e-08
</pre>

These metrics tell us:
- **Assembly Name**: chr2_hifi.asm.bp.p_ctg
- **Unique K-mers**: 40 k-mers found only in the assembly
- **Total K-mers**: 28,858,272 k-mers found in both assembly and reads
- **Quality Value (QV)**: 70.89 - this is excellent! QV is a logarithmic measure of error probability
- **Error Rate**: 8.15e-08 - extremely low, indicating very high accuracy

For comparison, a typical Oxford Nanopore assembly might show much lower values:
```
# Example: yeast assembly with nanopore reads
# assembly        132080  12207077        31.94   0.000639731
```

#### Completeness Assessment
Let's check the assembly completeness:

```bash
cat merqury_out.completeness.stats
```

<pre>
chr2_hifi.asm.bp.p_ctg  all     18214072        18235904        99.8803
</pre>

These statistics show:
- **Covered Bases**: 18,214,072 bases are supported by k-mers from our input reads
- **Total Bases**: 18,235,904 bases in our assembly
- **Completeness**: 99.88% - this is exceptional! It means almost every base in our assembly is supported by the input reads

For perspective, a typical Oxford Nanopore assembly might show lower completeness:
```
# Example: yeast assembly with nanopore reads
# assembly        all     8305853 8543157 97.2223
```

### Summary of Quality Assessment

Our assembly shows excellent quality metrics:
1. Very high completeness (99.88%)
2. Exceptional base accuracy (QV = 70.89)
3. Extremely low error rate (8.15e-08)

These metrics indicate that our HiFi-based assembly of chromosome 2 is highly accurate and nearly complete, which is exactly what we'd expect from HiFi data.
