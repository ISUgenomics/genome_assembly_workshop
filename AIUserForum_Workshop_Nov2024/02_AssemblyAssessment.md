---
title: "Assembly and Assessment"
authors:
  - "Viswanathan Satheesh"
  - "Sivanandan Chudalayandi"
date: "2024-11-20"
---

## Run `hifiasm` on the Hifi reads

```bash
mkdir 06_Assembly

time hifiasm -o 06_Assembly/chr2_hifi.asm -t 20 -m 10 01_Data/AT_HiFi_chr2.fastq.gz
```

Time:
<pre>
real    3m48.387s
</pre>

### Converting GFA to fa

```bash
awk '/^S/ { print ">"$2; print $3 }' 06_Assembly/chr2_hifi.asm.bp.p_ctg.gfa > 06_Assembly/chr2_hifi.asm.bp.p_ctg.fa
```

The number of contigs: 
```bash
grep ">" -c 06_Assembly/chr2_hifi.asm.bp.p_ctg.fa
```
<pre>
63
</pre>



### Busco analysis with `compleasm`:

BUSCO (Benchmarking Universal Single-Copy Orthologs) is a bioinformatics tool used to assess the completeness of genome assemblies, gene sets, or transcriptomes. It does this by comparing the sequences against a database of highly conserved, single-copy orthologs expected to be present in most species within a particular lineage.

Going to use the `embryophyta_odb10` database for *Arabidopsis thaliana*.

```bash
software=/project/gif_vrsc_workshop/software/
```
Busco not installed on Atlas.

```bash
time $software/compleasm_kit/compleasm.py run -t 20 \
  -l eukaryota -L 01_Data/busco_downloads/lineages/embryophyta_odb10/ \
  -a 06_Assembly/chr2_hifi.asm.bp.p_ctg.fa -o 07_Compleasm 
```
<pre>
real    4m55.897s
</pre>

#### Compleasm output

<pre>
 # S:21.96%, 56
 # D:0.00%, 0
 # F:1.57%, 4
 # I:0.00%, 0
 # M:76.47%, 195
 # N:255
</pre>

## Do not run
### Compleasm : Full genome assembly

```bash
time $software/compleasm_kit/compleasm.py run -t 20 \
  -l eukaryota -L 01_Data/busco_downloads/lineages/embryophyta_odb10/ \
  -a references/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -o 07_Compleasm/
```

<pre>
real    13m26.798s
</pre>

#### Compleasm output

<pre>
 S:90.59%, 231
D:9.41%, 24
F:0.00%, 0
I:0.00%, 0
M:0.00%, 0
N:255
</pre>

### Genome assembly evaluation with k-mers

Determining the optimal k-mer length (k) for analyzing a given set of sequencing reads
```bash
software=/project/gif_vrsc_workshop/software

$software/merqury/best_k.sh 19000000
```

<pre>
genome: 19000000
tolerable collision rate: 0.001
17.0719
</pre>

Evaluation with `merqury`.
```bash
time $software/meryl-1.4.1/bin/meryl \
  k=17 count output 08_AT_HiFi.meryl 01_Data/AT_HiFi_chr2.fastq.gz
# real    0m22.237s

mkdir 09_Merqury_Output_HiFi && cd 09_Merqury_Output_HiFi
merqury.sh ../08_AT_HiFi.meryl \
  ../06_Assembly/chr2_hifi.asm.bp.p_ctg.fa merqury_out
```
<pre>
real    0m22.237s
></pre>

```bash
cat merqury_out.qv
```

<pre>
chr2_hifi.asm.bp.p_ctg  40      28858272        70.8866 8.15344e-08

# Example: yeast assembly with nanopore reads
# assembly        132080  12207077        31.94   0.000639731
</pre>

The file `merqury_out.qv` contains quality value (QV) statistics for the genome assembly:

- **`chr2_hifi.asm.bp.p_ctg`**: The name of the assembly being evaluated.
- **`40`**: The number of k-mers uniquely found only in the assembly, indicating how much of the assembly is unique compared to the input read set.
- **`28858272`**: The total number of k-mers found in both the assembly and the read set.
- **`70.8866`**: The estimated quality value (QV) of the assembly, indicating the percentage of k-mers shared between the assembly and the reference read set.
- **`8.15344e-08`**: The estimated base error rate, representing the likelihood of an error occurring at a given base position.


```bash
cat merqury_out.completeness.stats
```

<pre>
chr2_hifi.asm.bp.p_ctg  all     18214072        18235904        99.8803

# Example: yeast assembly with nanopore reads
# assembly        all     8305853 8543157 97.2223
</pre>

The file `merqury_out.completeness.stats` contains statistics on the completeness of the genome assembly:

- **`chr2_hifi.asm.bp.p_ctg`**: The name of the assembly being evaluated.
- **`all`**: Refers to the entire genome or all contigs in the assembly.
- **`18214072`**: Number of bases that are present in both the k-mer database and the assembly (covered bases).
- **`18235904`**: Total number of bases in the assembly.
- **`99.8803`**: Percentage completeness, indicating that 99.88% of the assembly is covered by the k-mers, which implies a high-quality assembly.
