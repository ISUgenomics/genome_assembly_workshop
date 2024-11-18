## Run `hifiasm` on the Hifi reads

```bash
module load hifiasm

mkdir 03_Assembly

time hifiasm -o 03_Assembly/chr2_hifi.asm -t 20 -m 10 04_Chr2Fastq/mapped_reads.chr2.filtlong.fastq.gz
```

Time:
<pre>
real    3m42.123s
</pre>

### Converting GFA to fa

```bash
awk '/^S/ { print ">"$2; print $3 }' 03_Assembly/chr2_hifi.asm.bp.p_ctg.gfa > 03_Assembly/chr2_hifi.asm.bp.p_ctg.fa
```

The number of contigs: 
```bash
grep ">" -c 03_Assembly/chr2_hifi.asm.bp.p_ctg.fa
```
<pre>
63
</pre>

### Busco analysis with `compleasm`:

BUSCO (Benchmarking Universal Single-Copy Orthologs) is used to assess the completeness and quality of genome assemblies and annotations by comparing them to a database of conserved single-copy orthologs.

Downloading the Arabidopsis database on ceres.

workdir:`/project/isu_gif_vrsc/satheesh/07_AIUserForum_GenomeAssembly_Workshop_Nov2024/busco`

Going to use the `embryophyta_odb10` database for *Arabidopsis thaliana*.

```bash
# requested a node before executing the following commands.
module load busco5
busco --download embryophyta_odb10
```
The data base was transferred to Atlas. The next set of steps will be executed on Atlas.

workdir:`/project/isu_gif_vrsc/satheesh/07_AIUserForum_GenomeAssembly_Workshop_Nov2024/busco/08_Compleasm`

```bash
# 
software=/project/gif_vrsc_workshop/software/

time $software/compleasm_kit/compleasm.py run -t 20 \
  -l eukaryota -L 08_Compleasm/busco/busco_downloads/lineages/embryophyta_odb10/ \
  -a 05_FilteredReadsAssembly/chr2.filtlong_p2000.asm.bp.p_ctg.fa -o 08_Compleasm 
```
<pre>
real    1m39.890s
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

### Compleasm: Full genome assembly

```bash
time $software/compleasm_kit/compleasm.py run -t 20 \
  -l eukaryota -L 08_Compleasm/busco/busco_downloads/lineages/embryophyta_odb10/ \
  -a references/Genome.FINAL.fasta -o 08_Compleasm
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
  k=17 count output AT_HiFi.meryl 04_Chr2Fastq/mapped_reads.chr2.filtlong.fastq.gz
# real    0m30.826s

mkdir 07_Merqury_Output_HiFi && cd 07_Merqury_Output_HiFi
$software/merqury/merqury.sh ../AT_HiFi.meryl \
  ../05_FilteredReadsAssembly/chr2.filtlong_p2000.asm.bp.p_ctg.fa merqury_out
```
<pre>
real    0m34.100s
></pre>

```bash
cat merqury_out.qv
```

<pre>
chr2.filtlong_p2000.asm.bp.p_ctg        183     28815334        64.2762 3.73577e-07
# assembly        132080  12207077        31.94   0.000639731
</pre>

```bash
cat merqury_out.completeness.stats
```

<pre>
chr2.filtlong_p2000.asm.bp.p_ctg        all     18214157        18235904        99.8807
# assembly        all     8305853 8543157 97.2223
</pre>