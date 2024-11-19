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

BUSCO (Benchmarking Universal Single-Copy Orthologs) is used to assess the completeness and quality of genome assemblies and annotations by comparing them to a database of conserved single-copy orthologs.

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

## Below, compleasm run is not tested
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
  k=17 count output 08_AT_HiFi.meryl 01_Data/AT_HiFi_chr2.fastq.gz
# real    0m22.237s

mkdir 09_Merqury_Output_HiFi && cd 08_Merqury_Output_HiFi
$software/merqury/merqury.sh ../08_AT_HiFi.meryl \
  ../06_Assembly/chr2_hifi.asm.bp.p_ctg.fa merqury_out
```
<pre>
real    0m22.237s
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