Date: 15 April, 2024
workdir: /work/gif3/satheesh/2024_GenomeAssemblyWorkshop

The data used in this tutorial was published along with this paper: [Engineered yeast genomes accurately assembled from pure and mixed samples](https://www.nature.com/articles/s41467-021-21656-9)

Downloading the data set:

**Whole genome Illumina sequencing of Y. lipolytica PO1f**

Time: ~15 mins
```bash
mkdir 01_Data && cd 01_Data
fasterq-dump SRR13512347
# spots read      : 1,580,326
# reads read      : 1,580,326
# reads written   : 1,580,326
```
Downloading the nanopore data. This took approx. 15 mins.

[Downloading the illumina data](https://www.ncbi.nlm.nih.gov/sra/SRX9923607).

Design: Genomic DNA isolated using a modified version of Promega's Wizard Genomic DNA Purification Kit. Library created using Illumina's Nextera DNA Flex Library Prep Kit.

Time: ~30 seconds 
```bash
time fasterq-dump --split-files SRR13512348
# spots read      : 846,718
# reads read      : 1,693,436
# reads written   : 1,693,436
# 
# real    0m30.169s
# user    0m9.672s
# sys     0m2.128s
cd ..
```

### Nanopore reads QC
Note: `nanoQCtrim` is DSL1 syntax. So, use Nextflow 22.10.x or earlier.
```bash
ml nextflow/22.10.4-wv2pxwv
ml singularity
export NXF_SINGULARITY_CACHEDIR=/work/gif3/satheesh/NXFContainers
time nextflow run isugifNF/nanoQCtrim --fastqs 01_Data/SRR13512347.fastq -profile nova,singularity
```
Error:
```
Error executing process > 'runNanoPlot (1)'

Caused by:
  Process `runNanoPlot (1)` terminated with an error exit status (1)

Command executed:

  NanoPlot -t 16 --huge --verbose --store -o SRR13512347 -p SRR13512347_  -f png --loglength --dpi 300 --plots {'kde','hex','dot'} --ti
tle SRR13512347" Nanopore Sequence" --N50 --fastq_rich SRR13512347.fastq
  NanoPlot -t 16 --huge --verbose --store -o SRR13512347 -p SRR13512347_ -f pdf --loglength --plots {'kde','hex','dot'} --title SRR1351
2347" Nanopore Sequence" --N50 --pickle SRR13512347/SRR13512347_NanoPlot-data.pickle

  ## Run Markdown generator
  nanoPlotMDGenerator.sh SRR13512347

Command exit status:
  1

Command output:



  If you read this then NanoPlot 1.32.0 has crashed :-(
  Please try updating NanoPlot and see if that helps...
  ```
  
Installing [NanoPlot](https://github.com/wdecoster/NanoPlot): 
```bash
conda activate nanopore

pip install NanoPlot

NanoPlot -t 10 --huge --verbose --store -o NanoPlot_out -p NanoPlot_out_  -f png --loglength --dpi 300 --plots {'kde','hex','dot'} --title ${NanoPlot}" Nanopore Sequence" --N50 --fastq 01_Data/SRR13512347.fastq

real    3m29.457s
user    3m18.076s
sys     0m2.436s
```

In the `NanoPlot` command, `--fastq-rich` switch did not work. There was the error that was obtained with `nanoQCtrim`. So, using `--fastq`. 

Installing [downpore](https://github.com/jteutenberg/downpore):

```bash
mkdir 02_Programs && cd 02_Programs

wget https://github.com/jteutenberg/downpore/releases/download/0.3.4/downpore.gz

unpig downpore.gz

chmod 755 downpore

export PATH=`pwd`:$PATH

## **Download Adapters**
wget https://github.com/jteutenberg/downpore/raw/master/data/adapters_front.fasta

wget https://github.com/jteutenberg/downpore/raw/master/data/adapters_back.fasta

cd ..

time downpore trim -i 01_Data/SRR13512347.fastq -f 02_Programs/adapters_front.fasta -b 02_Programs/adapters_back.fasta --himem true --num_workers 10 > 01_Data/SRR13512347_trimmed.fastq

# real    0m51.758s
# user    1m36.859s
# sys     0m15.032s
```
Genome assembly with [flye]():

```bash
time /work/gif3/satheesh/programs/Flye/bin/flye \
  --nano-raw 01_Data/SRR13512347_trimmed.fastq \
  --out-dir 03_FlyeOut \
  --genome-size 20m \
  --threads 10 \
  -i 4
```
Error:
```
[2024-04-16 15:32:18] ERROR: The input contain reads with duplicated IDs. Make sure all reads have unique IDs and restart. The first pr
oblematic ID was: SRR13512347.33623
[2024-04-16 15:32:19] ERROR: Command '['flye-modules', 'assemble', '--reads', '/work/gif3/satheesh/2024_GenomeAssemblyWorkshop/01_Data/
SRR13512347_trimmed.fastq', '--out-asm', '/work/gif3/satheesh/2024_GenomeAssemblyWorkshop/03_FlyeOut/00-assembly/draft_assembly.fasta',
 '--config', '/work/gif3/satheesh/programs/Flye/flye/config/bin_cfg/asm_raw_reads.cfg', '--log', '/work/gif3/satheesh/2024_GenomeAssemb
lyWorkshop/03_FlyeOut/flye.log', '--threads', '10', '--genome-size', '20000000', '--min-ovlp', '1000']' returned non-zero exit status 1
.
[2024-04-16 15:32:19] ERROR: Pipeline aborted
```
To address this problem, we add an incremental number to the readnames. If the original header is: `@SRR13512347.33623 33623 length=3310_(left)` and `@SRR13512347.33623 33623 length=3310_(left)`, the modified header will look like: `@SRR13512347.33623_1579993 33623 length=3310_(left)` and `@SRR13512347.33623_1579996 33623 length=3310_(right)`.

```bash
awk '{if (NR % 4 == 1) printf("%s_%d%s\n", $1, ++i, substr($0, length($1) + 1)); else print $0;}' 01_Data/SRR13512347_trimmed.fastq > 01_Data/SRR13512347_trimmed_readname_edited.fastq
```
We run the `flye` assembly again with the modified data set. 

```bash
time /work/gif3/satheesh/programs/Flye/bin/flye \
  --nano-raw 01_Data/SRR13512347_trimmed_readname_edited.fastq \
  --out-dir 03_FlyeOut \
  --genome-size 20m \
  --threads 10 \
  -i 4

# Total length:   29598046
# Fragments:      46
# Fragments N50:  3322132
# Largest frg:    4185809
# Scaffolds:      0
# Mean coverage:  113
# 
# [2024-04-16 17:38:46] INFO: Final assembly: /work/gif3/satheesh/2024_GenomeAssemblyWorkshop/03_FlyeOut/assembly.fasta
# 
# real    76m25.927s
# user    662m32.529s
# sys     5m54.406s
```

Using 36 cpus: will it reduce the time taken to perform the assembly?
```bash
time /work/gif3/satheesh/programs/Flye/bin/flye \
  --nano-raw 01_Data/SRR13512347_trimmed_readname_edited.fastq \
  --out-dir 03b_FlyeOut_36cpus \
  --genome-size 20m \
  --threads 36 \
  -i 4


# Total length:   29506120
# Fragments:      38
# Fragments N50:  3290781
# Largest frg:    4185842
# Scaffolds:      0
# Mean coverage:  110
# 
# [2024-04-16 18:34:10] INFO: Final assembly: /work/gif3/satheesh/2024_GenomeAssemblyWorkshop/03_FlyeOut_36cpus/assembly.fasta
# 
# real    35m34.244s
# user    640m14.014s
# sys     6m32.266s
```
There is a reduction in the number of contigs and the total length of the assembly. Why is this the case?

BUSCO analysis:

```bash
conda activate busco

export AUGUSTUS_CONFIG_PATH=/work/gif3/satheesh/mambaforge/envs/busco/config/

DB=/work/gif3/satheesh/2024_Zengyi_Shao_YeastMutation/05_2016_2_sample11_FlyeAssembly/04_BUSCO/busco_downloads/lineages/fungi_odb10/

time busco -i 03_FlyeOut/assembly.fasta -l $DB -o 03_FlyeOut/busco -m genome --augustus -c 36 -r

# --------------------------------------------------
# |Results from dataset                             |
# --------------------------------------------------
# |C:92.5%[S:59.9%,D:32.6%],F:0.4%,M:7.1%,n:758     |
# |701    Complete BUSCOs (C)                       |
# |454    Complete and single-copy BUSCOs (S)       |
# |247    Complete and duplicated BUSCOs (D)        |
# |3      Fragmented BUSCOs (F)                     |
# |54     Missing BUSCOs (M)                        |
# |758    Total BUSCO groups searched               |
# --------------------------------------------------

# real    4m31.051s
# user    64m39.279s
# sys     17m57.883s
```
For the second assembly:

```bash
time busco -i 03b_FlyeOut_36cpus/assembly.fasta -l $DB -o 03b_FlyeOut_36cpus/busco -m genome --augustus -c 36 -r

# --------------------------------------------------
# |Results from dataset                             |
# --------------------------------------------------
# |C:92.8%[S:60.3%,D:32.5%],F:0.5%,M:6.7%,n:758     |
# |703    Complete BUSCOs (C)                       |
# |457    Complete and single-copy BUSCOs (S)       |
# |246    Complete and duplicated BUSCOs (D)        |
# |4      Fragmented BUSCOs (F)                     |
# |51     Missing BUSCOs (M)                        |
# |758    Total BUSCO groups searched               |
# --------------------------------------------------
# 
# real    4m16.042s
# user    65m2.869s
# sys     17m56.090s
```
This second assembly also showed a similar BUSCO completeness to the first. 

Going for a third assembly with `40X` coverage to see if there is an improvement in the assembly. 

```bash
time /work/gif3/satheesh/programs/Flye/bin/flye \
  --nano-raw 01_Data/SRR13512347_trimmed_readname_edited.fastq \
  --out-dir 03c_FlyeOut_40x \
  --genome-size 20m \
  --asm-coverage 40 \
  --threads 36 \
  -i 4

# Total length:   29612714
# Fragments:      18
# Fragments N50:  3322145
# Largest frg:    4214499
# Scaffolds:      0
# Mean coverage:  118

# real    27m53.323s
# user    431m20.782s
# sys     5m48.169s
```
This assembly is down to 18 contigs at 40x coverage and slightly longer than the other two assemblies. So, I am going to use the `--scaffold` parameter. 

```bash
conda activate busco

export AUGUSTUS_CONFIG_PATH=/work/gif3/satheesh/mambaforge/envs/busco/config/

DB=/work/gif3/satheesh/2024_Zengyi_Shao_YeastMutation/05_2016_2_sample11_FlyeAssembly/04_BUSCO/busco_downloads/lineages/fungi_odb10/

time busco -i 03c_FlyeOut_40x/assembly.fasta -l $DB -o 03c_FlyeOut_40x/busco -m genome --augustus -c 36 -r

#  --------------------------------------------------
#  |Results from dataset                             |
#  --------------------------------------------------
#  |C:92.9%[S:59.8%,D:33.1%],F:0.4%,M:6.7%,n:758     |
#  |704    Complete BUSCOs (C)                       |
#  |453    Complete and single-copy BUSCOs (S)       |
#  |251    Complete and duplicated BUSCOs (D)        |
#  |3      Fragmented BUSCOs (F)                     |
#  |51     Missing BUSCOs (M)                        |
#  |758    Total BUSCO groups searched               |
#  --------------------------------------------------
# 
# real    4m48.820s
# user    66m41.229s
# sys     30m32.154s
# ```
```
There is no significant improvement in the assembly completeness.

```bash
time /work/gif3/satheesh/programs/Flye/bin/flye \
  --nano-raw 01_Data/SRR13512347_trimmed_readname_edited.fastq \
  --out-dir 03d_FlyeOut_30x \
  --genome-size 20m \
  --asm-coverage 30 \
  --threads 36 \
  -i 4

# Total length:   29586930
# Fragments:      21
# Fragments N50:  3351456
# Largest frg:    4185720
# Scaffolds:      0
# Mean coverage:  118

# real    26m3.774s
# user    425m22.490s
# sys     5m33.477s
```
This assembly gives 21 contigs. I still need to run Busco to check its completeness. Before that let me run a `Shasta` assembly to see what I get. 

```bash
mkdir 04_Shasta && cd 04_Shasta
ln -s ../01_Data/SRR13512347_trimmed_readname_edited.fastq
time seqtk seq -a SRR13512347_trimmed_readname_edited.fastq > SRR13512347_trimmed_readname_edited.fasta

# real    0m6.861s
# user    0m0.730s
# sys     0m1.841s

time /work/gif3/satheesh/programs/shasta-Linux-0.10.0 --input SRR13512347_trimmed_readname_edited.fasta --config Nanopore-May2022
```

Going for an hybrid assembly with `masurca`. 

