Date: 28 May, 2024
workdir: /project/isu_gif_vrsc/satheesh/03_20240426_GenomeAssemblyWorkshop/2024_GenomeAssemblyWorkshop_S288C/07_HybridAssembly

`Filtlong` on Hifi reads to get 50x data. 
```bash
ln -s ../01_Data/SRR18210286_pass.fastq
ln -s ../01_Data/SRR17374239_*fastq.gz .
ln -s ../05_NanoporeReadAssembly/SRR17374240_800Mb.fastq.gz
```
Going to try to run `verkko` on Hifi and nanopore reads. 
```bash
conda activate /project/isu_gif_vrsc/satheesh/conda_envs
mamba install -c conda-forge -c bioconda -c defaults verkko

time verkko -d work --hifi SRR18210286_pass.fastq --nano SRR17374240_800Mb.fastq.gz
# real    221m54.331s
conda deactivate

# Installing compleasm
mamba create -n compleasm -c conda-forge -c bioconda compleasm
conda activate compleasm
cd work/
time compleasm run -a assembly.fasta -o compleasm_out -t 36 --autolineage
# 17 min

# S:49.46%, 1057
# D:49.93%, 1067
# F:0.05%, 1
# I:0.00%, 0
# M:0.56%, 12
# N:2137

# Run with -l saccharomycetes
time compleasm run -a assembly.fasta -o compleasm_out -t 36 -l saccharomycetes
cat work/compleasm_out/summary.txt
## lineage: saccharomycetes_odb10
# S:49.46%, 1057
# D:49.93%, 1067
# F:0.05%, 1
# I:0.00%, 0
# M:0.56%, 12
# N:2137
```

gfastats: <https://github.com/vgl-hub/gfastats>
```bash
mkdir 01_gfastats && cd 01_gfastats
ln -s ../work/assembly.fasta .
/project/isu_gif_vrsc/satheesh/programs/gfastats/build/bin/gfastats assembly.fasta > assembly.stats
cd ..
```

`purge_haplotigs`:
```bash
conda activate /project/isu_gif_vrsc/satheesh/conda_envs
mamba install -c bioconda purge_haplotigs
conda deactivate
mkdir 02_Purge_Haplotigs && cd 02_Purge_Haplotigs
# salloc -N1 -n20 -t 4:00:00
ml minimap2
# preparation
time minimap2 -t 20 -ax map-pb ../work/assembly.fasta ../SRR18210286_pass.fastq --secondary=no |samtools sort -m 1G -o aligned.bam -T tmp.ali
# 40m38.659s

# Step 1
conda activate /project/isu_gif_vrsc/satheesh/conda_envs
time purge_haplotigs hist -b aligned.bam -g ../work/assembly.fasta -t 20
# 0m41.795s

# Manual step see histogram generated
purge_haplotigs cov -i aligned.bam.200.gencov -l 5 -m 50 -h 150 -o coverage_stats.csv -j 80 -s 80
# Analysis finished successfully! Contig coverage stats saved to 'coverage_stats.csv'.
```
To set the parameters for purge_haplotigs, we need to identify suitable cutoffs from the histogram. Here’s how you can estimate the values:

-l / -low: This is the read depth low cutoff. Typically, this value is set to include the low-coverage peak (repetitive regions). From the histogram, this looks to be around 5.
-m / -mid: This is the low point between the haploid and diploid peaks. In the absence of a clear diploid peak, this can be approximated from the histogram. A reasonable value might be around 50.
-h / -high: This is the read depth high cutoff. This should be set to a value that includes the diploid peak and filters out high-coverage regions, which often represent duplications or repetitive elements. This can be set to around 150 based on the tail end of the histogram.

Need to add the histogram.

Step 3: 
```bash
time purge_haplotigs  purge  -g ../work/assembly.fasta  -c coverage_stats.csv
# PURGE HAPLOTIGS HAS COMPLETED SUCCESSFULLY!
# real    0m39.447s
cd ..
```
Need to go through the results.

Redundans: 
```bash
mkdir 03_Redundans && cd 03_Redundans

# Running in slurm script
/project/isu_gif_vrsc/satheesh/programs/redundans/redundans.py -v -i ../SRR17374239_?.fastq.gz -l ../SRR18210286_pass.fastq ../SRR17374240_800Mb.fastq.gz -f ../work/assembly.fasta -o run_short_long_populatescaffold --minimap2scaffold --populateScaffolds

time compleasm run -a run_short_long_populatescaffold/scaffolds.reduced.fa -o compleasm_out -t 20 -l saccharomycetes
# real    5m30.545s

# S:78.52%, 1678
# D:20.87%, 446
# F:0.05%, 1
# I:0.00%, 0
# M:0.56%, 12
# N:2137
```