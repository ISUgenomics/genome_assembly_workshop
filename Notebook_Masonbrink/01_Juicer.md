# Introduction to Hi-C Data and Genome Scaffolding with Juicer, 3D-DNA, and Juicebox

##### What is Hi-C?

* Hi-C is a chromosome conformation capture technique that measures the 3D spatial organization of genomes.
* It captures physical proximity between different regions of the genome by crosslinking, digesting, and ligating DNA, followed by paired-end sequencing.
* Hi-C reads represent pairs of DNA fragments that were close together in the nucleus, even if they are distant in the linear genome sequence.
* Hi-C data can reveal:
    * Chromosome territories
    * Topologically associating domains (TADs)
    * Fine-scale contacts useful for genome assembly

##### Why Use Hi-C for Genome Scaffolding?

* Proximity information from Hi-C contacts can be used to:
    * Order contigs along chromosomes
    * Orient contigs correctly
    * Detect misassemblies in initial genome drafts

* Hi-C offers long-range linking information (>1 Mb) that is difficult to obtain from short-read sequencing alone.

##### How Juicer and 3D-DNA Use Hi-C Data for Scaffolding

**Juicer**

* Juicer is a pipeline for processing raw Hi-C reads into a contact map.
    * Aligns Hi-C read pairs to the draft genome assembly (e.g., using BWA or other aligners).
    * Filters and deduplicates valid Hi-C contacts.
    * Generates .hic files — compressed, indexed contact maps.
    * Provides quality control statistics on the data (e.g., contact matrices at various resolutions).

**3D-DNA**

* 3D-DNA is a pipeline that uses Juicer outputs to scaffold and improve genome assemblies.
    * Takes the initial assembly and Hi-C contact map as input.
    * Identifies and breaks potential misassemblies based on inconsistent Hi-C signals.
    * Clusters contigs into chromosomes using contact frequencies.
    * Orders and orients contigs within each chromosome.
    * Outputs a new assembly with chromosome-length scaffolds.
* Optionally, the assembly can be manually curated using Juicebox Assembly Tools to correct or refine scaffolding.


# Test run Juicer on arabidopsis 

### Start your interactive compute node

```bash
srun -N 1 -p interactive --ntasks-per-node=36 -t 5:00:00 --pty bash
```


### Setting up the references folder
There are a few folders that you will need to get ready before submitting juicer. <br>
You will need a reference folder that has your genome.fasta and its BWA-mem index
```bash
#move to our workshop directory
cd /90daydata/shared/

#create your own folder according to your user name
mkdir $USER ; cd $USER/
mkdir 01_Juicer; cd 01_Juicer

mkdir references; cd references

cp /project/scinet_workshop2/Bioinformatics_series/wk2_workshop/day2/02_Files/Genome.fasta  .

ml bwa
bwa index Genome.fasta

cd ../
```
It is best to avoid using any special characters in your fasta headers/names. This can cause issues in juicer. If you have a larger genome >700mb that is heavily fragmented >2000 contigs, I would suggest removing anything shorter than 5kb from the assembly before running juicer. The reason is that Juicebox will be dependent on the power of your personal computer, and lots of small contigs will result in lag times for loading different zoom levels of your genome. For smaller genomes with less than 500 contigs the placement of many little contigs will be less cumbersome. Another reason is that typically these little contigs do not have great resolution in juicebox and may result in erroneous scaffolding.

### Create a chromosome sizes file 
```bash
ml samtools
samtools faidx references/Genome.fasta 
cut -f 1,2 references/Genome.fasta.fai >chrom.sizes
```

### Create the fastq folder
Be aware of how I named the fastq files.  Juicer requires that they be named with this extension "_R1.fastq" and "_R2.fastq". Any other naming scheme will fail and they must also be unzipped. 
```bash
mkdir fastq; cd fastq/
cp /project/scinet_workshop2/Bioinformatics_series/wk2_workshop/day2/02_Files/*fastq .

cd ..
```

### Create the splits folder
Typically juicer will create this folder for you and split the fastq files. It doesnt always work automatically, and you can manipulate the number of jobs your juicer run will submit. What we are aiming for is to get the maximum number of jobs that you can get to finish within your selected queue's time. For example, getting all of your jobs to run in the 1hr queue, instead of 48 hrs in a long node.   
```bash
/90daydata/shared/rick.masonbrink

#Compute how many reads we have
wc -l 2MAtHiCDedup_R1.fastq
2000000 2MAtHiCDedup_R1.fastq
2000000/4=500,000 reads

#Total length of reads
500,000 *2 reads *150bp length = 150,000,000 bp

#Compute the genome size
cp /project/scinet_workshop2/Bioinformatics_series/wk2_workshop/day2/02_Files/new_Assemblathon.pl .
./new_Assemblathon.pl references/Genome.fasta


#Coverage is reads multiplied by read length, then divided by genome size.  Note this is a minimal amount of reads that I have already filtered and deduplicated to save time for our workshop. Deduplication of reads is the most time consuming part of this step. 
150,000,000/192,720,089 = 0.7783x coverage

mkdir splits; cd splits
split -a 3 -l 240000 -d --additional-suffix=_R1.fastq ../fastq/2MAtHiCDedup_R1.fastq &
split -a 3 -l 240000 -d --additional-suffix=_R2.fastq ../fastq/2MAtHiCDedup_R2.fastq &
cd ..
```

### Run juicer
```bash
/90daydata/shared/rick.masonbrink

ml juicer;JUICER juicer.sh -d /90daydata/shared/rick.masonbrink/01_Juicer -p chrom.sizes -s none -z references/Genome.fasta -t 8 --assembly
```


### Desired results
A merged_nodups.txt file is all that you'll need to run the next step of the pipeline
```
ls -lrth 

total 321M
-rw-r-----. 1 rick.masonbrink proj-vrsc  317 Apr 28 12:40 header
-rw-r-----. 1 rick.masonbrink proj-vrsc  14M Apr 28 12:45 merged1.txt
-rw-r-----. 1 rick.masonbrink proj-vrsc  13M Apr 28 12:45 merged30.txt
-rw-r-----. 1 rick.masonbrink proj-vrsc 128M Apr 28 12:45 merged_dedup.bam
-rw-r-----. 1 rick.masonbrink proj-vrsc 1.9K Apr 28 12:45 inter.txt
-rw-r-----. 1 rick.masonbrink proj-vrsc 7.1K Apr 28 12:45 inter_hists.m
-rw-r-----. 1 rick.masonbrink proj-vrsc 1.9K Apr 28 12:45 inter_30.txt
-rw-r-----. 1 rick.masonbrink proj-vrsc 7.1K Apr 28 12:45 inter_30_hists.m
-rw-r-----. 1 rick.masonbrink proj-vrsc 167M Apr 28 12:45 merged_nodups.txt
```

**Juicer2 Output Files — Basic Definitions**

| File | Definition | Use |
|:-----|:-----------|:----|
| `header` | Text file with metadata about the Hi-C run (e.g., genome version, parameters). | Used internally by Juicer; small summary of the processing run. |
| `merged1.txt` | List of Hi-C read pairs at 1 bp resolution (raw, not normalized). | Can be used to build low-resolution contact maps. |
| `merged30.txt` | List of Hi-C read pairs at 30 bp resolution (binning shortens file size). | Used for making higher-level, binned contact maps more efficiently. |
| `merged_dedup.bam` | BAM file of aligned, deduplicated Hi-C reads (valid pairs only). | Useful for visualization in genome browsers like IGV, JBROWSE, etc. |
| `inter.txt` | Contact statistics between contigs/scaffolds (raw form). | Basic QC: shows how reads link different parts of the genome. |
| `inter_hists.m` | MATLAB script with histograms of Hi-C contact distributions. | Helps visualize contact decay with distance; used for QC. |
| `inter_30.txt` | Same as `inter.txt`, but based on reads binned at 30 bp resolution. | Faster/lighter version for QC and distance plotting. |
| `inter_30_hists.m` | MATLAB script plotting histograms for 30 bp binned contacts. | Visual QC of contact maps at lower resolution. |
| `merged_nodups.txt` | Full list of valid Hi-C contacts, deduplicated, in tab-separated text format. | **Main input** for 3D-DNA scaffolding and for building `.hic` files for Juicebox visualization. |

### Troubleshooting common problems

##### Deduplication is not finishing
```
If this occurs, it is likely due to a low complexity region that has too many reads and thus too much depth for deduplication and becomes a memory hog without any progress. 

You can create a blacklist for these regions in your genome. Typically you can find those regions by running a repeat finder on your genome, and then masking those regions so reads will not map there. Then before you run 3ddna you can swap your genome back to the unmasked version.  Tandem repeats and low complexity repeats are usually the culprit, ribosomal DNA may cause problems as well.  
```

##### My juicer script wont submit any jobs
```
Juicer has not been modified suitably to use your HPC system. If this is not the issue, then it is likely that you must make -q and -l match the names of your queue's. We are using the CPU version of Juicer today, so we should not see any of these issues.
```

##### My alignments are not finishing
```
Create a larger number of split fastq files within splits/.
```



