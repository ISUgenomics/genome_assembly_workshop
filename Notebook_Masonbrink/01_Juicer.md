# Run juicer on arabidopsis 

There are a few folders that you will need to get ready before submitting juicer
You will need a reference folder that has your genome.fasta and its BWA-mem index

### Setting up the references folder
```
/work/gif3/masonbrink/USDA/02_JuicerTutorial

mkdir references ; cd references
ln -s  /work/gif3/satheesh/2024_GenomeAssemblyWorkshop/08_Arabidopsis_HiFi_HiC_Illumina/data_to_share/AT.asm.hic.p_ctg.fasta Genome.fasta
bwa index Genome.fasta

It is best to avoid using any special characters in your fasta headers/names.  This can cause issues in juicer. If you have a larger genome >700mb that is heavily fragmented >2000 contigs, I would suggest removing anything shorter than 5kb from the assembly before running juicer. The reason is that Juicebox will be dependent on the power of your personal computer, and lots of small contigs will result in large lag times for loading different zoom levels of your genome. For smaller genomes with less than 500 contigs the placement of many little contigs will be less cumbersome. Another reason is that typically these little contigs do not have great resolution in juicebox and may result in erroneous scaffolding.
```

### Create a chromosome sizes file 
```
/work/gif3/masonbrink/USDA/02_JuicerTutorial
ml bioawk
bioawk -c fastx '{print $name"\t"length($seq)}' references/Genome.fasta >chrom.sizes
```

### Create the fastq folder
Be aware of how I named the fastq files.  Juicer requires that they be named with this extension "_R1.fastq" and "_R2.fastq". Any other naming scheme will fail and they must also be unzipped. 
```
/work/gif3/masonbrink/USDA/02_JuicerTutorial
mkdir fastq; cd fastq/
cp  /work/gif3/satheesh/2024_GenomeAssemblyWorkshop/08_Arabidopsis_HiFi_HiC_Illumina/data_to_share/AT_HiC_1.fastq.gz AtHic_R1.fastq.gz
cp  /work/gif3/satheesh/2024_GenomeAssemblyWorkshop/08_Arabidopsis_HiFi_HiC_Illumina/data_to_share/AT_HiC_2.fastq.gz AtHic_R2.fastq.gz

gunzip *
```

### Create the splits folder
Typically juicer will create this folder for you and split the fastq files. It doesnt always work automatically, and you can manipulate the number of jobs your juicer run will submit. What we are aiming for is to get the maximum number of jobs that you can get to finish within your selected queue's time. For example, getting all of your jobs to run in the 1hr queue, instead of 48 hrs in a long node.   
```
/work/gif3/masonbrink/USDA/02_JuicerTutorial/fastq

#compute how many reads we have
wc -l AtHic_R1.fastq
281915000 
281915000/4=70,478,750 reads

#genome size 
192,720,089bp

#coverage is reads multiplied by read length, then divided by genome size
(70,478,750 *300)/192,720,089 = 109.7x

The

mkdir splits; cd splits
split -a 3 -l 10000000 -d --additional-suffix=_R1.fastq ../fastq/AtHic_R1.fastq &
split -a 3 -l 10000000 -d --additional-suffix=_R2.fastq ../fastq/AtHic_R1.fastq &

```


### Run juicer
I have kept using 1.5.7, as our 1.6 module has problems
```
/work/gif3/masonbrink/USDA/02_JuicerTutorial

ml juicer/1.5.7;ml bwa; juicer.sh -d /work/gif3/masonbrink/USDA/02_JuicerTutorial -p chrom.sizes -s none -z references/Genome.fasta  -q nova -Q 2:00:00 -l nova -L 12:00:00 -t 8 

```


### Desired results
A merged_nodups.txt file is all that you'll need to run the next step of the pipeline
```
ls -lrth /work/gif3/masonbrink/USDA/02_JuicerTutorial/aligned

-rw-r--r--. 1 remkv6 its-hpc-nova-gif 695M Jun 19 17:05 merged_sort.txt
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 2.8M Jun 19 17:06 opt_dups.txt
-rw-r--r--. 1 remkv6 its-hpc-nova-gif  46M Jun 19 17:06 dups.txt
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 647M Jun 19 17:06 merged_nodups.txt
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 2.2G Jun 19 17:09 abnormal.sam
-rw-r--r--. 1 remkv6 its-hpc-nova-gif    0 Jun 19 17:09 unmapped.sam
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 1.9K Jun 19 17:10 inter.txt
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 8.0K Jun 19 17:10 inter_hists.m
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 1.9K Jun 19 17:10 inter_30.txt
-rw-r--r--. 1 remkv6 its-hpc-nova-gif 7.8K Jun 19 17:10 inter_30_hists.m
-rw-r--r--. 1 remkv6 its-hpc-nova-gif    1 Jun 19 17:10 collisions.txt
-rw-r--r--. 1 remkv6 its-hpc-nova-gif  19M Jun 19 17:12 inter_30.hic
-rw-r--r--. 1 remkv6 its-hpc-nova-gif  20M Jun 19 17:12 inter.hic
drwxr-sr-x. 2 remkv6 its-hpc-nova-gif    2 Jun 19 17:13 inter_30_contact_domains

```

### Troublshooting common problems

##### Deduplication is not finishing
```
If this occurs, it is likely due to a low complexity region that has too many reads and thus too much depth for deduplication and becomes a memory hog without any progress. 

You can create a blacklist for these regions in your genome. Typically you can find those regions by running a repeat finder on your genome, and then masking those regions so reads will not map there. Then before you run 3ddna you can swap your genome back to the unmasked version.  Tandem repeats and low complexity repeats are usually the culprit, ribosomal DNA may cause problems as well.  
```

##### My juicer script wont submit any jobs
```
Juicer has not been modified suitably to use your HPC system. If this is not the issue, then it is likely that you must make -q and -l match the names of your queue's.   
```

##### My alignments are not finishing
```
Create a larger number of split fastq files within splits/.
```



