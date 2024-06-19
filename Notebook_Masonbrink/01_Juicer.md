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