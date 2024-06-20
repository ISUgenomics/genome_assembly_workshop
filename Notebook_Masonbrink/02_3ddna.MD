# run 3ddna on arabidopsis genome

3d-dna is taking your HiC mapping information and converting it into connections between contigs, with contigs that are closer in the genome biologically having more HiC signal between them than different chromosome entities. 


### Set up the environment
```

#starting directory
/work/gif3/masonbrink/USDA/02_JuicerTutorial/3Ddna/

#clone the 3ddna current repository
git clone https://github.com/aidenlab/3d-dna.git

# navigate to cloned folder
cd 3d-dna

# softlink the files from the juicer run
ln -s ../../aligned/merged_nodups.txt
ln -s ../../references/Genome.fasta
```

### run 3ddna

There are many ways to optimize 3ddna.   Those options can be seen with sh run-asm-pipeline.sh --help . In my experience I have found the manual scaffolding to be much easier when I just use the initial scaffolding that 3ddna produces, the Genome.0.hic and Genome.0.assembly files that are produced. So I typically use the default parameters.
```
*****************************************************
3D de novo assembly: version 180114

USAGE: ./run-asm-pipeline.sh [options] <path_to_input_fasta> <path_to_input_mnd>

DESCRIPTION:
This is a script to assemble draft assemblies (represented in input by draft fasta and deduplicated list of alignments of Hi-C reads to this fasta as produced by the Juicer pipeline) into chromosome-length scaffolds. The script will produce an output fasta file, a Hi-C map of the final assembly, and a few supplementary annotation files to help review the result in Juicebox.

ARGUMENTS:
path_to_input_fasta                     Specify file path to draft assembly fasta file.
path_to_input_mnd                       Specify path to deduplicated list of alignments of Hi-C reads to the draft assembly fasta as produced by the Juicer pipeline: the merged_nodups file (mnd).

OPTIONS:
-m|--mode haploid/diploid                       Runs in specific mode, either haploid or diploid (default is haploid).
-i|--input input_size                   Specifies threshold input contig/scaffold size (default is 15000). Contigs/scaffolds smaller than input_size are going to be ignored.
-r|--rounds number_of_edit_rounds                       Specifies number of iterative rounds for misjoin correction (default is 2).
-s|--stage stage                                        Fast forward to later assembly steps, can be polish, split, seal, merge and finalize.
-h|--help                       Shows this help. Type --help for a full set of options.
*****************************************************


# default run of 3ddna
module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load openjdk/17.0.5_8-5pdrwz2;module load parallel/20220522-sxcww47;cd 3d-dna ;bash run-asm-pipeline.sh Genome.fasta merged_nodups.txt 

# run of 3ddna that will attempt to scaffold all contigs >1000bp. default has anything smaller than 15000 skipped.
module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load openjdk/17.0.5_8-5pdrwz2;module load parallel/20220522-sxcww47;cd 3d-dna ;bash run-asm-pipeline.sh Genome.fasta merged_nodups.txt -i 1000

# run of 3ddna that allows a large amount of variation in repeat coverage to occur, thus fewer splits in scaffolds in the final output
module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load openjdk/17.0.5_8-5pdrwz2;module load parallel/20220522-sxcww47;cd 3d-dna ;bash run-asm-pipeline.sh Genome.fasta merged_nodups.txt --editor-repeat-coverage 20

# run of 3ddna that allows a larger variation of repeat coverage and also trys to reassemble unmapped reads with Lastz, if installed
module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load openjdk/17.0.5_8-5pdrwz2;module load parallel/20220522-sxcww47;cd 3d-dna ;bash run-asm-pipeline.sh Genome.fasta merged_nodups.txt -m diploid --editor-repeat-coverage 20
```

### Juicebox evaluation 

Installation (none, just a download)
```
Download the juicebox software that is appropriate for your personal computer. https://github.com/aidenlab/Juicebox/wiki/Download
There is an online version, but if your genome is large or you have lots of reads, it will be very slow to load between zooms.  
```

Loading your genome into juicebox
```
File, Open, select your HiC file.  Click the plot, otherwise the next step will not work. 
Assembly, Import Map Assembly, Click your .assembly file. This may take a few minutes to load. 
```

Modification of contigs in juicebox
```
Selecting contigs
Flipping contigs
Moving contigs
Scaffolds/Chromosomes
Telomeres and Centromeres

shift select a contig and then shift select another contig and you will lag.
```

Saving your work (do it often!)
```
Juicebox can randomly have bugs and close.  Make sure that you are regularly exporting your assembly file.   
Assembly, Export Assembly  -->  the program will automatically rename your .assembly file to .assembly.review. Name it how you like. 
```

Reloading your unfinished genome. 
```
One extra step more is required.  After loading HiC, clicking the plot, and loading the original unmodified .assembly file, you need to load your modified assembly file after.   Assembly, Load Modified Assembly, select your latest genome.assembly.review file. 
```


### 3ddna Finalization 

This is where 3ddna takes all of your modifications in juicebox and recreates the scaffolded genome.  If I am sure that I am done with juicebox, I go ahead and sort the output so that the largest scaffolds will begin with the lowest numbers.   This adds 500bp gaps between each contig you scaffolded.   -s is indicating that we have already done the juicboxing and want to try to have the cpu see if it can place some any of the remaining smaller contigs. 

```
module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load jdk;module load parallel; bash run-asm-pipeline-post-review.sh --sort-output -s seal -i 500 -r Genome.0.review.assembly  Genome.fasta merged_nodups.txt
```



