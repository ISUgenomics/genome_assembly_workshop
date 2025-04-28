# run 3ddna on arabidopsis genome

3d-dna is taking your HiC mapping information and converting it into connections between contigs, with contigs that are closer in the genome biologically having more HiC signal between them than different chromosome entities. 


### Set up the environment
```bash
#starting directory
/90daydata/shared/rick.masonbrink

#clone the 3ddna current repository
git clone https://github.com/aidenlab/3d-dna.git

# navigate to cloned folder
cd 3d-dna

# softlink the files from the juicer run
ln -s ../aligned/merged_nodups.txt
ln -s ../references/Genome.fasta
```


### Run 3ddna

There are many ways to optimize 3ddna. Those options can be seen with sh run-asm-pipeline.sh --help . In my experience I have found the manual scaffolding with Juicebox to be much easier when I use the initial scaffolding that 3ddna produces: Genome.0.hic and Genome.0.assembly. So I typically use the default parameters.
```bash
ml dna_3d; bash run-asm-pipeline.sh -h 


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
```

### Different ways to run 3d-dna

**Default run of 3ddna**
```bash
module load dna_3d;module load parallel;module java;bash run-asm-pipeline.sh Genome.fasta merged_nodups.txt 
```
**3ddna run that will scaffold all contigs >1000bp** <br>
Default skips anything smaller than has anything smaller than 15,000 skipped.
```bash
module load dna_3d;module load parallel;module load java;bash run-asm-pipeline.sh -i 1000 Genome.fasta merged_nodups.txt 
```
**3ddna run that allows a large amount of variation in repeat coverage** <br>
The fewer splits in scaffolds in the final output <br>
```bash
module load dna_3d;module load parallel; module load java;bash run-asm-pipeline.sh Genome.fasta merged_nodups.txt --editor-repeat-coverage 20
```
**3ddna run that allows a larger variation of repeat coverage**
```bash
module load dna_3d;module load parallel; module load java;bash run-asm-pipeline.sh --editor-repeat-coverage 20 Genome.fasta merged_nodups.txt 
```
**3Ddna run that trys to assemble unmapped reads with Lastz, if installed**
```bash
module load dna_3d;module load parallel;module load java; bash run-asm-pipeline.sh -m diploid Genome.fasta merged_nodups.txt
```

### Juicebox evaluation 

**Installation (none, just a download)**
```
Download the juicebox software that is appropriate for your personal computer. https://github.com/aidenlab/Juicebox/wiki/Download
There is an online version, but if your genome is large or you have lots of reads, it will be very slow to load between zooms.  
```

**Loading your genome into juicebox**
```
File, Open, select your HiC file.  Click the plot, otherwise the next step will not work. 
Assembly, Import Map Assembly, Click your .assembly file. This may take a few minutes to load. 
```

**Modification of your assembly in Juicebox**
* Selecting contigs
* Flipping contigs
* Moving contigs
* Scaffolds/Chromosomes
* Telomeres and Centromeres

**Juicebox normalization**
* None -- raw HiC contacts, can be noisy
* Coverage --  Normalizes each row and column by the total number of contacts -- fast but less precise for editing
* Coverage Square Root -- Similar to coverage, but normalizes by square root of the row and column sums
* Balanced -- Full matrix balancing, adjusts matices so each row and column sums to the same total. Best for accurate contig ordering and scaffolding. 


### How to see and fix misjoins and inversions
https://gitlab.com/wtsi-grit/rapid-curation/-/blob/main/Interpreting_HiC_Maps_guide.pdf <br>
<br>
Tip:  If you shift-select a contig and then shift-select another contig and you will lag. Always deselect the contig first by selecting outside of the box.
<br>

**Saving your work (do it often!)**
Juicebox can randomly have bugs and close. Make sure that you are regularly exporting your assembly file.   <br>
Assembly, Export Assembly  -->  the program will automatically rename your .assembly file to .assembly.review. Name it how you like. <br>

**Reloading your unfinished assembly into juicebox for further correction.** 
One extra step is required. After loading HiC, clicking the plot, and loading the original unmodified .assembly file, you need to load your modified assembly file after. Assembly, Load Modified Assembly, select your latest genome.assembly.review file. 

### 3ddna Finalization 

This is where 3ddna takes all of your modifications in juicebox and recreates the scaffolded genome. If I am sure that I am done with juicebox, I go ahead and sort the output so that the largest scaffolds will begin with the lowest numbers. This adds 100bp gaps between each contig you scaffolded. The '-s' is to set the stage of assembly to start at. We can put  indicating that we have already done the juiceboxing and want to try to have the cpu see if it can place any of the remaining smaller contigs. 

```bash
module load dna_3d;module load java;module load parallel; bash run-asm-pipeline-post-review.sh --sort-output -s seal -i 100 -r Genome.0.review.assembly Genome.fasta merged_nodups.txt
```



