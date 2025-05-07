# run 3D-dna on arabidopsis genome

3d-dna is taking your HiC mapping information and converting it into connections between contigs, with contigs that are closer in the genome biologically having more HiC signal between them than different chromosome entities. 


### Set up the environment
```bash
#starting directory
/90daydata/shared/rick.masonbrink

#clone the 3D-dna current repository. If 3D-dna is installed as a module on your system you may not need to clone the repo. 
git clone https://github.com/aidenlab/3d-dna.git

# navigate to cloned folder, or make named folder
cd 3d-dna

# softlink the files from the juicer run
ln -s ../aligned/merged_nodups.txt
ln -s ../references/Genome.fasta
```


### Run 3D-dna

There are many ways to optimize 3D-dna. Those options can be seen with sh run-asm-pipeline.sh --help . In my experience I have found the manual scaffolding with Juicebox to be much easier when I use the initial scaffolding that 3D-dna produces: Genome.0.hic and Genome.0.assembly. So I typically use the default parameters.
```bash

#for basic commands
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

#for advanced parameters
ml dna_3d; bash run-asm-pipeline.sh --help
```

### Different ways to run 3d-dna

**Set the temp directory for java before you run 3D-dna**
```bash
export _JAVA_OPTIONS=-Djava.io.tmpdir=${TMPDIR}
```
**Default run of 3D-dna**
```bash
module load dna_3d;module load parallel;module load java;bash run-asm-pipeline.sh Genome.fasta merged_nodups.txt 
```
**3D-dna run that will scaffold all contigs >1000bp** <br>
Default skips anything smaller than has anything smaller than 15,000 skipped.
```bash
module load dna_3d;module load parallel;module load java;bash run-asm-pipeline.sh -i 1000 Genome.fasta merged_nodups.txt 
```
**3D-dna run that allows a larger variation of repeat coverage**
```bash
module load dna_3d;module load parallel; module load java;bash run-asm-pipeline.sh --editor-repeat-coverage 20 Genome.fasta merged_nodups.txt 
```
**3D-dna run that trys to assemble unmapped reads with Lastz, if installed**
```bash
module load dna_3d;module load parallel;module load java; bash run-asm-pipeline.sh -m diploid Genome.fasta merged_nodups.txt
```
**3D-dna run that only uses reads with a mapping quality greater than 30**
```bash
module load dna_3d;module load parallel;module load java; bash run-asm-pipeline.sh -q 30 Genome.fasta merged_nodups.txt
```
**3D-dna run that increases stringency for misjoin detection for less fragmentation**
```bash
module load dna_3d;module load parallel;module load java; bash run-asm-pipeline.sh -q 30 --editor-coarse-stringency 90 --editor-repeat-coverage 20 --splitter-coarse-stringency 80 --splitter-coarse-resolution 250000 --editor-coarse-resolution 250000 --editor-fine-resolution 25000 --editor-coarse-region 500000 Genome.fasta merged_nodups.txt
```
**3D-dna run that was the best for this arabidopsis dataset**
```bash
module load dna_3d;module load parallel;module load java; bash run-asm-pipeline.sh --editor-coarse-stringency 90 --editor-coarse-resolution 125000 --editor-coarse-region  250000 --editor-repeat-coverage 30 Genome.fasta merged_nodups.txt
```

### Juicebox evaluation 

**Installation (none, just a download)** <br>
Download the juicebox software that is appropriate for your personal computer. https://github.com/aidenlab/Juicebox/wiki/Download <br>
There is an online version, but if your genome is large or you have lots of reads, it will be very slow to load between zooms.  <br>


**Loading your genome into juicebox** <br>
File, Open, select your HiC file.  Click the plot, otherwise the next step will not work. <br>
Assembly, Import Map Assembly, Click your .assembly file. This may take a few minutes to load. <br>

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

Tip:  If you shift-select a contig and then shift-select another contig and you will lag. Always deselect the contig first by selecting outside of the box.
<br>

**Saving your work (do it often!)** <br>
Juicebox can randomly have bugs and close. Make sure that you are regularly exporting your assembly file.   <br>
Assembly, Export Assembly  -->  the program will automatically rename your .assembly file to .assembly.review. Name it how you like. <br>

**Reloading your unfinished assembly into juicebox for further correction.** <br>
One extra step is required. After loading HiC, clicking the plot, and loading the original unmodified .assembly file, you need to load your modified assembly file after. Assembly, Load Modified Assembly, select your latest genome.assembly.review file. 

### 3D-dna Finalization 

This is where 3D-dna takes all of your modifications in juicebox and recreates the scaffolded genome. If I am sure that I am done with juicebox, I go ahead and sort the output so that the largest scaffolds will begin with the lowest numbers. This adds 100bp gaps between each contig you scaffolded. The '-s' is to set the stage of assembly to start at. We can put indicating that we have already done the juiceboxing and want to try to have the cpu see if it can place any of the remaining smaller contigs. 

```bash

mkdir 3ddnaFinalize
cd 3ddnaFinalize/

module load dna_3d;module load java;module load parallel; bash run-asm-pipeline-post-review.sh -s seal -r Genome.0.review.assembly Genome.fasta merged_nodups.txt

../../02_Files/new_Assemblathon.pl Genome_HiC.fasta 

---------------- Information for assembly 'Genome_HiC.fasta' ----------------
                                         Number of scaffolds       1546
                                     Total size of scaffolds  192727589
                                            Longest scaffold   32344268
                                           Shortest scaffold      15223
                                 Number of scaffolds > 1K nt       1546 100.0%
                                Number of scaffolds > 10K nt       1546 100.0%
                               Number of scaffolds > 100K nt         25   1.6%
                                 Number of scaffolds > 1M nt          5   0.3%
                                Number of scaffolds > 10M nt          5   0.3%
                                          Mean scaffold size     124662
                                        Median scaffold size      32910
                                         N50 scaffold length   22789912
                                          L50 scaffold count          4
                                                 scaffold %A      31.43
                                                 scaffold %C      18.57
                                                 scaffold %G      18.42
                                                 scaffold %T      31.57
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0
                Percentage of assembly in scaffolded contigs      69.8%
              Percentage of assembly in unscaffolded contigs      30.2%
                      Average number of contigs per scaffold        1.0
Average length of break (>25 Ns) between contigs in scaffold        500
                                           Number of contigs       1561
                              Number of contigs in scaffolds         21
                          Number of contigs not in scaffolds       1540
                                       Total size of contigs  192720089
                                              Longest contig   22927795
                                             Shortest contig      15223
                                   Number of contigs > 1K nt       1561 100.0%
                                  Number of contigs > 10K nt       1561 100.0%
                                 Number of contigs > 100K nt         37   2.4%
                                   Number of contigs > 1M nt         13   0.8%
                                  Number of contigs > 10M nt          6   0.4%
                                            Mean contig size     123459
                                          Median contig size      32976
                                           N50 contig length    9415973
                                            L50 contig count          7
                                                   contig %A      31.44
                                                   contig %C      18.57
                                                   contig %G      18.42
                                                   contig %T      31.57
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```



