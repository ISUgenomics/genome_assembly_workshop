

## RepeatModeler:

* Models transposon repeat families	
```bash
module load repeatmodeler/2.0.5
BuildDatabase -engine ncbi -name ATNDB chr2.fa 
RepeatModeler -database ATNDB -engine ncbi -threads 36 >& rm.log&
```
## Repeat masking:

```bash
RepeatMasker -pa 6 -gff -nolow -engine ncbi -lib consensi.fa chr2.fa
```

## Braker:

```bash
module load braker
module load augustus
cp -r /software/el9/apps/augustus/3.5.0/config/ config
AUGUSTUS_CONFIG_PATH="$PWD/config/"

braker.pl --threads 20 \
--AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
--species=Athaliana \
--genome=chr2.fa.masked \
--bam=SRR26158835_chr2.bam \
--prot_seq=chr2_proteins.fasta \
--gff3
```
