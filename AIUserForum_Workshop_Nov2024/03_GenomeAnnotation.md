

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
#SBATCH -N1
#SBATCH -n24
#SBATCH -t 12:00:00
#SBATCH -o "Braker_stdout"
#SBATCH -e "Braker_stderr"
# loading modules
module load braker
module load augustus
# copying config directory to a place where you have write permission
cp -r /apps/spack-managed/gcc-11.3.1/augustus-3.5.0-iuybkzt5zpfde5fk6rscgs34ihug7htm/config .
AUGUSTUS_CONFIG_PATH="$PWD/config/"
# Braker run
braker.pl --threads 20 \
--AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
--species=Athaliana \
--genome=chr2.fa.masked \
--bam=SRR26158835_chr2.bam \
--prot_seq=chr2_proteins.fasta \
--gff3
```

## Extracting proteins and transcripts
```bash
gffread braker.gff3 -g 02_References/chr2.fa -y proteins.fa

gffread braker.gff3 -g 02_References/chr2.fa -w transcripts.fa
```
