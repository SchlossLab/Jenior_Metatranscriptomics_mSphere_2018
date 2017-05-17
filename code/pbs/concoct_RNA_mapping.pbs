#!/bin/sh
#PBS -N concoct_RNA_mapping
#PBS -q fluxod
#PBS -l qos=flux
#PBS -l nodes=1:ppn=4,mem=32GB
#PBS -l walltime=300:00:00
#PBS -j oe
#PBS -V
#PBS -M mljenior@umich.edu
#PBS -A pschloss_fluxod


echo "ncpus-2.pbs"
cat $PBS_NODEFILE
qstat -f $PBS_JOBID

cd $PBS_O_WORKDIR

NCPUS=`wc -l $PBS_NODEFILE | awk '{print $1}'`


# Load samtools
module load samtools

# Align paired-end reads to annotated CONCOCT clusters
#/nfs/turbo/pschloss/bin/bowtie2 -f --fr -p 8 -x ${cluster}_gene_db \
#        -1 /nfs/turbo/pschloss/matt/metatranscriptomes_HiSeq/final_reads/${metatranscriptome}.read1.pool.cut.trim.silva_filter.mus_filter.fasta \
#        -2 /nfs/turbo/pschloss/matt/metatranscriptomes_HiSeq/final_reads/${metatranscriptome}.read2.pool.cut.trim.silva_filter.mus_filter.fasta \
#        -S ${cluster}_DNA.${metatranscriptome}_RNA.pair.sam
#samtools view -bS ${cluster}_DNA.${metatranscriptome}_RNA.pair.sam > ${cluster}_DNA.${metatranscriptome}_RNA.pair.bam
#samtools sort ${cluster}_DNA.${metatranscriptome}_RNA.pair.bam ${cluster}_DNA.${metatranscriptome}_RNA.pair.sort
#rm ${cluster}_DNA.${metatranscriptome}_RNA.pair.sam

# Align orphaned reads to Full length C. difficile 630 genome and plasmid
#/nfs/turbo/pschloss/bin/bowtie2 -f -p 8 -x ${cluster}_gene_db \
#        -U /nfs/turbo/pschloss/matt/metatranscriptomes_HiSeq/final_reads/${metatranscriptome}.orphan.pool.cut.trim.silva_filter.mus_filter.fasta \
#        -S ${cluster}_DNA.${metatranscriptome}_RNA.orphan.sam
#samtools view -bS ${cluster}_DNA.${metatranscriptome}_RNA.orphan.sam > ${cluster}_DNA.${metatranscriptome}_RNA.orphan.bam
#samtools sort ${cluster}_DNA.${metatranscriptome}_RNA.orphan.bam ${cluster}_DNA.${metatranscriptome}_RNA.orphan.sort
#rm ${cluster}_DNA.${metatranscriptome}_RNA.orphan.sam

# Merge alignments and convert to idxstats format
#samtools merge ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.bam ${cluster}_DNA.${metatranscriptome}_RNA.pair.sort.bam ${cluster}_DNA.${metatranscriptome}_RNA.orphan.sort.bam
#samtools sort ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.bam ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort

# Screen duplicate headers
samtools view -h -o ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.sam ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.bam
python /nfs/turbo/pschloss/matt/home/scripts/rmdup_SQ.py ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.sam ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.rmsq.sam
samtools view -bS ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.rmsq.sam > ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.rmsq.bam
rm ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.sam ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.rmsq.sam

# Remove PCR duplicates and sort
java -Xmx2g -jar /nfs/turbo/pschloss/matt/bin/picard-tools-1.119/MarkDuplicates.jar \
        INPUT=${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.rmsq.bam \
        OUTPUT=${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.rmsq.rmdup.bam \
        METRICS_FILE=${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.rmsq.rmdup.metrics \
        AS=TRUE \
        VALIDATION_STRINGENCY=LENIENT \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        REMOVE_DUPLICATES=TRUE
samtools index ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.rmsq.rmdup.bam

# Convert to human-readable format
samtools idxstats ${cluster}_DNA.${metatranscriptome}_RNA.sort.merge.sort.rmsq.rmdup.bam > ${cluster}_DNA.${metatranscriptome}_RNA.final.txt

# Format and normalize idxstats files
python /nfs/turbo/pschloss/matt/home/scripts/idxstats/pool_unmapped_idxstats.py ${cluster}_DNA.${metatranscriptome}_RNA.final.txt
python /nfs/turbo/pschloss/matt/home/scripts/idxstats/normalize_idxstats.py ${cluster}_DNA.${metatranscriptome}_RNA.final.pool.txt 50
python /nfs/turbo/pschloss/matt/home/scripts/idxstats/get_ko.py ${cluster}_DNA.${metatranscriptome}_RNA.final.pool.norm.txt


echo "qsub working directory absolute is"
echo $PBS_O_WORKDIR
exit
