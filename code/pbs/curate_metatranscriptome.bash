#!/bin/bash


## Set up output directories
mkdir /mnt/EXT/Schloss-data/matt/metatranscriptomes_HiSeq/processing/
mkdir /mnt/EXT/Schloss-data/matt/metatranscriptomes_HiSeq/curated/


## Assign correct Nextera XT primer pair that was used (5' + 3' reverse complement)
F_primer=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
if [ $1 = 'cefoperazone_630' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='infected'
elif [ $1 = 'cefoperazone_mock' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='mock'
elif [ $1 = 'clindamycin_630' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='infected'
elif [ $1 = 'clindamycin_mock' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='mock'
elif [ $1 = 'streptomycin_630' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='infected'
elif [ $1 = 'streptomycin_mock' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATCGTTTCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='mock'
elif [ $1 = 'conventional' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATAAGGCCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='mock'
elif [ $1 = 'germfree' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATATCCACTCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='infected'
fi

## Pool reads from the same lane
cd /mnt/EXT/Schloss-data/matt/metatranscriptomes_HiSeq/fastq/$1
pool_job_id=$(qsub -v transcriptome=$1 /mnt/EXT/Schloss-data/matt/metatranscriptomes_HiSeq/pbs/pool.pbs | sed 's/\..*$//')
echo $1 read pooling: $pool_job_id

## Curate and filter contaminating reads
cd /mnt/EXT/Schloss-data/matt/metatranscriptomes_HiSeq/processing/
trimming_job_id=$(qsub -v transcriptome=$1,forward=$F_primer,reverse=$R_primer -W depend=afterok:$pool_job_id /mnt/EXT/Schloss-data/matt/metatranscriptomes_HiSeq/pbs/trimming.pbs | sed 's/\..*$//')
echo $1 quality trimming: $trimming_job_id
filter_job_id=$(qsub -v transcriptome=$1,group=$infection,count=$2 -W depend=afterok:$trimming_job_id /mnt/EXT/Schloss-data/matt/metatranscriptomes_HiSeq/pbs/filter.pbs | sed 's/\..*$//')
echo $1 filtering reads: $filter_job_id


echo