#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N blastp
module load GNU/4.8.3
module load BLAST+/2.3.0

cd $PBS_O_WORKDIR

tblastn -query $seq  -db $DB -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 1 > $output

qstat -f ${PBS_JOBID}
