#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -V

 
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -m ae
#PBS -q batch
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
NP=$(wc -l $PBS_NODEFILE | awk '{print $1}')
echo NP is $NP

#########################################################################


for i in `cat /wtdai/Benchmark4/Domain.ID`
do

	/SPR/Script/Profile.bash  "$i".fasta

done

#########################################################################
