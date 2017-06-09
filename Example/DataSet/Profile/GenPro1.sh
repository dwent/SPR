#!/bin/bash 
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -m ae
#PBS -q batch
#echo Working directory is $PBS_O_WORKDIR
#cd $PBS_O_WORKDIR
#NP=$(wc -l $PBS_NODEFILE | awk '{print $1}')
#echo NP is $NP

#########################################################################


for i in `cat ../../ex.id`
do

	../../Profile.bash  "$i".fasta

done

#########################################################################
