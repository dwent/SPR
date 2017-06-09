#!/bin/bash 
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -m ae
#PBS -q batch 


#########################


pdbpath=/daiwentao/Database/DockingBenchmark/SingleDomain
outpath=/daiwentao/Database/DockingBenchmark/Solvation

for i in `cat /daiwentao/Database/DockingBenchmark/Domain.ID`
do
	
	/SPR/Script/indep/wesol2  "$pdbpath"/"$i".ent  |  tee  "$outpath"/"$i".sol

done
