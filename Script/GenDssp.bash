#!/bin/bash


pdbpath=/daiwentao/Database/DockingBenchmark/SingleDomain
targpath=/daiwentao/Database/DockingBenchmark/Dssp


for i in `cat /daiwentao/Database/DockingBenchmark/Domain.ID`
do

	/Script/DSSP/dsspcmbi  "$pdbpath"/"$i".ent  "$targpath"/"$i".dssp

done

