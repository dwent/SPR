#!/bin/bash 


#################################################

targpath=/daiwentao/Database/DockingBenchmark/SingleDomain

for i in `cat /daiwentao/Database/DockingBenchmark/Domain.ID`
do
	
	cp "$targpath"/"$i".ent .
	./naccess  "$i".ent
	
	rm "$i".ent
	rm "$i".log
		
done

#################################################
