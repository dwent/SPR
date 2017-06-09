#!/bin/bash 


#################################################

targpath=../pdb/

for i in `cat ../../ex.id`
do
	
	cp "$targpath"/"$i".ent .
	./naccess  "$i".ent
	
	rm "$i".ent
	rm "$i".log
		
done

#################################################
