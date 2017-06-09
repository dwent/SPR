#!/bin/bash

path=../fasta

for i in `cat ../../ex.id`
do

	cp  "$path"/"$i".fasta  ./

done
