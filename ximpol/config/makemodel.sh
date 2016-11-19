#!/bin/bash
for var in "$@"
do
	   filename=${var##*/}
	   echo $filename
	   header=${filename%%.*}
	   sed "s/10qedon.txt/$filename/" < template_herx1.py > herx1_pulsar${header/qed/}.py
done
