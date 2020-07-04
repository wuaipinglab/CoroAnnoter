#!/bin/bash
meme=/home/waplab/meme/bin/meme

input_path=$1

for i in $(ls ${input_path}/*.TRS.fasta)
do
echo $i
samplename=`basename $i|sed 's/.TRS.fasta//g'`
echo ${samplename}
meme ${input_path}/${samplename}.TRS.fasta  -dna -mod zoops -objfun classic  -nmotifs 5 -minw 6 -maxw 8 -minsites 5 -nostatus -time 18000  -markov_order 0 -oc ${input_path}/${samplename}
done


