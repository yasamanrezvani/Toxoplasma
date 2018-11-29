#!/bin/bash



## This script performs quality control on fastq files 

module load fastqc/0.10.1

FASTQ_FILE_DIR="/project/umb_kourosh_zarringhalam/toxo_new/Genome/FASTQ_run6_toxoplasma/"
OUTPUT_DIR="/project/umb_kourosh_zarringhalam/toxo_new/Genome/fastqc_run6/"


file="/project/umb_kourosh_zarringhalam/toxo_new/Genome/FASTQ_run6_toxoplasma/sample.list.txt"

cd ${FASTQ_FILE_DIR}


for item in `ls *fastq`
do
echo ${item}

fastqc -t 5 ${item} -o ${OUTPUT_DIR}

done

