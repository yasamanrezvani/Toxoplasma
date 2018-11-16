#!/bin/bash                                                                                                                                                                           

#### This script removes the adaptor of fastq files (R1 and R2)                                                                                                                         
#                                                                                                                                                                                     
###########  loaoding Modules ####################################################                                                                                                                                                                 
module unload trim_galore cutadapt
module load trim_galore/0.4.2 cutadapt/1.9


###########  Job Submission Parameter ############################################

#PROCESSORS=4
#MEMORY="1024"
#DURATION="24:00"
#QUEUE="long"


############  Input and Output Paths #############################################

FASTQ_FILE_PATH="/project/umb_kourosh_zarringhalam/toxo_new/Genome/"
FASTQ_FILE_DIR="Fastq_run3_toxoplasm"
#OUTPUT_DIR=${FASTQ_FILE_PATH}${FASTQ_FILE_DIR}/
OUTPUT_DIR="/project/umb_kourosh_zarringhalam/toxo_new/Genome/trim_bsub_out"

cd ${FASTQ_FILE_PATH}${FASTQ_FILE_DIR}/

############  Sample Names in each ILLUMINA RUN  ##################################                                                                                                     

for item in `ls`; do echo ${item%_R*} ; done | uniq  | sed '$d' > sample.list.txt


####################### trimming ##################################################                                                                                                       
##### Adaptors #######
AD1="AGATCGGAAGAGCACACGTCTGAACTCCAGTC"
AD2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

file="sample.list.txt"
while IFS= read -r line
do
    R1file="$line"_R1.fastq
    R2file="$line"_R2.fastq
    echo ${R1file} ${R2file}
    echo trim_galore --paired -a ${AD1} -a2 ${AD2} $R1file $R2file
    bsub -J "trim" -q long -n 8 -o "trim.out" -e "trim.err" trim_galore --paired -a $AD1 -a2 $AD2 $R1file  $R2file -o $OUTPUT_DIR
    
done < "$file"

####################################################################################
#cd ${OUTPUT_DIR}                                                                                                                                                                    
#exec bash   


