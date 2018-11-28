

#bsub -J "tophat" -q long -W 18:00 -e "tophat.err" -o "tophat.out" tophat -p 8 -G ToxoDB35.gtf -o thout ToxoDB35 B2-P11-intra-bio-rep-2_S5_R1_val_1.fq B2-P11-intra-bio-rep-2_S5_R2_val_2.fq

# This script maps the reads to the reference genome using tophat ###################

#tophat -p 8 -G genes.gtf -o C1_R1_thout genome Cond2_Rep3_R1.fastq Cond2_Rep3_R2.fastq

module unload tophat
module unload tophat/2.0.1
module load tophat/2.0.14
module load samtools/1.4.1
module load bowtie2/2.3.4.3

###########  Job Submission Parameter ###############################################

PROCESSORS=8
MEMORY="49152"       # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32\
DURATION="24:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days
QUEUE="long"        # short = max 4 hours;   long = max 30 days

############  Input and Output Paths ################################################

BOWTIE_INDEX="/project/umb_kourosh_zarringhalam/toxo_new/Genome/bowtie_index_toxoplasma/ToxoDB35"
GTF_FILE="/project/umb_kourosh_zarringhalam/toxo_new/Genome/bowtie_index_toxoplasma/ToxoDB35.gtf"
SAMPLE_FILE="/project/umb_kourosh_zarringhalam/toxo_new/Genome/Fastq_run3_toxoplasm/sample.list.txt"
FASTQ_TRIMED_FILES="/project/umb_kourosh_zarringhalam/toxo_new/Genome/trim_run3/"
OUTPUT_DIR="/project/umb_kourosh_zarringhalam/toxo_new/Genome/tophat_run3/"


for sample in $(cat ${SAMPLE_FILE}); do
  SAMPLE_R1="${sample}"_R1_val_1.fq
  SAMPLE_R2="${sample}"_R2_val_2.fq
  echo ${SAMPLE_R1} ${SAMPLE_R2} 
  echo ${sample}_thout
  bsub -J "th_run3" -q long -W 24:00 -n 8 -R "rusage[mem=49152] span[hosts=1]" tophat -p ${PROCESSORS} -G ${GTF_FILE} -o ${OUTPUT_DIR}${sample}_thout ${BOWTIE_INDEX} ${FASTQ_TRIMED_FILES}${SAMPLE_R1} ${FASTQ_TRIMED_FILES}${SAMPLE_R2}
  
  
done









