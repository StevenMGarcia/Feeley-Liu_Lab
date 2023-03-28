#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o ~/log_cellranger/                        #-- output directory (fill in)
#$ -e ~/log_cellranger/                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=50G
##$ -l scratch=50G
#$ -l h_rt=72:00:00
#$ -m ea                           #--email when done
#$ -M steven.garcia@ucsf.edu        #--email


ulimit -n 16000

echo "Sample run: " $1



module load CBI
module load cellranger/5.0.1

transcriptome=/../refdata-gex-GRCh38-2020-A

fastq=/../$1
ls $transcriptome

cellranger count \
--id="$1"_cellranger \
--fastqs=$fastq/lane1/,$fastq/lane2/,$fastq/lane3/,$fastq/lane4/  \
--sample=$1 \
--transcriptome=$transcriptome \
--jobmode=/../sge.template
#--chemistry=SC5P-PE

#cp -r "$1"_SARSwhole ../$1

qstat -j $JOB_ID
