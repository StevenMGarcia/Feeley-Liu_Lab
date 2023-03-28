#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=30G
#$ -l scratch=100G
#$ -l h_rt=72:00:00
#$ -pe smp 4
#$ -m ea                           #--email when done


sample=$1
bampath=~/$sample

if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

cd "$TMPDIR"

#source ../miniconda3/bin/activate velo

#velocyto --help
velocyto run -b $bampath/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
	-o "$sample"_veloc \
	-@ 8 \
	-m ~/resources/GRCh38_rmsk.gtf \
	$bampath/outs/possorted_genome_bam.bam \
	/../refdata-gex-GRCh38-2020-A/genes/genes.gtf 

cp -r "$sample"_veloc /$bampath/

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
