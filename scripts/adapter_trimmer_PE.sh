#!/bin/bash -login
#PBS -l walltime=02:00:00,nodes=1:ppn=1,mem=24Gb		## 2h x 2 Gb
#mdiag -A ged	
#PBS -m abe			#send email to myself
#PBS -N T_Trim		#give name to the job


module load Trimmomatic 

#temp=$(basename "$R1_INPUT")
#new_dir=${temp%_R1_*}

#new_dir=$(basename "$R1_INPUT" | cut -f 1,3 -d "_")
#mkdir ${new_dir}
#cd ${new_dir}

cd $PBS_O_WORKDIR

java -jar $TRIM/trimmomatic PE -threads 1 -phred33 ${R1_INPUT} ${R2_INPUT} ${output_pe1} ${output_se1} ${output_pe2} ${output_se2} ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:20


qstat -f ${PBS_JOBID}

