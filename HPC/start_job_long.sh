# /bin/bash !

fp=`pwd` 
job=${fp##*/}

#echo $job

bsub -W 120:00 -R "rusage[mem=12000]" -J $job "./DAC > experiment.log"
