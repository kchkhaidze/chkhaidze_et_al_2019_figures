#!/bin/bash

set -e
set -o pipefail

for id in {1..3}
do
echo "
#!/bin/bash

#BSUB -J "patient_$id"
#BSUB -P DMPPHDAAR
#BSUB -W 48:00
#BSUB -o out_organoids_$id
#BSUB -e out_organoids_$id
#BSUB -n 10
#BSUB -R "span[hosts=1]"

DATADIR=/scratch/DMP/EVGENMOD/kchkhaidze/simulations/CHESS.cpp/

module load java
module load R/3.5.0
module load cmake/3.9.6
module load gcc/7.2.0-abi-compat

export R_LIBS_USER=/scratch/DMP/EVGENMOD/kchkhaidze/R/3.5.0:$R_LIBS_USER

Rscript run_inference_sc_organoids.R $id 500 20 10 500" | bsub

done
