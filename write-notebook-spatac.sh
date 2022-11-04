#!/bin/bash

SERVERPORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()') ; echo $SERVERPORT

RUNNING_DIR=/lustre/scratch/kiviaho/spatac

module load CUDA/11.2
module load anaconda
source activate spatac


jupyter notebook \
--notebook-dir=${RUNNING_DIR} --NotebookApp.token='passatac' \
--ip=0.0.0.0 --port=$SERVERPORT --no-browser --allow-root
