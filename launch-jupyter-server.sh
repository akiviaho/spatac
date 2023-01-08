#!/bin/bash

SERVERPORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()') ; echo ${SERVERPORT}

RUNNING_DIR=/lustre/scratch/kiviaho/spatac

module load anaconda
source activate squidpy


jupyter notebook \
--notebook-dir=${RUNNING_DIR}  \
--ip=0.0.0.0 --port=$SERVERPORT --no-browser --allow-root
