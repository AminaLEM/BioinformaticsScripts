#!/bin/bash

module load deeplexicon
#. /beegfs/biosw/deeplexicon/1.0.0/venv/bin/activate
deeplexicon.py \
        -p /prj/Isabel_ONT_rRNA/YeastData+Results/Data/FAST5/RNA814001/workspace \
        -f multi \
        -m /beegfs/biosw/deeplexicon/1.0.0/deeplexicon/models/resnet20-final.h5 > demulti_outputlong3.tsv

