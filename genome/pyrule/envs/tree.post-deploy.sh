#!/bin/bash

cd $CONDA_PREFIX/lib
    git clone https://github.com/motu-tool/fetchMGs.pl.git
    git checkout v1.2

    echo `realpath fetchMGs.pl/fetchMGs.pl` '$@' > ../bin/fetchMGs
    chmod +x ../bin/fetchMGs
