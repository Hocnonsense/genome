#!/bin/bash

cd $CONDA_PREFIX/lib || {
    echo "Error: Failed to change to $CONDA_PREFIX/lib" >&2
    exit 1
}

    (
        git clone https://github.com/motu-tool/fetchMGs.pl.git
        cd fetchMGs.pl && git checkout v1.2
    ) || {
        echo "Error: Failed to clone or checkout fetchMGs.pl repository" >&2
        exit 1
    }

    # Create wrapper script
    echo "$(realpath fetchMGs.pl/fetchMGs.pl)" '$@' > ../bin/fetchMGs
    chmod +x ../bin/fetchMGs
