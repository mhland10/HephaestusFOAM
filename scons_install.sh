#!/bin/bash

# --- Dependencies ---
sudo apt update
sudo apt install -y \
    g++ gfortran scons python3-dev python3-numpy python3-yaml \
    libboost-all-dev libeigen3-dev \
    git