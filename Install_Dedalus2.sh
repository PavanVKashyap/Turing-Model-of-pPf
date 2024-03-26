#!/usr/bin/env bash

## Install Dedalus v2

echo "Installing Dedalus v2" 

conda create -n dedalus2

conda activate dedalus2

if [ "$(uname)" == "Darwin" ]; then
    # Do something under Mac OS X platform
    echo "Additional configuration being added in Apple system"
    conda config --env --set subdir osx-64    
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Do something under GNU/Linux platform
    ## NOthing special needs to be done
    echo "No additional settings necessary for Linux system"
fi

conda env config vars set OMP_NUM_THREADS=1
conda env config vars set NUMEXPR_MAX_THREADS=1

conda install -c conda-forge dedalus=2.2207.3

## Activate the dedalus environment
conda activate dedalus2

## Install PyQt
echo "Installing PyQt5"
conda install pyqt

## Install sympy
echo "Installing sympy"
conda install sympy

