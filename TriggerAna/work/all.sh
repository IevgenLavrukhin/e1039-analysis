#!/bin/bash

nJob=$1
nEvent=$2

./gridsub.sh DY_zvtx_spin 1 $nJob $nEvent 0 
./gridsub.sh DY_zvtx_sea 1 $nJob $nEvent 1
