#!/bin/zsh
## run the awk script and generate a test input set
DATAPATH=$HOME/Projects/urqmd-extended
AWKSCRIPT=./src/pt-phi-grab.awk
awk -f $AWKSCRIPT < $DATAPATH/test-analysis-data/pilot-run-0-partial/urqmd_output/* > flow-test-input.dat
