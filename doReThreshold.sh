#!/bin/bash
rdaFile=$1

newFile=NewThreshold/$(basename $rdaFile | sed 's/.rda/__ReThreshold_SOX10.rda/')

Rscript --no-save reThreshold.R $rdaFile
Rscript --no-save reassignThresholds.R $rdaFile
Rscript --no-save getCellTypeCombinations.R $newFile
