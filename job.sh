#!/bin/bash

cd /work/clas12/users/$USER/largex-eic
echo "root -q -b macro/dis-5x41/analysis_resolution_SD.C" | ./container/shell.sh
echo "root -q -b macro/dis-5x41/postprocess_resolution_SD.C" | ./container/shell.sh
echo DONE
