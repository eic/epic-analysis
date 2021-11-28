#!/bin/bash

cd /work/clas12/users/$USER/eic3/largex-eic
echo "root -q -b macro/dis-5x41/analysis_resolution.C" | ./container/shell.sh
echo "root -q -b macro/dis-5x41/postprocess_resolution.C" | ./container/shell.sh
echo DONE
