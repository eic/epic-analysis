#!/bin/bash

# software paths
export LARGEXASYM_HOME=$(dirname $(realpath $0))
export BRUFIT=${LARGEXASYM_HOME}/deps/brufit
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${BRUFIT}/lib

# PROOF logs (your path may be different)
jobdir=$(echo $LARGEXASYM_HOME | sed 's,'"$HOME/"',,' | sed 's,\/,-,g')
export PROOF_LOG=${HOME}/.proof/${jobdir}/last-lite-session

# print results
env|grep --color -w LARGEXASYM_HOME
env|grep --color -w BRUFIT
env|grep --color -w LD_LIBRARY_PATH
env|grep --color -w PROOF_LOG

# brufit alias
alias brufit="root $BRUFIT/macros/LoadBru.C"
