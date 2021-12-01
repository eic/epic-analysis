#!/bin/bash

# Create fast simulation resolutions plotting scripts and submit to ifarm
METHOD="Ele" # Switch this to select reconstruction method from {"Ele","JB","DA"}
script="$PWD/macro/analysis_resolution.C"
postScript="$PWD/macro/postprocess_resolution.C"
submitScript="$PWD/submit.sh"
jobScript="$PWD/job.sh"
out="$PWD/macro/"
cd datarec
for file in *-xm25.config
do
    energies=`echo $file | sed "s;.config;;g"`
    config="${METHOD}_${energies}"
    mkdir -p $out/$config
    cp $script $out/$config/
    newscript=${out}${config}/*.C
    eleIn=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*" | sed "s;x.*;;g"`
    beamIn=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*" | sed "s;.*x;;g"`
    xAng=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*-x[0-9]*" | sed "s;.*-x;;g"`
    xAngM=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*-xm[0-9]*" | sed "s;.*-xm;;g"`
    echo "file=${file}"
    echo "energies=${energies}"
    echo "config=${config}"
    echo "newscript=${newscript}"
    echo "eleIn=${eleIn}"
    echo "beamIn=${beamIn}"
    echo "xAng=${xAng}"
    echo "xAngM=${xAngM}"

    # Modify analysis script
    sed -i "s;datarec/dis-5x41;datarec/${energies};g" $newscript
    sed -i "s;Ele_dis-5x41;${config};g" $newscript
    sed -i "s;eleBeamEn=5;eleBeamEn=${eleIn};g" $newscript
    sed -i "s;ionBeamEn=41;ionBeamEn=${beamIn};g" $newscript
    if [ $xAng ]; then
        sed -i "s;crossingAngle=25;crossingAngle=-${xAng};g" $newscript
    fi
    if [ $xAngM ]; then
        sed -i "s;crossingAngle=25;crossingAngle=${xAngM};g" $newscript
    fi
    sed -i "s;Ele;${METHOD};g" $newscript

    # Postprocessor
    cp $postScript $out/$config
    sed -i "s;dis-5x41;${config};g" $out/$config/postprocess*.C
	sed -i "s;testheader;${eleIn}x${beamIn}GeV;g" $out/$config/postprocess*.C

    # And job scripts
    cp $submitScript $out/$config
    cp $jobScript $out/$config
    sed -i "s;dis-5x41;${config};g" $out/$config/*.sh

    # And submit
	sbatch $out/$config/submit.sh 
    echo --------------------
 
done
cd ..
echo DONE
