#!/bin/bash

# Create fast simulation resolutions plotting scripts and submit to ifarm

script="$PWD/macro/analysis_resolution_ELE.C"
postScript="/$PWD/macro/postprocess_resolution_ELE.C"
submitScript="$PWD/submit.sh"
jobScript="$PWD/job.sh"
out="$PWD/macro/"
cd datarec
for file in *-xm25.config
do
    config=`echo $file | sed "s;.config;;g"`
    mkdir -p $out/$config
    cp $script $out/$config/
    newscript=${out}${config}/*.C
    eleIn=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*" | sed "s;x.*;;g"`
    beamIn=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*" | sed "s;.*x;;g"`
    xAng=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*-x[0-9]*" | sed "s;.*-x;;g"`
    xAngM=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*-xm[0-9]*" | sed "s;.*-xm;;g"`
    echo "file=${file}"
    echo "config=${config}"
    echo "newscript=${newscript}"
    echo "eleIn=${eleIn}"
    echo "beamIn=${beamIn}"
    echo "xAng=${xAng}"
    echo "xAngM=${xAngM}"

    sed -i "s;dis-5x41;${config};g" $newscript
    sed -i "s;eleBeamEn=5;eleBeamEn=${eleIn};g" $newscript
    sed -i "s;ionBeamEn=41;ionBeamEn=${beamIn};g" $newscript
    if [ $xAng ]; then
        sed -i "s;crossingAngle=25;crossingAngle=-${xAng};g" $newscript
    fi
    if [ $xAngM ]; then
        sed -i "s;crossingAngle=25;crossingAngle=${xAngM};g" $newscript
    fi
    cp $postScript $out/$config
    sed -i "s;dis-5x41;${config};g" $out/$config/*.C
	sed -i "s;testheader;${eleIn}x${beamIn}GeV;g" $out/$config/*.C
    cp $submitScript $out/$config
    cp $jobScript $out/$config
    sed -i "s;dis-5x41;${config};g" $out/$config/*.sh
	sbatch $out/$config/submit.sh 
    echo --------------------
 
done
cd ..
echo DONE
