#!/bin/bash
# Create full simulation resolution plotting scripts and submit to ifarm
script='/work/clas12/users/mfmce/largex-eic/macro/analysis_resolution_SD_Full.C'
postScript='/work/clas12/users/mfmce/largex-eic/macro/postprocess_resolution_SD.C'
submitScript='/work/clas12/users/mfmce/largex-eic/submit.sh'
jobScript='/work/clas12/users/mfmce/largex-eic/job.sh'
out='/work/clas12/users/mfmce/largex-eic/macro/canyonlands/'
cd /work/clas12/users/mfmce/largex-eic/datarec/canyonlands/5x41/
for file in /work/clas12/users/mfmce/largex-eic/datarec/canyonlands/*/*.config
do
    config=`echo $file | sed "s;.config;;g" | grep -Eo 'canyonlands/.*' | sed 's;canyonlands/;;g' | sed 's;/files;;g'`
    mkdir -p $out/$config
    cp $script $out/$config/
    newscript=${out}${config}/*.C
    eleIn=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*" | sed "s;x.*;;g"`
    beamIn=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*" | sed "s;.*x;;g"`
    xAng=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*-x[0-9]*" | sed "s;.*-x;;g"`
    xAngM=`echo $file | grep -Eo "[0-9][0-9]*x[0-9]*-xm[0-9]*" | sed "s;.*-xm;;g"`
    echo "name=${name}"
    echo "file=${file}"
    echo "config=${config}"
    echo "newscript=${newscript}"
    echo "eleIn=${eleIn}"
    echo "beamIn=${beamIn}"
    echo "xAng=${xAng}"
    echo "xAngM=${xAngM}"
    echo --------------------

    sed -i "s;datarec/dis-5x41.config;${file};g" $newscript
    sed -i "s;full_dis-5x41;full_${config};g" $newscript
    sed -i "s;eleBeamEn=5;eleBeamEn=${eleIn};g" $newscript
    sed -i "s;ionBeamEn=41;ionBeamEn=${beamIn};g" $newscript
    if [ $xAng ]; then
        sed -i "s;crossingAngle=0;crossingAngle=${xAng};g" $newscript
    fi
    if [ $xAngM ]; then
        sed -i "s;icrossingAngle=0;crossingAngle=-${xAngM};g" $newscript
    fi
    cp $postScript $out/$config
    sed -i "s;out/resolution;out/full_${config};g" $out/$config/*.C
    if [ $xAng ]; then
	sed -i "s;testheader;${eleIn}x${beamIn} ${xAng};g" $out/$config/*.C
    fi
    if [ $xAngM ]; then
        sed -i "s;testheader;${eleIn}x${beamIn} -${xAngM};g" $out/$config/*.C
    fi
    cp $submitScript $out/$config
    cp $jobScript $out/$config
    sed -i "s;dis-5x41;canyonlands/${config};g" $out/$config/*.sh
    
    # submit job
    sbatch $out/$config/submit.sh
    
done
cd ..
echo DONE
#cd /work/clas12/users/mfmce/largex-eic
#echo "root -q -b macro/analysis_resolution_SD.C" | ./container/shell.sh
#echo "root -q -b macro/postprocess_resolution_SD.C" | ./container/shell.sh
#echo DONE
