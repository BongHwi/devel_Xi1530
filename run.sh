#!/bin/bash

###########################
## Author : Beomkyu Kim  ##
## email  : kimb@cern.ch ##
###########################

if [[ "$#" != 2 ]]
then
echo "Wrong usage"
##echo "usage : $0 <taskname> <Period> <full|terminate|download|merge>"
echo "usage : $0 <full|terminate>"
exit 0;
fi

##taskname : task name
taskname=Xi1530
#taskname=$1
##periods : remove some periods if you don't want to run all of them
#periods=LHC10d
#periods=LHC16k
periods=$1
method=$2

basepath=/alice/cern.ch/user/$(echo $USER|perl -pe's|.|$&/$&|')/
DOWN_DIR=~/Desktop
DATA_DIR=$DOWN_DIR${basepath}/${taskname}${periods};
currentdir=$PWD


function download {
cd $DOWN_DIR
perl ${ALICE_PHYSICS}/PWGUD/DIFFRACTIVE/macros/alien_cp.pl ${basepath}/${taskname}${periods}/ root_archive.zip AnalysisResults.root
}

function merge_list {
outname=${1:-temp_out.root}
np=${2:-1} # Number Of Process
tag=tmp-$outname-$(date +%s)-$RANDOM
mkdir -p $tag
xargs -n25 | perl -ne'print "'$tag'/${.}.root $_"' | xargs  -P$np -L1   hadd -f
#xargs -n25 | perl -ne'print "${.}.root $_"' | xargs -P$np -I% bash -c 'echo hadd -f '$tag/'%'
hadd -f $outname $tag/*.root
echo $outname
rm $tag/*.root && rmdir $tag
ls
pwd
}
export -f merge_list


######################
#  MERGE Run by Run
#####################
function merge_RunByRun {
cd $DATA_DIR
rm AnalysisResults*.root
#ls -d out/*| xargs -P4 -I% bash -c 'find % -name AnalysisResults.root | merge_list 'AnalysisResults_${taskname}${periods}'_$(basename %).root'
ls -d out/* | xargs -P4 -I% bash -c 'hadd -f AnalysisResults_'${taskname}${periods}'_$(basename %).root $(find % -name AnalysisResults.root)'
hadd -f  AnalysisResults_${taskname}${periods}.root AnalysisResults_${taskname}${periods}_*.root
echo "#------------------------------"
echo "#      Merged Files "
cd - > /dev/null
rm -r  $DOWN_DIR/${taskname}${periods}
mkdir -p $DOWN_DIR/${taskname}${periods}
mv $DATA_DIR/AnalysisResults*.root  $DOWN_DIR/${taskname}${periods}/
echo "#   ./${taskname}${periods}/AnalysisResults_${taskname}${periods}.root"
echo "#------------------------------"
}


if [ $method = "full" ]
then
istart=(`seq 0 5 100`)
iend=(`seq 5 5 105`)
for ijob in "0"
do
for i in ${periods}
do
root -l -b -q run.C\(\"${taskname}\",\"${i}\",\"full\",${istart[$ijob]},${iend[$ijob]}\,\"grid\"\)
rm ${taskname}*
rm *.d *.so
rm *.root
rm *.xml
rm myAnalysis.C
rm stderr stdout
done
##sleep 300
##./resubmit_alien.sh
done

elif [ $method = "download" ]
then
download
elif [ $method = "merge" ]
then
merge_RunByRun
elif [ $method = "terminate" ]
then
download
merge_RunByRun
fi

##cd $currentdir

