#!/bin/bash

jobname=$1
do_sub=$2

njobs=$3
nevents=$4

if [ $do_sub == 1 ]; then
  echo "grid!"
else
  echo "local!"
fi

echo "njobs=$njobs"
echo "nevents=$nevents"

macros=/e906/app/users/yuhw/e1039-analysis/TargetSim

sed "s/nevents=NAN/nevents=$nevents/"             $macros/gridrun.sh > $macros/gridrun_new.sh 
#sed -i "s/gap=NAN/gap=$gap/"                      $macros/gridrun_new.sh
chmod +x $macros/gridrun_new.sh

if [ $do_sub == 1 ]; then
work=/pnfs/e906/persistent/users/yuhw/TargetSim/$jobname
else
work=$macros/scratch/$jobname
fi

mkdir -p $work
chmod -R 01755 $work

cd $macros
tar -czvf $work/input.tar.gz *.C *.cfg *.opts *.so
cd -

for (( id=1; id<=$njobs; id++ ))
do  
  mkdir -p $work/$id/log
  mkdir -p $work/$id/out
  chmod -R 01755 $work/$id

  rsync -av $macros/gridrun_new.sh $work/$id/gridrun_new.sh

  cmd="jobsub_submit --grid"
  cmd="$cmd -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/e1039/e1039-sl7:latest\"'"
  cmd="$cmd --append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)'"
  cmd="$cmd --use_gftp --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE -e IFDHC_VERSION --expected-lifetime='short'"
  cmd="$cmd --mail_never"
  cmd="$cmd -L $work/$id/log/log.txt"
  cmd="$cmd -f $work/input.tar.gz"
  cmd="$cmd -d OUTPUT $work/$id/out"
  cmd="$cmd file://`which $work/$id/gridrun_new.sh`"

  if [ $do_sub == 1 ]; then
    echo $cmd
    $cmd
  else
    mkdir -p $work/$id/input
    rsync -av $work/input.tar.gz $work/$id/input
    cd $work/$id/
    $work/$id/gridrun_new.sh | tee $work/$id/log/log.txt
    cd -
  fi
done
