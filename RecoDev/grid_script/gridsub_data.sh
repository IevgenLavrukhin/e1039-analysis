#!/bin/bash

dir_macros=$(dirname $(readlink -f $BASH_SOURCE))
LIFE_TIME=long # short (3h), medium (8h) or long (23h)

jobname=$1
do_sub=$2
run_list=$3
nevents=$4
dst_mode=${5:-'split'} # 'split' or 'single'

echo "njobs=$njobs"
echo "nevents=$nevents"

 
if [ $do_sub == 1 ]; then
    echo "Grid mode."
    if ! which jobsub_submit &>/dev/null ; then
	echo "Command 'jobsub_submit' not found."
	echo "Forget 'source /e906/app/software/script/setup-jobsub-spinquest.sh'?"
	exit
    fi
    work=/pnfs/e906/persistent/users/$USER/RecoDataDev/$jobname
    ln -sf /pnfs/e906/persistent/users/$USER/RecoDataDev data
else
    echo "Local mode."
    work=$dir_macros/scratch/$jobname
fi

mkdir -p $work
chmod -R 01755 $work

cd $dir_macros/..
tar -czvf $work/input.tar.gz  support RecoE1039Data.C
cd -

for run_num in $(cat $run_list) ; do  
  #location of the data in your persistent area
  data_dir="/pnfs/e906/persistent/users/$USER/DstRun"
  declare -a data_path_list=()
  if [ $dst_mode = 'single' ] ; then
      data_path_list=( $data_dir/$(printf 'run_%06d_spin.root' $run_num) )
  else # 'split'
      data_path_list=( $(find $data_dir -name $(printf 'run_%06d_spill_*_spin.root' $run_num) ) )
  fi

  for data_path in ${data_path_list[*]} ; do
    data_file=$(basename $data_path)
    job_name=${data_file%'.root'}
    mkdir -p $work/$job_name/log
    mkdir -p $work/$job_name/out
    chmod -R 01755 $work/$job_name
    
    rsync -av $dir_macros/gridrun_data.sh $work/$job_name/gridrun_data.sh

    if [ $do_sub == 1 ]; then
      cmd="jobsub_submit --grid"
      cmd="$cmd -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/e1039/e1039-sl7:latest\"'"
      cmd="$cmd --append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)'"
      cmd="$cmd --use_gftp --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE -e IFDHC_VERSION --expected-lifetime='$LIFE_TIME'"
      cmd="$cmd --mail_never"
      cmd="$cmd -L $work/$job_name/log/log.txt"
      cmd="$cmd -f $work/input.tar.gz"
      cmd="$cmd -d OUTPUT $work/$job_name/out"
      cmd="$cmd -f $data_path"
      cmd="$cmd file://`which $work/$job_name/gridrun_data.sh` $nevents $run_num $data_file"
      echo "$cmd"
      $cmd
    else
      mkdir -p $work/$job_name/input
      rsync -av $work/input.tar.gz $data_path  $work/$job_name/input
      cd $work/$job_name/
      $work/$job_name/gridrun_data.sh $nevents $run_num $data_file | tee $work/$job_name/log/log.txt
      cd -
    fi
  done
done 2>&1 | tee log_gridsub.txt
