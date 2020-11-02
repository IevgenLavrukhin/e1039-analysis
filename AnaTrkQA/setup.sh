DIR_TOP=$(dirname $(readlink -f $BASH_SOURCE))

source /e906/app/software/osg/software/e1039/this-e1039.sh
#source /e906/app/software/osg/users/$USER/e1039/core/this-e1039.sh
export LD_LIBRARY_PATH=$DIR_TOP/inst/lib:$LD_LIBRARY_PATH

function cmake-ana-TrkQA {
    cmake -DCMAKE_INSTALL_PREFIX=$DIR_TOP/inst $DIR_TOP/src
    ret=$?
    test $ret -eq 0 && echo "OK, run 'make install'."
}

if [ ${HOSTNAME:0:13} != 'spinquestgpvm' ] ; then
    echo "!!CAUTION!!"
    echo "This AnaTrkQA package does not support your computer ($HOSTNAME)."
    echo "It might not run properly even if you follow 'README.md'."
    echo
fi
