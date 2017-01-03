#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
# resolve $SOURCE until the file is no longer a symlink
while [ -h "$SOURCE" ]; do
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  # if $SOURCE was a relative symlink, we need to resolve it relative to
  # the path where the symlink file was located
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
CUR_DIR=$(pwd)

cd ${DIR}/others

# check TOFCS package
if [ -d 'TFOCS' ]; then
  git pull
else
  git clone git@github.com:cvxr/TFOCS.git
fi

# get FPC_AS at the following address:
# http://www.caam.rice.edu/~optimization/L1/FPC_AS/request-for-downloading-fpc_as.html

# get glmnet_matlab at the following address:
# http://web.stanford.edu/~hastie/glmnet_matlab/download.html


exit 0;
