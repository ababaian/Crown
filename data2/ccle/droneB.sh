#!/bin/bash
# droneB.sh
#

# This script-layer is neccesary to launch a screen session
# on each ec2-machine. The pipeline is run within that session
# and the output is logged. This allows 'looking in' on sessions
# as they are running.

# Commands to run on server-side
# ===============================================================

SCRIPTPATH=$1

SCRIPT=$(basename $1)

shift # drop first (TASK or SCRIPT variable)

# Download pipeline / droneB's function
  aws s3 cp $SCRIPTPATH ./

  chmod 777 *.sh

# open screen; run gather.sh function. -L logged
  screen -Ldmt sh ~/$SCRIPT $@


# ===============================================================
