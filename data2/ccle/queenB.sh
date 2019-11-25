#!/bin/bash
# queenB.sh
# 20180814 build
# EC2 Launch / Control Script
#

# 1. queenB script is initialized locally and input files
#    are parsed ready for cluster analaysis
# 2. queenB launches instances, logs in to it and runs the
#    droneB.sh script remotely.
# 3. The droneB script is executed on the instance and it
#    launches a `screen` on the instance and loads and 
#    starts to perform the $TASK (gather.sh) script.
# 4. TASK script should include a instance shut-down
#    command to close instance upon completion.
#

# Control Panel =========================
# Amazon AWS S3 Home URL
S3URL='s3://crownproject/ccle'

# EC2 TASK Script - script for droneB to execute
TASK="$S3URL/scripts/hgr1_align_v4.ccle.sh"

# Parameter file:
# Each line of PARAMETERS will be input to STDIN of
# the droneB script which can then be used to run the
# TASK script.
# i.e. bash droneB.sh <line_N_of_PARAMETERS>
# PARAMETERS="tcga0_input.txt"
PARAMETERS=$1

# EC2 Set-up
instanceTYPE='c4.xlarge'
imageID='ami-0b375c9c58cb4a7a2' #AMI Crown CCLE/SRA

devNAME='/dev/sda1' # /dev/sda1 for Crown-AMI
volSIZE='300' # in Gb

# Maximum number of EC2 instances to run simultaniously
MAX_INSTANCE='45'

# Number of instances to launch
#COUNT=2 # predetermined number
COUNT=$(awk 'END{print NR}' $PARAMETERS) # for each input argument

# Security
keyNAME='CrownKey'
keyPATH="/home/ec2-user/.ssh/CrownKey.pem"
secGROUP='crown-group'

# Input Validation ======================
# Test input file for duplicate sample ID (which will cause file collisions downstream)
# exit with error 1 if collision exists

SAMPLEIDS=$(cut -f1 $1 | uniq -d )

if [ -n "$SAMPLEIDS" ]
then
  # If there are duplicate sample IDs exit
  echo "Error 02 - Duplicate Sample ID detected in input"
  echo " I hope you know what you're doing"
  echo " script will not exit in this version"
  echo " re-run with unique ID or outputs will overwrite"
  echo ""
  #echo "$SAMPLEIDS"
  #exit 2 
fi

# Script Core ===========================

for ITER in $(seq 1 $COUNT)
do

  # Count the total number of instances
  running_instance=$(aws ec2 describe-instances |\
    grep '"Name": "running"' - |\
    wc -l - | cut -f1 -d' ' - )

  while [ $running_instance -ge $MAX_INSTANCE ];
  do

  	echo "Reached max of $running_instance, wait."
  	echo ''

  	# If there are greater then 25 instances
  	# wait for instance to close
  	sleep 5m

  	running_instance=$(aws ec2 describe-instances |\
    grep '"Name": "running"' - |\
    wc -l - | cut -f1 -d' ' - )
  done

  # Extract Parameters/Arguments ----------

  ARGS=$(sed -n "$ITER"p $PARAMETERS | sed 's/\t/ /g' - )

  echo "Launch instance # $ITER"
  date
  echo "Instance Type: $instanceTYPE"
  echo "AMI Image: $imageID"
  echo "Run Script: $TASK"
  echo "Parameters: $ARGS"

  # Launch an instance --------------------
  # NOTE: each iteration of the for loop launches one instance
  # therefore each loop launches only one instance
  aws ec2 run-instances --image-id $imageID --count 1 \
   --instance-type $instanceTYPE --key-name $keyNAME \
   --block-device-mappings DeviceName=$devNAME,Ebs={VolumeSize=$volSIZE} \
   --security-groups $secGROUP > launch.tmp

  # Another alternative is to use --user-data droneB.sh 
  # which will run at instance boot-up
  # passing arguments to it may be challenging

  # Retrieve instance ID
  instanceID=$(cat launch.tmp | \
    egrep -o -e 'InstanceId[":/A-Za-z0-9_ \\-]*' - |\
    cut -f2 -d' ' - | xargs)

  echo "Instance ID: $instanceID"


  # Add a few minute wait here to allow for Public DNS to be assigned
  # otherwise ssh doesn't work
  # CURRENT BOTTLENECK
  sleep 180s

  # Check that the instance is 'running' and not pending
  pending_instance=$(aws ec2 describe-instances |\
   grep '"Name": "pending"' - |\
   wc -l - | cut -f1 -d' ' - )

  while [ $pending_instance -gt 0 ];
  do

  	echo "Launchign instance not ready yet, chill"
  	# Instance still pending
  	sleep 180s

    pending_instance=$(aws ec2 describe-instances |\
     grep '"Name": "pending"' - |\
     wc -l - | cut -f1 -d' ' - )
  done


  # Retrieve public DNS
  aws ec2 describe-instances --instance-ids $instanceID > launch2.tmp

  pubDNS=$(cat launch2.tmp | \
    egrep -o -m 1 -e 'PublicDnsName[.":/A-Za-z0-9_ \\-]*' - |\
    cut -f2 -d' ' - | xargs)

  echo "Public DNS: $pubDNS"

  # Access the instance -------------------

  LOGIN="ubuntu@$pubDNS" 

  ssh -i $keyPATH \
    -o StrictHostKeyChecking=no \
    $LOGIN 'bash -s' < droneB.sh $TASK $(echo $ARGS)

  # Cleanup
  rm *.tmp

  echo ''
  echo ''

done

# end of script