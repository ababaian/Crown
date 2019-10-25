#!/bin/bash
# queenB.sh
#
# EC2 Launch / Control Script
#

# Control Panel =========================
# EC2 Run Script - script for droneB to execute
TASK="s3://crownproject/scripts/1kg_align_v0.sh"

# Parameter file, each line is given to a droneB to execute
# gather.sh by
PARAMETERS="1kg_runs_1.txt"

# EC2 Set-up
instanceTYPE='c4.2xlarge'
imageID='ami-66129306' #AMI

devNAME='/dev/sda1' # /dev/sda1 for Crown-AMI
volSIZE='200' # in Gb

# Number of instances to launch
#COUNT=2 # predetermined number
COUNT=$(wc -l $PARAMETERS | cut -f 1 -d' ' ) # for each input argument

# Security
keyNAME='XXXX'
keyPATH="XXXXXXX.pem"
secGROUP='XXXXXX'

# =======================================

for ITER in $(seq 1 $COUNT)
do

  # Extract Parameters/Arguments ----------

  ARGS=$(sed -n "$ITER"p $PARAMETERS | sed 's/\t/ /g' - )

  echo "Launch instance # $ITER"
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
  sleep 180s

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
