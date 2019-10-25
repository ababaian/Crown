#!/bin/bash
# Modified - AB
# Usage: <input.bam> <AWS_ID>
# 
# dependencies, jq and parallel:
# sudo dnf install jq
# sudo dnf install parallel
# sudo pip install awscli

# Input File to upload
INFILE=$1
AWSID=$2

# Control Panel ==================================
#byteSize=1048576 # 1 Mb chunks
byteSize=4194304 # 4 Mb chunks
# ================================================


# Split the input file and prefix with part
split --bytes=$byteSize --verbose $INFILE part

# count the number of files that begin with "part"
fileCount=$(ls -1 | grep "^part" | wc -l)
echo "Total parts to upload: " $fileCount

# get the list of part files to upload.  Edit this if you chose a different prefix in the split command
files=$(ls | grep "^part")

echo "--------------------------------------"
# initiate multipart upload connection to glacier
init=$(aws glacier initiate-multipart-upload --account-id $AWSID --part-size $byteSize --vault-name crown --archive-description "$INFILE")
echo $init
echo ''


location=$(echo $init | egrep -o -e 'location[":/A-Za-z0-9_ \\-]*' - | cut -f2 -d' ' - | xargs)
uploadId=$(echo $init | egrep -o -e 'uploadId[":/A-Za-z0-9_ \\-]*' - | cut -f2 -d' ' - | xargs)

echo "Location: "$location
echo "UploadID: "$uploadId

echo "---------------------------------------"
# xargs trims off the quotes
# jq pulls out the json element titled uploadId
# uploadId=$(echo $init | jq '.uploadId' | xargs)


# create upload commands to be run in parallel and store in commands.txt
i=0
for f in $files 
  do
     fileSize=$(ls -l $f | cut -f5 -d' ' -) # file size
    
     byteStart=$((i*byteSize))
     byteEnd=$((i*byteSize+fileSize-1))
     byteRange=$(echo "'bytes $byteStart-$byteEnd/*'")

     echo "aws glacier upload-multipart-part --body $f --range $byteRange  --account-id $AWSID --vault-name crown --upload-id $uploadId" > cmd.txt
     sh cmd.txt ; rm cmd.txt

     i=$(($i+1))
     
  done

# treeHash Variable
TREEHASH=$(./treehash.py $INFILE)

fileSize=$(ls -l $INFILE | cut -f5 -d' ' - )

echo "List Active Multipart Uploads:"
echo "Verify that a connection is open:"
aws glacier list-multipart-uploads --account-id $AWSID --vault-name crown

# end the multipart upload
#aws glacier abort-multipart-upload --account-id $AWSID --vault-name crown --upload-id $uploadId

rm part*

COMPLETE=$(aws glacier complete-multipart-upload --checksum $TREEHASH --archive-size $fileSize --upload-id $uploadId --account-id $AWSID --vault-name crown)

echo $COMPLETE > $INFILE.glacier
