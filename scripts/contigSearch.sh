#!/bin/bash
# contigSearch.sh
#

INPUT='NA12_contigs.fa'

# Initialize local BLAT server 
	#gfServer start <host> <port> <genome>
	gfServer start glitch 2666 hg19.2bit &
	sleep 60s


Line=$(wc -l $INPUT | cut -f1 -d' ')
Count='1'

# Fasta Loop ==============================================
while [ $Count -lt $Line ]
do

# Count and Count+1 to extract
Count1=$((Count + 1))

# Extract the input
sed -n "$Count,$Count1"p $INPUT > iteration.fa

# gfClient <host> <port> <seqDir> <input.fa> <output.psl>
gfClient glitch 2666 . iteration.fa tmp.psl

# Check if a match exists
	# If tmp.psl has 5 lines => no matches
	# If tmp.psl has <5 lines => match exists

if [ '5' -lt $(wc -l tmp.psl | cut -f1 -d' ') ]
then
	# Matches found in Genome
	# Just continue to the next iteration
	echo $Count
else
	# No matches found in Genome
	#
	echo $Count
	cat iteration.fa >> notInGenome.mfa
	
fi

# Count up to the next odd
Count=$((Count + 2))

done
# =========================================================

# Stop BLAT server
	#gfServer stop <host> <port>
	#gfServer stop glitch 2666
	kill %1
