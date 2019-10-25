#!/bin/bash
# Converts a gapped ALIGN and its CONSERVATION to a wig file
ALIGN='lsu_vert_align.fa' # Gapped Fasta Alignment
CONS='lsu_vert_cons.csv' # JalView Consensus Output
OUTPUT='VertCons_LSU.wig'

# Start position of SSU (1003657) or 5.8S (1006623) or LSU (1007935)

# Convert ALIGN.fa into one nucleotide per row

sed 1d $ALIGN |
sed 's/\./\-\n/g' - |
sed 's/[ATGCU]/\N\n/g' - |
sed '/^\s*$/d' - > align.csv


sed 2d $CONS |
sed 's/,/\n/g' - |
sed 1d - |
head -n -1 - > cons.csv
 
paste align.csv cons.csv |
grep -e '^N' - |
cut -f 2 - |
sed 's/ //g' - > cons.values


# SSU Alignments

	#echo "fixedStep chrom=chr13 start=1003657 step=1" > $OUTPUT
	#cat cons.values >> $OUTPUT
	# rm align.csv cons.csv cons.values 

# LSU Alignments
# there are differences between hgr and riboZone LSU

	# Delete in hgr (28S) relative to riboZone
	# 40 nucleotides removed
	sed 4978,4980d cons.values |
	sed 4963d - |
	sed 4879d - |
	sed 4297,4303d - |
	sed 4260d - |
	sed 4182d - |
	sed 3730d - |
	sed 3539d - |
	sed 3348,3349d - |
	sed 3332,3335d - |
	sed 3319d - |
	sed 3274d - |
	sed 3137d - |
	sed 3119,3120d - |
	sed 939,942d - |
	sed 417,419d - |
	sed 407,409d - |
	sed 293,295d - > cons.d.values

	# Insert in hgr (28S) relative to riboZone
	sed 3472i0.1 cons.d.values |
	sed 3454i0.1 - |
	sed 3441i0.1 - |
	sed 3327i0.1 - |
	sed 3269i0.1 - > cons.di.values

	# 5.8S subunit (1-157 nt)
	sed -i '158ifixedStep chrom=chr13 start=1007935 step=1' cons.di.values
	sed -i '1ifixedStep chrom=chr13 start=1006623 step=1' cons.di.values
	mv cons.di.values $OUTPUT

	rm align.csv cons.csv cons.values cons.d.values


# NOTE: There seems to be about a 2bp mis-alignment with the above script
# but it's kind of hit and miss. Don't take the bp value directly (or refer
# back to riboZone) but the general consensus trend does hold pretty well.
