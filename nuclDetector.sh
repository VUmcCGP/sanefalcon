export SCRIPT_NUCDEC=~/sanefalcon/nuclDetector.py

export SCRIPT_PYTHON=/illumina/diagnostics/bin/Python-2.7.3/bin/python

export INDIR=`readlink -f $1`

for CHROM in $(seq 22 -1 1)
do
	echo "$SCRIPT_PYTHON $SCRIPT_NUCDEC $INDIR/merge.$CHROM $INDIR/nucl_ex4.$CHROM" | qsub
done

