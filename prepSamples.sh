INDIR=$1
OUTDIR=$2
mkdir $OUTDIR

export DIR_FETALFRAC=~/fetalfrac
export SCRIPT_RETRO=$DIR_FETALFRAC/retro.py

export DIR_SCRIPTS=/illumina/diagnostics/bin
export SCRIPT_PYTHON=/illumina/diagnostics/bin/Python-2.7.3/bin/python
export SCRIPT_SAMTOOLS=$DIR_SCRIPTS/samtools-0.1.19/samtools

for SAMPLE in `find $INDIR -name "*.bam"`
do
	SHORT=`echo $SAMPLE | rev | cut -d"/" -f1 | rev`
	OUTFILE="$OUTDIR/${SHORT//".sort.bam"/}"
	echo "IN:$SAMPLE, OUT:$OUTFILE"
	
	for ARG_TASKID in `seq 1 22` # or "X"
		do
			$SCRIPT_SAMTOOLS view \
			$SAMPLE \
			chr$ARG_TASKID \
			-F 20 -q 1 | \
			$SCRIPT_PYTHON $SCRIPT_RETRO | \
					awk '{print $4}' \
					    > $OUTFILE.$ARG_TASKID.start.fwd &

			$SCRIPT_SAMTOOLS view \
			$SAMPLE \
			chr$ARG_TASKID \
			-f 16 -F 4 -q 1 | \
			$SCRIPT_PYTHON $SCRIPT_RETRO | \
					awk '{print $4+50}' \
					    > $OUTFILE.$ARG_TASKID.start.rev
		done
done
