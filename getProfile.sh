export SCRIPT_GETPRO=~/Workspace/sanefalcon/getProfile.py
export SCRIPT_GETPROSUB=~/Workspace/sanefalcon/getProfileSubmit.sh

export SCRIPT_PYTHON=python

export INDIR=`readlink -f $1`
export NUCLDIR=`readlink -f $2`
export OUTDIR=`readlink -f $3`

mkdir $OUTDIR

for SAMPLE in `find $INDIR -name "*.1.start.fwd"`
do
	export SIMPLE=${SAMPLE//.1.start.fwd/}
	export SIMPLER=`echo $SIMPLE | rev | cut -d"/" -f1 | rev`
	qsub -t 1-22:1 -V $SCRIPT_GETPROSUB
done
