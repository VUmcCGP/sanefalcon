export SCRIPT_REFFF=~/sanefalcon/getRefFF.py

export SCRIPT_PYTHON=/illumina/diagnostics/bin/Python-2.7.3/bin/python

for SAMPLE in `find $1 -name "*.22.start.fwd"`
do
	#NAME=`echo $SERIES | rev | cut -d"/" -f1 | rev`
	#echo ${SAMPLE//.22.start.fwd/}
	$SCRIPT_PYTHON $SCRIPT_REFFF ${SAMPLE//.22.start.fwd/}
done

