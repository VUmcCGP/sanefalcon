SCRIPT_PYTHON=python
SCRIPT_PREDICT=predict.py

FILE_CORR=$1
FILE_OUTPUT=$2

LC_NUMERIC=C \
	awk -F"," \
	'{ for(i=1;i<=NF;i++){ fetal[i]+=$i } } END { for(i=1;i<=NF;i++) {printf "%s,",fetal[i]}; print "";}' \
	$FILE_OUTPUT.*.rev.np $FILE_OUTPUT.*.fwd.np \
		> $FILE_OUTPUT.cnp


LC_NUMERIC=C \
	awk -F"," \
	'{ for(i=1;i<=NF;i++){ fetal[i]+=$i } } END { for(i=1;i<=NF;i++) {printf "%s,",fetal[i]}; print "";}' \
	$FILE_OUTPUT.*.irev.np $FILE_OUTPUT.*.ifwd.np \
		>> $FILE_OUTPUT.cnp

$SCRIPT_PYTHON $SCRIPT_PREDICT $FILE_CORR $FILE_OUTPUT.cnp > $FILE_OUTPUT.ff
