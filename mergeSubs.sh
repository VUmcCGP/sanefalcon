INDIR=`readlink -f $1`

for i in $(seq 1 22)
do
	echo Working on Chr: $i
	echo "sort -n -m $INDIR/*/merge.$i > $INDIR/merge.$i" | qsub
done
