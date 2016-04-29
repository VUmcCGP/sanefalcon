INDIR=`readlink -f $1`
SUB=$2
for CHROM in $(seq 22 -1 1)
do
	echo $CHROM $SUB
	COLLECT=`ls $INDIR/*/merge.$CHROM | grep -v "/$SUB/merge"`
	#echo $COLLECT
	sort -n -m $COLLECT > $INDIR/$SUB/anti.$CHROM #" | qsub
done

exit

for SUB in `find $INDIR/* -maxdepth 0 -type d`
do
	for CHROM in $(seq 22 -1 1)
	do
		echo $CHROM $SUB
		#find $INDIR/* -maxdepth 0 -type d | grep -v $SUB
		COLLECT=`ls $INDIR/*/merge.$CHROM | grep -v "^$SUB"`
		#echo FOR:  $SUB/anti.$CHROM COLLECT: $COLLECT
		#echo $COLLECT
		sort -n -m $COLLECT > $SUB/anti.$CHROM #" | qsub
	done
done

