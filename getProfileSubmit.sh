CHROM=${SGE_TASK_ID}

$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex5.$CHROM $SIMPLE.$CHROM.start.fwd 0 $OUTDIR/$SIMPLER.$CHROM.fwd 
$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex5.$CHROM $SIMPLE.$CHROM.start.rev 1 $OUTDIR/$SIMPLER.$CHROM.rev


$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex5.$CHROM $SIMPLE.$CHROM.start.fwd 1 $OUTDIR/$SIMPLER.$CHROM.ifwd 
$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex5.$CHROM $SIMPLE.$CHROM.start.rev 0 $OUTDIR/$SIMPLER.$CHROM.irev

