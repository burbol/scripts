STARTTIME=$(date +%s)
#command block that takes time to complete...
#........
ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."