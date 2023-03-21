
listBranchNo=$1
echo $listBranchNo

for (( i=1; i<=7; i++ ))
do
    echo $i
    DIR=$(sed -n ${i}p $listBranchNo)
    echo $DIR
    # for (( j=0; j<4; j++ ))
    # do
    #     folder1="${inputFiles[$i]}/spikeFileLong.csv"
    #     folder2="${inputFiles[$i]}/correlation.csv"
    #     matlab -nodisplay -r "initialize_network('$DIR')"
    # done

done
