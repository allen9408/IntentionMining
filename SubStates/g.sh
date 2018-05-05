#!/bin/bash
outputFile=../result.dot

function fsub()
{
	# $1 = state index	
	# $2 = depth
	local curState=$1
	local curDepth=$2
	local nextDepth=$[curDepth+1]
	# echo depth: $curDepth
	# echo $curState
	# echo $nextDepth
	local stateNum=0
	for j in `ls images/*.png`;do
		local num=$(echo $j | grep -o '_' | wc -l)
		# echo $num
		if [ $num -eq $nextDepth ];then
			stateNum=$((stateNum+1))
			# echo -e $j", \c"
			local graphName=${j/images\//}
			graphName=${graphName/SubStates_/}
			graphName=${graphName/.png/}
			local checkState=$curState"_"
			# echo $'\t'$graphName[shapefile=\"SubStates\/$j\"]';' >> $outputFile
			if [[ ${graphName/${checkState}//} != $graphName ]];then
				for (( k=1;k<=$curDepth;k++));do
					echo -e $'\t''\c' >> $outputFile
				done
				# check whether it's the deepest level
				local nextCheck=$graphName"_"
				local haveDeeper=false;
				for k in `ls images/*.png`;do
					if [[ ${k/${nextCheck}//} != $k ]];then
						# echo ${k/${nextCheck}//}
						haveDeeper=true;
						break
					fi
				done
				# echo $haveDeeper
				if [[ $haveDeeper = true ]]; then
					echo $curState -\> $graphName';' >> $outputFile
				else
					echo $curState -\> $graphName[style = bold, color = red]';' >> $outputFile

				fi
				fsub $graphName $nextDepth;
			fi
		fi
	done
	# echo depth: $curDepth have $stateNum states
}


maxDepth=1;
logfile=~/Desktop/lambda_$1.txt
echo digraph g > $outputFile
echo \{ >> $outputFile
for i in `ls *.dot`;do
	# echo $i
	outfile=${i/dot/png}
	dot -T png $i -o images/$outfile
	# get depth of the state
	num=$(echo $outfile | grep -o '_' | wc -l)
	# echo $num
	if [ $num -gt $maxDepth ];then
		# echo $num
		maxDepth=$[num]
	fi
	# echo $num
done
echo Max depth: $maxDepth >> $logfile
for (( i=1;i<=$maxDepth;i++ ));do
	stateNum=0
	for j in `ls images/*.png`;do	
		num=$(echo $j | grep -o '_' | wc -l)
		# echo $num
		if [ $num -eq $i ];then
			stateNum=$((stateNum+1))
			# echo -e $j", \c"
			graphName=${j/images\//}
			graphName=${graphName/SubStates_/}
			graphName=${graphName/.png/}
			echo $'\t'$graphName[shapefile=\"SubStates\/$j\", label = \"\"]';' >> $outputFile
		fi
	done
	echo depth: $i has $stateNum states >> $logfile
done
depth=1
index=0
for i in `ls images/*.png`;do
	
	num=$(echo $i | grep -o '_' | wc -l)
	if [ $num -eq $depth ];then
		graphName=${i/images\//}
		graphName=${graphName/SubStates_/}
		graphName=${graphName/.png/}
		array[$index]=$graphName
		# index=`expr $index+1`
		index=$((index+1))
		# echo $index
	fi
done
for state in ${array[@]};do
	fsub $state $depth;
done


echo } >> $outputFile
cd ..
dot -T pdf result.dot -o ../result.pdf
