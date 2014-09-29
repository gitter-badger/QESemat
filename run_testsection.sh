#!/bin/bash

#print first 3 lines from README file
head -n3 README.md 

echo "arguments: 1 - NuAnu, 2 - flavour, 3 - c_or_v, 4 - M_A, 7 - dMA, 6-8 - \"mixname\" Mode \"formula\", 9-10 - E_nu_min E_nu_max"
PROGRAM=./Test_Section
nu_arr=$1
flv_arr=$2
corv_arr=$3
m_a_arr=$4
dMA_arr=$5
mixname=$6
mode=$7
formula=$8
e0=$9
e1=${10}
echo $mixname "::::" $formula
#numbers to symbols:
nus="na"
flvs="emt"
corvs="cva"
dMAs="Ll0hH"

IFS=,
for nu in $nu_arr; do
	for flv in $flv_arr; do
	  	for corv in $corv_arr; do
	  		for m_a in $m_a_arr; do
	  			for dMA in $dMA_arr; do
					#find symbol position of $nu in $nus
					nui=`expr index $nus $nu`
					#find symbol position of $nu in $nus
					flvi=`expr index $flvs $flv`
					#find symbol position of $flv in $flvs
					corvi=`expr index $corvs $corv`
					dMAi=`expr index $dMAs $dMA`
					outfile=../output/"$mixname"_"$nu$flv"_"$corv"_"$m_a"_"$dMA".dat
					logfile=../output/"$mixname"_"$nu$flv"_"$corv"_"$m_a"_"$dMA".log
					echo $PROGRAM "$outfile" $nui $flvi $corvi $m_a $dMAi "$mixname" $mode "$formula" $e0 $e1 ">" $logfile
					$PROGRAM "$outfile" $nui $flvi $corvi $m_a $dMAi "$mixname" $mode "$formula" $e0 $e1 #> $logfile
				done
			done
		done
	done
done
