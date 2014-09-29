#!/bin/bash

#print first 3 lines from README file
head -n3 README.md 

echo "arguments: 1-2 - \"path\" \"flux.sng\", 3 - NuAnu, 4 - flavour, 5 - c_or_v, 6 - M_A, 7 - dMA, 8-10 - \"mixname\" Mode \"formula\", 11-12 - P_lep_min P_lep_max"
QESEMAT=./QESemat
fluxpath=$1
fluxname=$2
nu_arr=$3
flv_arr=$4
corv_arr=$5
m_a_arr=$6
dMA_arr=$7
mixname=$8
mode=$9
formula=${10}
p0=${11}
p1=${12}
echo $mixname "::::" $formula
#numbers to symbols:
nus="na"
flvs="emt"
corvs="cva"
dMAs="Ll0hH"

fluxfile="$fluxpath""$fluxname"
echo "arg1=" $fluxfile
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
					outfile=../output/"${fluxname%.*}"_"$mixname"_"$nu$flv"_"$corv"_"$m_a"_"$dMA".dat
					logfile=../output/"${fluxname%.*}"_"$mixname"_"$nu$flv"_"$corv"_"$m_a"_"$dMA".log
					echo $QESEMAT "$outfile" "$fluxfile" $nui $flvi $corvi $m_a $dMAi "$mixname" $mode "$formula" $p0 $p1 ">" $logfile
					$QESEMAT "$outfile" "$fluxfile" $nui $flvi $corvi $m_a $dMAi "$mixname" $mode "$formula" $p0 $p1 #> $logfile
				done
			done
		done
	done
done
