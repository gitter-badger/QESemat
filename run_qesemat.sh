#!/bin/bash

#print first 3 lines from README file
head -n3 README.md 

echo "arguments: 1 - flux.sng, 2 - NuAnu, 3 - flavour, 4 - c_or_v, 5 - M_A,  6-7 - \"mixture\" \"formula\", 8 - P_lep_min, 9 - P_lep_max"
QESEMAT=./Test_Section
fluxfile=$1
nu_arr=$2
flv_arr=$3
corv_arr=$4
m_a_arr=$5
mixname=$6
formula=$7
p0=$8
p1=$9
echo $mixname "::::" $formula
#numbers to symbols:
nus="na"
flvs="emt"
corvs="cv"

echo "arg1=" $fluxfile
IFS=,
for nu in $nu_arr; do
  for flv in $flv_arr; do
	  	for corv in $corv_arr; do
	  		for m_a in $m_a_arr; do
				#find symbol position of $nu in $nus
				nui=`expr index $nus $nu`
				#find symbol position of $nu in $nus
				flvi=`expr index $flvs $flv`
				#find symbol position of $flv in $flvs
				corvi=`expr index $corvs $corv`
				outfile="${fluxfile%.*}"_"$mixname"_"$nu$flv"_"$corv"_$m_a.dat
				logfile="${fluxfile%.*}"_"$mixname"_"$nu$flv"_"$corv"_$m_a.log
				echo $QESEMAT "$outfile" "$fluxfile" $nui $flvi $corvi $m_a "$mixname" "$formula" $p0 $p1 ">" $logfile
				$QESEMAT "$outfile" "$fluxfile" $nui $flvi $corvi $m_a "$mixname" "$formula" $p0 $p1 #> $logfile
			done
		done
	done
done
