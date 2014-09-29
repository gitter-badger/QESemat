#!/bin/bash

#print first 3 lines from README file
head -n3 README.md 

echo "arguments: 1 - NuAnu, 2 - flavour, 3 - c_or_v, 4 - M_A,  5-7 - \"mixname\" Mode \"formula\", 8-9 - E_nu_min E_nu_max"
PROGRAM=./Test_Section
nu_arr=$1
flv_arr=$2
corv_arr=$3
m_a_arr=$4
mixname=$5
mode=$6
formula=$7
e0=$8
e1=$9
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
				outfile=../output/"$mixname"_"$nu$flv"_"$corv"_$m_a.dat
				logfile=../output/"$mixname"_"$nu$flv"_"$corv"_$m_a.log
				echo $PROGRAM "$outfile" $nui $flvi $corvi $m_a "$mixname" $mode "$formula" $e0 $e1 ">" $logfile
				$PROGRAM "$outfile" $nui $flvi $corvi $m_a "$mixname" $mode "$formula" $e0 $e1 #> $logfile
			done
		done
	done
done
