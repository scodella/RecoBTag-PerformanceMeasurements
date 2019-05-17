#!/bin/sh
echo ""
cd $1/
echo "WORKINGDIRECTORY" $PWD
eval `scramv1 runtime -sh`
echo ""
export TEMPLATEVARIABLE=$2
export PUWEIGHTING=$3
export KINWEIGHTING=$4
export MACRONAME=$5
export SELECTION=$6
export DATARANGEINDEX=$(($7+1))
echo ""
echo "TEMPLATEVARIABLE " $TEMPLATEVARIABLE
echo "PUWEIGHTING      " $PUWEIGHTING
echo "KINWEIGHTING     " $KINWEIGHTING
echo "SELECTION        " $SELECTION
echo "MACRONAME        " $MACRONAME
echo "DATARANGEINDEX   " $DATARANGEINDEX
echo ""
root -l -b -q 'RunPtRelAnalyzer.C'
rm core.* 
echo ""
   
