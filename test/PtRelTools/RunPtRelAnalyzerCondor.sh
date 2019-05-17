#!/bin/sh
#
echo ""
cd $1/
echo "WORKINGDIRECTORY" $PWD
eval `scramv1 runtime -sh`
#
export TEMPLATEVARIABLE=$2
export PUWEIGHTING=$3
export KINWEIGHTING=$4
export MACRONAME=$5
export SELECTION=$6
export DATARANGEINDEX=$(($7+1))
#
if [ "$PUWEIGHTING" == "None" ]; then
    export PUWEIGHTING=""
fi
if [ "$KINWEIGHTING" == "None" ]; then
    export KINWEIGHTING=""
fi
if [ "$SELECTION" == "None" ]; then
    export SELECTION=""
fi
#
echo ""
echo "TEMPLATEVARIABLE " $TEMPLATEVARIABLE
echo "PUWEIGHTING      " $PUWEIGHTING
echo "KINWEIGHTING     " $KINWEIGHTING
echo "MACRONAME        " $MACRONAME
echo "SELECTION        " $SELECTION
echo "DATARANGEINDEX   " $DATARANGEINDEX
#
echo ""
root -l -b -q RunPtRelAnalyzer.C
rm core.* 
echo ""
   
