#!/bin/csh

if ($#argv == 1) root -l -b -q CompilePtRelAnalyzer.C

setenv WORKINGDIRECTORY $PWD

setenv TEMPLATEVARIABLE 'PtRel' 
#setenv TEMPLATEVARIABLE 'System8'  

setenv PUWEIGHTING '_PSRun20182018Ultimate'

setenv KINWEIGHTING ''
#setenv KINWEIGHTING '_KinPtBinsCentral'

#setenv SELECTION ''
#setenv SELECTION '_TrgEmul'
setenv SELECTION '_TrgConf'

#setenv DATARANGEINDEX '-999'

mkdir -p Weights/PileUp
mkdir -p Templates
mkdir -p Templates/Histograms

mkdir -p JobOutput
mkdir -p JobOutput/out
mkdir -p JobOutput/err
mkdir -p JobOutput/log

#setenv MACRONAME 'EventCounterQCDMuQCDX'
setenv MACRONAME 'EventCounterQCDX' 
setenv NJOBS 1

cp RunPtRelAnalyzerCondor.sub RunPtRelAnalyzerCondor_$MACRONAME.sub

sed -i 's|WORKINGDIRECTORY|'${WORKINGDIRECTORY}'|g' RunPtRelAnalyzerCondor_$MACRONAME.sub 
sed -i 's|TEMPLATEVARIABLE|'${TEMPLATEVARIABLE}'|g' RunPtRelAnalyzerCondor_$MACRONAME.sub  
sed -i 's|PUWEIGHTING|'${PUWEIGHTING}'|g'           RunPtRelAnalyzerCondor_$MACRONAME.sub  
sed -i 's|KINWEIGHTING|'${KINWEIGHTING}'|g'         RunPtRelAnalyzerCondor_$MACRONAME.sub  
sed -i 's|MACRONAME|'${MACRONAME}'|g'               RunPtRelAnalyzerCondor_$MACRONAME.sub  
sed -i 's|SELECTION|'${SELECTION}'|g'               RunPtRelAnalyzerCondor_$MACRONAME.sub  
sed -i 's|NJOBS|'${NJOBS}'|g'                       RunPtRelAnalyzerCondor_$MACRONAME.sub  

condor_submit RunPtRelAnalyzerCondor_$MACRONAME.sub

rm RunPtRelAnalyzerCondor_$MACRONAME.sub

