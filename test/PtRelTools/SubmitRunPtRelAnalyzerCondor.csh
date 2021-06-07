#!/bin/csh

if ($#argv == 1) root -l -b -q CompilePtRelAnalyzer.C

setenv WORKINGDIRECTORY $PWD

setenv TEMPLATEVARIABLE 'PtRel' 
#setenv TEMPLATEVARIABLE 'System8'  

#setenv PUWEIGHTING 'None'
#setenv PUWEIGHTING '_PSRun20182018Ultimate'
#setenv PUWEIGHTING '_PSRun2018ABC2018Ultimate'
#setenv PUWEIGHTING '_PSRun2018D2018Ultimate'
#setenv PUWEIGHTING '_PSRun2017Moriond18'
#setenv PUWEIGHTING '_PSRun2017UL17'
#setenv PUWEIGHTING '_PSRun2018UL18'
setenv PUWEIGHTING '_PSRun2016UL16'

setenv KINWEIGHTING 'None'
#setenv KINWEIGHTING '_KinPtBinsCentral'
#setenv KINWEIGHTING '_KinEtaAfterPtBinsCentral'

#setenv SELECTION 'None'
#setenv SELECTION '_TrgEmul'
setenv SELECTION '_TrgConf'

mkdir -p Weights/PileUp
mkdir -p Templates
mkdir -p Templates/Histograms

mkdir -p JobOutput
mkdir -p JobOutput/out
mkdir -p JobOutput/err
mkdir -p JobOutput/log

#setenv MACRONAME 'EventCounterQCDMuQCDX'
#setenv MACRONAME 'EventCounterQCDMu' 
#setenv NJOBS     1

#setenv MACRONAME 'ComputeBTaggingWorkingPointsDeepCSV'
#setenv MACRONAME 'ComputeBTaggingWorkingPointsDeepJet'
#setenv NJOBS     1

setenv MACRONAME 'ProduceHistogramsBTagMu'
setenv NJOBS     5 
#setenv NJOBS     3 

#setenv MACRONAME 'ProduceHistogramsQCDMu'
#setenv NJOBS     15

#setenv MACRONAME 'ProduceLightHistogramsJetHT'
#setenv NJOBS     5
#setenv NJOBS     3

#setenv MACRONAME 'ProduceLightHistogramsQCD'
#setenv NJOBS     22
#setenv NJOBS     24

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

