#!/bin/csh

if ($#argv == 1) root -l -b -q CompilePtRelAnalyzer.C

setenv WORKINGDIRECTORY $PWD

setenv TEMPLATEVARIABLE 'PtRel'
#setenv TEMPLATEVARIABLE 'System8' 

#setenv PUWEIGHTING ''
setenv PUWEIGHTING '_PSRun2018Prompt18'

#setenv KINWEIGHTING ''
#setenv KINWEIGHTING '_KinPtBinsCentral'
setenv KINWEIGHTING '_KinEtaAfterPtBinsCentral'
#setenv KINWEIGHTING '_KinPtInEtaBinsCentral'

#setenv SELECTION ''
#setenv SELECTION '_TrgEmul'
setenv SELECTION '_TrgConf'

#setenv MACRONAME 'ComputePileUpWeightsPSQCDMuQCDXJetHT'
#setenv MACRONAME 'ComputePileUpWeights_PSVRun2016BJuly15:JetHTQCDMuQCDX'
#setenv MACRONAME 'ComputePileUpWeightsTestPV'
setenv MACRONAME 'BuildTemplatesAll'
#setenv MACRONAME 'ComputeKinematicWeightsQCDMuJetHTQCDXAll'
#setenv MACRONAME 'CompareDataToMC_anyEta'
#setenv MACRONAME 'ComputePtRelScaleFactorsDeep:_Central:_LightTemplatesRatio'
#setenv MACRONAME 'ComputePtRelScaleFactorsDeep:_Central:_LightTemplatesRatio_bTempRatioCorr'
#setenv MACRONAME 'ComputePtRelScaleFactorsDeep:_Central:_LightTemplatesRatio_cJets_bTempRatioCorr'
#setenv MACRONAME 'ComputePtRelScaleFactorsDeep:_All:_LightTemplatesRatio_cJets_bTempRatioCorr'
#setenv MACRONAME 'PlotBTagPerformancePrompt18'
#setenv MACRONAME 'PlotBTagPerformanceKinEta'
#setenv MACRONAME 'AnalyzeSystematics_LightTemplatesRatio_cJets_bTempRatioCorr'
#setenv MACRONAME 'StoreScaleFactors_LightTemplatesRatio_cJets_bTempRatioCorr'
#setenv MACRONAME 'ProduceLightHistogramsQCD' # For testing

setenv DATARANGEINDEX '0'
echo ""

cd $WORKINGDIRECTORY/
echo "WORKINGDIRECTORY" $PWD
echo ""

eval `scramv1 runtime -csh`
#source /afs/cern.ch/project/eos/installation/cms/etc/setup.csh

echo "TEMPLATEVARIABLE " $TEMPLATEVARIABLE
echo "PUWEIGHTING      " $PUWEIGHTING
echo "KINWEIGHTING     " $KINWEIGHTING
echo "SELECTION        " $SELECTION
echo "MACRONAME        " $MACRONAME
echo "DATARANGEINDEX   " $DATARANGEINDEX
echo ""

mkdir -p Weights/KinematicWeights
mkdir -p Weights/PileUp
mkdir -p Tables
mkdir -p Plots

root -l -b -q 'RunPtRelAnalyzer.C'
 
   
