#!/bin/csh

if ($#argv == 1) root -l -b -q CompilePtRelAnalyzer.C

setenv WORKINGDIRECTORY $PWD

setenv TEMPLATEVARIABLE 'PtRel'
#setenv TEMPLATEVARIABLE 'System8' 

#setenv PUWEIGHTING ''
#setenv PUWEIGHTING '_PSRun20182018Ultimate'
#setenv PUWEIGHTING '_PSRun2018ABC2018Ultimate'
#setenv PUWEIGHTING '_PSRun2018D2018Ultimate'
#setenv PUWEIGHTING '_PSRun2017Moriond18'
#setenv PUWEIGHTING '_PSRun2017UL17'
#setenv PUWEIGHTING '_PSRun2018UL18'
setenv PUWEIGHTING '_PSRun2016UL16APV'
#setenv PUWEIGHTING '_PSRun2016UL16'

setenv KINWEIGHTING ''
#setenv KINWEIGHTING '_KinPtBinsCentral'
#setenv KINWEIGHTING '_KinEtaAfterPtBinsCentral'
#setenv KINWEIGHTING '_KinPtInEtaBinsCentral'

#setenv SELECTION ''
#setenv SELECTION '_TrgEmul'
setenv SELECTION '_TrgConf'

setenv MACRONAME 'ComputePileUpWeightsPSQCDMuQCDXJetHT'
#setenv MACRONAME 'ComputePileUpWeights_PSVRun2016BJuly15:JetHTQCDMuQCDX'
#setenv MACRONAME 'ComputePileUpWeightsTestPV'
#setenv MACRONAME 'BuildTemplates'
#setenv MACRONAME 'BuildTemplatesAll'
#setenv MACRONAME 'ComputeKinematicWeightsQCDMuAll'
#setenv MACRONAME 'ComputeKinematicWeightsQCDMuJetHTQCDXAll'
#setenv MACRONAME 'CompareDataToMC_anyEta'
#setenv MACRONAME 'ComputePtRelScaleFactorsDeep:_Central:_LightTemplatesRatio_bTempRatioCorr'
#setenv MACRONAME 'ComputePtRelScaleFactorsDeep:_Central:_LightTemplatesRatio_cJets_bTempRatioCorr'
#setenv MACRONAME 'ComputePtRelScaleFactorsDeep:_All:_LightTemplatesRatio_cJets_bTempRatioCorr'
#setenv MACRONAME 'PlotBTagPerformanceNonoalgo'
#setenv MACRONAME 'PlotBTagPerformanceEfficiency2018'
#setenv MACRONAME 'AnalyzeSystematics_LightTemplatesRatio_cJets_bTempRatioCorr'
#setenv MACRONAME 'StoreScaleFactors_LightTemplatesRatio_cJets_bTempRatioCorr'
#setenv MACRONAME 'ProduceHistogramsQCDMu' # For testing
#setenv MACRONAME 'ProduceLightHistogramsJetHT' # For testing
#setenv MACRONAME 'ComputeBTaggingWorkingPointsDeepCSV'

setenv DATARANGEINDEX '0'
echo ""

cd $WORKINGDIRECTORY/
echo "WORKINGDIRECTORY" $PWD
echo ""

#eval `scramv1 runtime -csh`
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
 
   
