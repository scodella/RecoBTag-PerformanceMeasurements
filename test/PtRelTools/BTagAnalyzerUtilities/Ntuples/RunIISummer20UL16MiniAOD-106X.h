TString TreeContentFlag = "";

// ttbar MC
TString EOSPathttbar = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL16/MC/XXX/";
const int nttbarRanges = 1;
TString ttbarRange[nttbarRanges] = { "hadronic" };
int nttbarTrees[nttbarRanges] = { 0 };
double CrossSectionttbar[nttbarRanges] = { 313.9 };
double GeneratedEventsttbar[nttbarRanges] = { 0 };

/*
// QCDMu-nonAPV MC
TString EOSPathQCDMu = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL16/MC/QCD_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13/";
const int nMonteCarloPtHatRanges = 15;
TString MonteCarloPtHatRange[nMonteCarloPtHatRanges] = { "Pt-1000toInf", "Pt-120to170", "Pt-15to20", "Pt-170to300:1of2", "Pt-170to300:2of2", "Pt-20to30", "Pt-300to470:1of2", "Pt-300to470:2of2", "Pt-30to50", "Pt-470to600", "Pt-50to80", "Pt-600to800", "Pt-800to1000:1of2", "Pt-800to1000:2of2", "Pt-80to120" };
int nMonteCarloTrees[nMonteCarloPtHatRanges] = { 321, 365, 93, 655, 655, 449, 700, 700, 497, 439, 314, 468, 910, 910, 437 };
double CrossSection[nMonteCarloPtHatRanges] = { 10.4305*0.15544, 469797.*0.05362, 1.27319E9*0.003, 117989.*0.07335, 117989.*0.07335, 5.58528E8*0.0053, 7820.25*0.10196, 7820.25*0.10196, 1.39803E8*0.01182, 645.528*0.12242, 1.92225E7*0.02276, 187.109*0.13412, 32.3486*0.14552, 32.3486*0.14552, 2.758420E6*0.03844 };
double GeneratedEvents[nMonteCarloPtHatRanges] = { 14243204, 19830214, 3966021, 36800197, 36800197, 28515475, 30103401, 30103401, 35481419, 20185257, 22279218, 19553333, 41615643, 41615643, 23407998 };
*/

// QCDMu-nonAPV MC
TString EOSPathQCDMu = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL16/MC/QCD_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13_fixTrg/";
const int nMonteCarloPtHatRanges = 15;
TString MonteCarloPtHatRange[nMonteCarloPtHatRanges] = { "Pt-1000toInf", "Pt-120to170", "Pt-15to20", "Pt-170to300:1of2", "Pt-170to300:2of2", "Pt-20to30", "Pt-300to470:1of2", "Pt-300to470:2of2", "Pt-30to50", "Pt-470to600", "Pt-50to80", "Pt-600to800", "Pt-800to1000:1of2", "Pt-800to1000:2of2", "Pt-80to120" };
int nMonteCarloTrees[nMonteCarloPtHatRanges] = { 321, 365, 93, 655, 655, 449, 700, 700, 497, 439, 314, 468, 910, 910, 437 };
double CrossSection[nMonteCarloPtHatRanges] = { 10.4305*0.15544, 469797.*0.05362, 1.27319E9*0.003, 117989.*0.07335, 117989.*0.07335, 5.58528E8*0.0053, 7820.25*0.10196, 7820.25*0.10196, 1.39803E8*0.01182, 645.528*0.12242, 1.92225E7*0.02276, 187.109*0.13412, 32.3486*0.14552, 32.3486*0.14552, 2.758420E6*0.03844 };
double GeneratedEvents[nMonteCarloPtHatRanges] = { 14202400, 19830214, 3966021, 36800197, 36800197, 28515475, 30041015, 30041015, 35481419, 20185257, 22279218, 19488055, 41677403, 41677403, 23407998 };

// BTagMu Data
TString EOSPathBTagMu = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL16/Data/BTagMu_Run2016-21Feb2020_UL16_fixTrg/";
const int nBTagMuRanges = 3;
TString BTagMuRangeName[nBTagMuRanges] = { "Run2016F", "Run2016G", "Run2016H" };
int nBTagMuTrees[nBTagMuRanges] = { 25, 428, 475 };

// JetHT Data
TString EOSPathJetHT = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL16/Data/JetHT_Run2016-21Feb2020_UL16/";
const int nJetRunRanges = 3;
TString JetRunRangeName[nJetRunRanges] = { "Run2016F", "Run2016G", "Run2016H" };
int nJetTrees[nJetRunRanges] = { 25, 429, 476 };

// QCD MC
TString EOSPathQCD = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL16/MC/QCD_TuneCP5_13TeV_pythia8_RunIISummer19UL16MiniAOD-106X_mcRun2_asymptotic_v13/";
const int nMCInclusivePtHatRanges = 24;
TString MCInclusivePtHatRange[nMCInclusivePtHatRanges] = { "Pt-2400to3200", "Pt-30to50", "Pt-80to120", "Pt-170to300:1of2", "Pt-170to300:2of2", "Pt-50to80", "Pt-1400to1800", "Pt-600to800:1of4", "Pt-600to800:2of4", "Pt-600to800:3of4", "Pt-600to800:4of4", "Pt-1800to2400", "Pt-800to1000:1of2", "Pt-800to1000:2of2", "Pt-1000to1400", "Pt-470to600:1of2", "Pt-470to600:2of2", "Pt-300to470:1of3", "Pt-300to470:2of3", "Pt-300to470:3of3", "Pt-120to170:1of2", "Pt-120to170:2of2", "Pt-3200toInf", "Pt-15to30" };
int nMCInclusiveTrees[nMCInclusivePtHatRanges] = { 101, 359, 547, 710, 710, 367, 186, 1767, 1767, 1767, 1767, 124, 966, 966, 572, 711, 711, 1407, 1407, 1407, 659, 659, 61, 279 };
double CrossSectionInclusive[nMCInclusivePtHatRanges] = { 0.00682981, 1.40932E8, 2762530., 117276., 117276., 1.92043E7, 0.84265, 186.9, 186.9, 186.9, 186.9, 0.114943, 32.293, 32.293, 9.4183, 648.2, 648.2, 7823., 7823., 7823., 471100., 471100., 0.000165445, 1.83741E9 };
double GeneratedEventsInclusive[nMCInclusivePtHatRanges] = { 1998000, 19995697, 29817598, 29931300, 29931300, 18694200, 5912000, 67845999, 67845999, 67845999, 67845999, 3000000, 39988000, 39988000, 19389000, 27487000, 27487000, 57804000, 57804000, 57804000, 29831600, 29831600, 799000, 19785997 };

// https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/SimGeneral/MixingModule/python/mix_2016_25ns_UltraLegacy_PoissonOOTPU_cfi.py
float PileUpScenario[] = { 1.00402360149e-05, 5.76498797172e-05, 7.37891400294e-05, 0.000110932895295, 0.000158857714773,
                           0.000368637432599, 0.000893114107873, 0.00189700774575, 0.00358880167437, 0.00636052573486,
                           0.0104173961179, 0.0158122597405, 0.0223785660712, 0.0299186888073, 0.0380275944896,
                           0.0454313901624, 0.0511181088317, 0.0547434577348, 0.0567906239028, 0.0577145461461,
                           0.0578176902735, 0.0571251566494, 0.0555456541498, 0.053134383488, 0.0501519041462,
                           0.0466815838899, 0.0429244592524, 0.0389566776898, 0.0348507152776, 0.0307356862528,
                           0.0267712092206, 0.0229720184534, 0.0193388653099, 0.0159602510813, 0.0129310510552,
                           0.0102888654183, 0.00798782770975, 0.00606651703058, 0.00447820948367, 0.00321589786478,
                           0.0022450422045, 0.00151447388514, 0.000981183695515, 0.000609670479759, 0.000362193408119,
                           0.000211572646801, 0.000119152364744, 6.49133515399e-05, 3.57795801581e-05, 1.99043569043e-05,
                           1.13639319832e-05, 6.49624103579e-06, 3.96626216416e-06, 2.37910222874e-06, 1.50997403362e-06,
                           1.09816650247e-06, 7.31298519122e-07, 6.10398791529e-07, 3.74845774388e-07, 2.65177281359e-07,
                           2.01923536742e-07, 1.39347583555e-07, 8.32600052913e-08, 6.04932421298e-08, 6.52536630583e-08,
                           5.90574603808e-08, 2.29162474068e-08, 1.97294602668e-08, 1.7731096903e-08, 3.57547932012e-09,
                           1.35039815662e-09, 8.50071242076e-09, 5.0279187473e-09, 4.93736669066e-10, 8.13919708923e-10,
                           5.62778926097e-09, 5.15140589469e-10, 8.21676746568e-10, 0.0, 1.49166873577e-09,
                           8.43517992503e-09, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0 };


