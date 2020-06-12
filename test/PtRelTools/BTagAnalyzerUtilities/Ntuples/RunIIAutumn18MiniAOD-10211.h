// EOS paths 
TString EOSPathQCD       = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/2018_Ultimate/MC/QCD_TuneCP5_13TeV_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15/";
TString EOSPathQCDMu     = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/2018_Ultimate/MC/QCD_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15/";
TString EOSPathBTagMu    = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/2018_Ultimate/Data/BTagMu/";
TString EOSPathJetHT     = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/2018_Ultimate/Data/JetHT/";

TString TreeContentFlag = "";

// BTag data
const int nBTagMuRanges = 6;
TString BTagMuRangeName[nBTagMuRanges] = {"BTagMu_Run2018A_17Sep2018-v1", "BTagMu_Run2018B_17Sep2018-v1", "BTagMu_Run2018C_17Sep2018-v1", "BTagMu_Run2018D-PromptReco-v2:1of3", "BTagMu_Run2018D-PromptReco-v2:2of3", "BTagMu_Run2018D-PromptReco-v2:3of3"};
int nBTagMuTrees[nBTagMuRanges] = {554, 262, 262, 1255, 1255, 1255};

// Jet data
const int nJetRunRanges = 6;
TString JetRunRangeName[nJetRunRanges] = {"JetHT_Run2018A_17Sep2018-v1", "JetHT_Run2018B_17Sep2018-v1", "JetHT_Run2018C_17Sep2018-v1", "JetHT_Run2018D-PromptReco-v2:1of3", "JetHT_Run2018D-PromptReco-v2:2of3", "JetHT_Run2018D-PromptReco-v2:3of3"};
int nJetTrees[nJetRunRanges]    = {554, 264, 262, 1254, 1254, 1254};

// QCD muon enriched 13 TeV
const int nMonteCarloPtHatRanges = 12;
TString MonteCarloPtHatRange[nMonteCarloPtHatRanges] = {    "Pt-15to20",      "Pt-20to30",       "Pt-30to50",       "Pt-50to80",       "Pt-80to120",   "Pt-120to170",   "Pt-170to300",             "Pt-300to470", "Pt-470to600",    "Pt-600to800",          "Pt-800to1000",     "Pt-1000toInf"};
double CrossSection[nMonteCarloPtHatRanges]          = {1.27319E9*0.003, 5.58528E8*0.0053, 1.39803E8*0.01182, 1.92225E7*0.02276, 2.758420E6*0.03844, 469797.*0.05362, 117989.*0.07335,         7820.25*0.10196,        645.528*0.12242,   187.109*0.13412,         32.3486*0.14552,       10.4305*0.15544};
double GeneratedEvents[nMonteCarloPtHatRanges]       = {4576065, 30612338, 29076487, 19921151, 25039361, 19116008, 35025896, 27158639, 19869770, 16363909, 14974995, 10410480};
int nMonteCarloTrees[nMonteCarloPtHatRanges]         = {129, 454, 509, 338, 471, 387, 647, 734, 544, 419, 466, 279};

// QCD inclusive 13 TeV 
const int nMCInclusivePtHatRanges = 17;                                       
TString MCInclusivePtHatRange[nMCInclusivePtHatRanges] = {"Pt-15to30", "Pt-30to50", "Pt-50to80",    "Pt-80to120", "Pt-120to170", "Pt-170to300",  "Pt-300to470",  "Pt-470to600", "Pt-600to800:1of3", "Pt-600to800:2of3", "Pt-600to800:3of3", "Pt-800to1000", "Pt-1000to1400", "Pt-1400to1800", "Pt-1800to2400", "Pt-2400to3200", "Pt-3200toInf"};
double CrossSectionInclusive[nMCInclusivePtHatRanges]  = {  1.83741E9,   1.40932E8,   1.92043E7,        2762530.,       471100.,       117276.,          7823.,         648.2,         186.9,         186.9,         186.9,         32.293,          9.4183,         0.84265,        0.114943,      0.00682981,    0.000165445};
double GeneratedEventsInclusive[nMCInclusivePtHatRanges] = {19451000, 18872000, 12909000, 29535000, 25255000, 29710000, 41744000, 17712000, 64061000, 64061000, 64061000, 9436000, 18485000,  0, 0, 0, 0};
int nMCInclusiveTrees[nMCInclusivePtHatRanges]           = {277, 265, 163, 402, 332, 478, 588, 247, 1573, 1573, 1573, 220, 458, 0, 0, 0, 0};

// https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3271/1/1/1.html
// https://github.com/cms-sw/cmssw/blob/CMSSW_10_4_X/SimGeneral/MixingModule/python/mix_2018_25ns_JuneProjectionFull18_PoissonOOTPU_cfi.py

float PileUpScenario[] = {    4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05,
			      3.508135e-05, 7.12608e-05, 0.0001400641, 0.0002663403, 0.0004867473,
			      0.0008469, 0.001394142, 0.002169081, 0.003198514, 0.004491138,
			      0.006036423, 0.007806509, 0.00976048, 0.0118498, 0.01402411,
			      0.01623639, 0.01844593, 0.02061956, 0.02273221, 0.02476554,
			      0.02670494, 0.02853662, 0.03024538, 0.03181323, 0.03321895,
			      0.03443884, 0.035448, 0.03622242, 0.03674106, 0.0369877,
			      0.03695224, 0.03663157, 0.03602986, 0.03515857, 0.03403612,
			      0.0326868, 0.03113936, 0.02942582, 0.02757999, 0.02563551,
			      0.02362497, 0.02158003, 0.01953143, 0.01750863, 0.01553934,
			      0.01364905, 0.01186035, 0.01019246, 0.008660705, 0.007275915,
			      0.006043917, 0.004965276, 0.004035611, 0.003246373, 0.002585932,
			      0.002040746, 0.001596402, 0.001238498, 0.0009533139, 0.0007282885,
			      0.000552306, 0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012,
			      0.0001238161, 8.96531e-05, 6.438087e-05, 4.585302e-05, 3.23949e-05,
			      2.271048e-05, 1.580622e-05, 1.09286e-05, 7.512748e-06, 5.140304e-06,
			      3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07,
			      5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07,
			      1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08,
			      6.119019e-08, 5.443693e-08, 4.85036e-08, 4.31486e-08, 3.822112e-08   };
