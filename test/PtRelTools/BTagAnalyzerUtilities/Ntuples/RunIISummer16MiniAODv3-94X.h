// EOS paths 
TString EOSPathQCD       = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/2016Legacy/MC/QCD_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3/";
TString EOSPathQCDMu     = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/2016Legacy/MC/QCD_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3/";
TString EOSPathBTagMu    = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/2016Legacy/Data/BTagMu/";
TString EOSPathJetHT     = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/2016Legacy/Data/JetHT/";

TString TreeContentFlag = "";

// BTag data
const int nBTagMuRanges = 7;
TString BTagMuRangeName[nBTagMuRanges] = {"BTagMu_Run2016B-17Jul2018_ver2-v1", "BTagMu_Run2016C-17Jul2018-v1", "BTagMu_Run2016D-17Jul2018-v1", "BTagMu_Run2016E-17Jul2018-v1", "BTagMu_Run2016F-17Jul2018-v1", "BTagMu_Run2016G-17Jul2018-v1", "BTagMu_Run2016H-17Jul2018-v1"};
int nBTagMuTrees[nBTagMuRanges] = {523, 175, 292, 248, 181, 427, 475};

// Jet data
const int nJetRunRanges = 7;
TString JetRunRangeName[nJetRunRanges] = {"JetHT_Run2016B-17Jul2018_ver2-v2", "JetHT_Run2016C-17Jul2018-v1", "JetHT_Run2016D-17Jul2018-v1", "JetHT_Run2016E-17Jul2018-v1", "JetHT_Run2016F-17Jul2018-v1", "JetHT_Run2016G-17Jul2018-v1", "JetHT_Run2016H-17Jul2018-v1"};
int nJetTrees[nJetRunRanges]    = {523, 174, 292, 248, 181, 427, 475};

// QCD muon enriched 13 TeV
const int nMonteCarloPtHatRanges = 12;
TString MonteCarloPtHatRange[nMonteCarloPtHatRanges] = {    "Pt-15to20",      "Pt-20to30",       "Pt-30to50",       "Pt-50to80",       "Pt-80to120",   "Pt-120to170",   "Pt-170to300",             "Pt-300to470", "Pt-470to600",    "Pt-600to800",          "Pt-800to1000",     "Pt-1000toInf"};
double CrossSection[nMonteCarloPtHatRanges]          = {1.27319E9*0.003, 5.58528E8*0.0053, 1.39803E8*0.01182, 1.92225E7*0.02276, 2.758420E6*0.03844, 469797.*0.05362, 117989.*0.07335,         7820.25*0.10196,        645.528*0.12242,   187.109*0.13412,         32.3486*0.14552,       10.4305*0.15544};
double GeneratedEvents[nMonteCarloPtHatRanges]       = {4141251, 31878737, 29809491, 19662174, 13750712, 11912228, 19789671, 19095603, 9847662, 7003339, 5035680, 7472857};
int nMonteCarloTrees[nMonteCarloPtHatRanges]         = {31, 233, 220, 145, 103, 102, 157, 184, 76, 106, 87, 101};

// QCD inclusive 13 TeV 
const int nMCInclusivePtHatRanges = 15;                                       
TString MCInclusivePtHatRange[nMCInclusivePtHatRanges] = {"Pt-15to30", "Pt-30to50", "Pt-50to80",    "Pt-80to120", "Pt-120to170", "Pt-170to300",  "Pt-300to470",  "Pt-470to600", "Pt-600to800", "Pt-800to1000", "Pt-1000to1400", "Pt-1400to1800", "Pt-1800to2400", "Pt-2400to3200", "Pt-3200toInf"};
double CrossSectionInclusive[nMCInclusivePtHatRanges]  = {  1.83741E9,   1.40932E8,   1.92043E7,        2762530.,       471100.,       117276.,          7823.,         648.2,         186.9,         32.293,          9.4183,         0.84265,        0.114943,      0.00682981,    0.000165445};
double GeneratedEventsInclusive[nMCInclusivePtHatRanges] = {39184456, 9253499, 9954369, 6696872, 6867421, 6958707, 3860489, 3959986, 3896412, 3992110, 2999068, 0, 0, 0, 0};
int nMCInclusiveTrees[nMCInclusivePtHatRanges]           = {281, 71, 73, 57, 53, 56, 32, 31, 43, 36, 23, 0, 0, 0, 0};

// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py
float PileUpScenario[] = {1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05};
