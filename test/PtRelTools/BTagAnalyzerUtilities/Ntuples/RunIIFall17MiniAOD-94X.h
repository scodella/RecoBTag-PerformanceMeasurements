// EOS paths 
TString EOSPathQCD       = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/CMSSW_9_4_1/MC/QCD_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10/";
TString EOSPathQCDMu     = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/CMSSW_9_4_1/MC/QCD_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10/";
TString EOSPathBTagMu    = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/CMSSW_9_4_1/Data/BTagMu/";
TString EOSPathJetHT     = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/CMSSW_9_4_1/Data/JetHT/";

TString TreeContentFlag = "_FatJets_Subjets";

// BTag data
const int nBTagMuRanges = 7;
TString BTagMuRangeName[nBTagMuRanges] = {"BTagMu_Run2017B_17Nov2017-v1", "BTagMu_Run2017C_17Nov2017-v1:1of2", "BTagMu_Run2017C_17Nov2017-v1:2of2", "BTagMu_Run2017D_17Nov2017-v1", "BTagMu_Run2017E_17Nov2017-v1", "BTagMu_Run2017F_17Nov2017-v1:1of2", "BTagMu_Run2017F_17Nov2017-v1:2of2"};
int nBTagMuTrees[nBTagMuRanges] = {244, 525, 525, 274, 285, 464, 464};

// Jet data
const int nJetRunRanges = 15;
TString JetRunRangeName[nJetRunRanges] = {"JetHT_Run2017B_17Nov2017-v1:1of2", "JetHT_Run2017B_17Nov2017-v1:2of2", 
					  "JetHT_Run2017C_17Nov2017-v1:1of4", "JetHT_Run2017C_17Nov2017-v1:2of4", "JetHT_Run2017C_17Nov2017-v1:3of4", "JetHT_Run2017C_17Nov2017-v1:4of4", 
					  "JetHT_Run2017D_17Nov2017-v1:1of2", "JetHT_Run2017D_17Nov2017-v1:2of2", 
					  "JetHT_Run2017E_17Nov2017-v1:1of3", "JetHT_Run2017E_17Nov2017-v1:2of3", "JetHT_Run2017E_17Nov2017-v1:3of3", 
					  "JetHT_Run2017F_17Nov2017-v1:1of4", "JetHT_Run2017F_17Nov2017-v1:2of4", "JetHT_Run2017F_17Nov2017-v1:3of4", "JetHT_Run2017F_17Nov2017-v1:4of4"};
int nJetTrees[nJetRunRanges]    = {242, 242, 531, 531, 531, 531, 271, 271, 435, 435, 435, 591, 591, 591, 591};

const int nMonteCarloPtHatRanges = 13;
TString MonteCarloPtHatRange[nMonteCarloPtHatRanges] = {    "Pt-15to20",      "Pt-20to30",       "Pt-30to50",       "Pt-50to80",       "Pt-80to120",   "Pt-120to170",   "Pt-170to300:1of2",   "Pt-170to300:2of2",             "Pt-300to470", "Pt-470to600",    "Pt-600to800",          "Pt-800to1000",     "Pt-1000toInf"};
double CrossSection[nMonteCarloPtHatRanges]          = {1.27319E9*0.003, 5.58528E8*0.0053, 1.39803E8*0.01182, 1.92225E7*0.02276, 2.758420E6*0.03844, 469797.*0.05362, 117989.*0.07335, 117989.*0.07335,         7820.25*0.10196,        645.528*0.12242,   187.109*0.13412,         32.3486*0.14552,       10.4305*0.15544};
double GeneratedEvents[nMonteCarloPtHatRanges]       = {5748924, 27932241, 28592087, 21074359, 23028753, 20177642, 46227426, 11174145, 1023477, 3658286, 7572383, 7404958};
int nMonteCarloTrees[nMonteCarloPtHatRanges]         = {120, 414, 488, 288, 351, 377, 797, 797, 202, 406, 69, 169, 180};


// QCD inclusive 13 TeV 
const int nMCInclusivePtHatRanges = 15;                                       
TString MCInclusivePtHatRange[nMCInclusivePtHatRanges] = {"Pt-15to30", "Pt-30to50", "Pt-50to80",    "Pt-80to120", "Pt-120to170", "Pt-170to300",  "Pt-300to470", "Pt-470to600", "Pt-600to800", "Pt-800to1000", "Pt-1000to1400", "Pt-1400to1800", "Pt-1800to2400", "Pt-2400to3200", "Pt-3200toInf"};
double CrossSectionInclusive[nMCInclusivePtHatRanges]  = {  1.83741E9,   1.40932E8,   1.92043E7,        2762530.,       471100.,       117276.,          7823.,         648.2,         186.9,         32.293,          9.4183,         0.84265,        0.114943,      0.00682981,    0.000165445};
double GeneratedEventsInclusive[nMCInclusivePtHatRanges] = {19792480, 19731387, 16260363, 28780949, 24746121, 29791318, 23834796, 26384434, 35796318, 30253528, 19465166, 0, 0, 0, 0};
int nMCInclusiveTrees[nMCInclusivePtHatRanges]           = {289, 260, 215, 420, 432, 478, 452, 575, 662, 590, 348, 0, 0, 0, 0};

// https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py
float PileUpScenario[] = { 3.39597497605e-05,
              6.63688402133e-06,
              1.39533611284e-05,
              3.64963078209e-05,
              6.00872171664e-05,
              9.33932578027e-05,
              0.000120591524486,
              0.000128694546198,
              0.000361697233219,
              0.000361796847553,
              0.000702474896113,
              0.00133766053707,
              0.00237817050805,
              0.00389825605651,
              0.00594546732588,
              0.00856825906255,
              0.0116627396044,
              0.0148793350787,
              0.0179897368379,
              0.0208723871946,
              0.0232564170641,
              0.0249826433945,
              0.0262245860346,
              0.0272704617569,
              0.0283301107549,
              0.0294006137386,
              0.0303026836965,
              0.0309692426278,
              0.0308818046328,
              0.0310566806228,
              0.0309692426278,
              0.0310566806228,
              0.0310566806228,
              0.0310566806228,
              0.0307696426944,
              0.0300103336052,
              0.0288355370103,
              0.0273233309106,
              0.0264343533951,
              0.0255453758796,
              0.0235877272306,
              0.0215627588047,
              0.0195825559393,
              0.0177296309658,
              0.0160560731931,
              0.0146022004183,
              0.0134080690078,
              0.0129586991411,
              0.0125093292745,
              0.0124360740539,
              0.0123547104433,
              0.0123953922486,
              0.0124360740539,
              0.0124360740539,
              0.0123547104433,
              0.0124360740539,
              0.0123387597772,
              0.0122414455005,
              0.011705203844,
              0.0108187105305,
              0.00963985508986,
              0.00827210065136,
              0.00683770076341,
              0.00545237697118,
              0.00420456901556,
              0.00367513566191,
              0.00314570230825,
              0.0022917978982,
              0.00163221454973,
              0.00114065309494,
              0.000784838366118,
              0.000533204105387,
              0.000358474034915,
              0.000238881117601,
              0.0001984254989,
              0.000157969880198,
              0.00010375646169,
              6.77366175538e-05,
              4.39850477645e-05,
              2.84298066026e-05,
              1.83041729561e-05,
              1.17473542058e-05,
              7.51982735129e-06,
              6.16160108867e-06,
              4.80337482605e-06,
              3.06235473369e-06,
              1.94863396999e-06,
              1.23726800704e-06,
              7.83538083774e-07,
              4.94602064224e-07,
              3.10989480331e-07,
              1.94628487765e-07,
              1.57888581037e-07,
              1.2114867431e-07,
              7.49518929908e-08,
              4.6060444984e-08,
              2.81008884326e-08,
              1.70121486128e-08,
			   1.02159894812e-08};
