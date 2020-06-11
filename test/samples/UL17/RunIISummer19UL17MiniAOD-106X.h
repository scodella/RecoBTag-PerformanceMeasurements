// BTagMu Data
TString EOSPathBTagMu = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL17/Data/BTagMu_Run2017-09Aug2019_UL2017/";
const int nBTagMuRanges = 5;
TString BTagMuRangeName[nBTagMuRanges] = { "Run2017F", "Run2017E", "Run2017D", "Run2017C", "Run2017B" };
int nBTagMuTrees[nBTagMuRanges] = { 585, 429, 271, 516, 243 };

// QCDMu MC
TString EOSPathQCDMu = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL17/MC/QCD_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6/";
const int nMonteCarloPtHatRanges = 13;
TString MonteCarloPtHatRangeName[nMonteCarloPtHatRanges] = { "Pt-20to30", "Pt-50to80", "Pt-600to800", "Pt-120to170", "Pt-800to1000", "Pt-170to300:1of2", "Pt-170to300:2of2", "Pt-470to600", "Pt-30to50:1of2", "Pt-30to50:2of2", "Pt-80to120", "Pt-1000toInf", "Pt-300to470" };
int nMonteCarloTrees[nMonteCarloPtHatRanges] = { 439, 359, 0, 24, 552, 973, 973, 27, 609, 609, 22, 444, 29 };
double CrossSection[nMonteCarloPtHatRanges] = { 5.58528E8*0.0053, 1.92225E7*0.02276, 187.109*0.13412, 469797.*0.05362, 32.3486*0.14552, 117989.*0.07335, 117989.*0.07335, 645.528*0.12242, 1.39803E8*0.01182, 1.39803E8*0.01182, 2.758420E6*0.03844, 10.4305*0.15544, 7820.25*0.10196 };
double GeneratedEvents[nMonteCarloPtHatRanges] = { 28365947, 20937742, 0, 648944, 16995944, 36918785, 36918785, 517382, 30992450, 30992450, 613257, 14719636, 494796 };

// JetHT Data
TString EOSPathJetHT = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL17/Data/JetHT_Run2017-09Aug2019_UL2017/";
const int nJetRunRanges = 5;
TString JetRunRangeName[nJetRunRanges] = { "Run2017E", "Run2017D", "Run2017F", "Run2017C", "Run2017B" };
int nJetTrees[nJetRunRanges] = { 0, 0, 0, 0, 0 };

// QCD MC
TString EOSPathQCD = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL17/MC/QCD_TuneCP5_13TeV_pythia8_RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6/";
const int nMCInclusivePtHatRanges = 21;
TString MCInclusivePtHatRangeName[nMCInclusivePtHatRanges] = { "Pt-80to120", "Pt-3200toInf", "Pt-1400to1800", "Pt-600to800:1of4", "Pt-600to800:2of4", "Pt-600to800:3of4", "Pt-600to800:4of4", "Pt-50to80", "Pt-170to300", "Pt-120to170", "Pt-300to470:1of3", "Pt-300to470:2of3", "Pt-300to470:3of3", "Pt-470to600", "Pt-15to30", "Pt-1000to1400", "Pt-2400to3200", "Pt-30to50", "Pt-800to1000:1of2", "Pt-800to1000:2of2", "Pt-1800to2400" };
int nMCInclusiveTrees[nMCInclusivePtHatRanges] = { 531, 47, 173, 1575, 1575, 1575, 1575, 291, 553, 482, 1152, 1152, 1152, 595, 268, 439, 78, 310, 861, 861, 109 };
double CrossSectionInclusive[nMCInclusivePtHatRanges] = { 2762530., 0.000165445, 0.84265, 186.9, 186.9, 186.9, 186.9, 1.92043E7, 117276., 471100., 7823., 7823., 7823., 648.2, 1.83741E9, 9.4183, 0.00682981, 1.40932E8, 32.293, 32.293, 0.114943 };
double GeneratedEventsInclusive[nMCInclusivePtHatRanges] = { 29640900, 800000, 5434800, 65128700, 65128700, 65128700, 65128700, 19989398, 29522100, 29951400, 57365500, 57365500, 57365500, 27559600, 19997400, 19967700, 1919400, 19082198, 39261600, 39261600, 2999700 };

float PileUpScenario[] = { }
