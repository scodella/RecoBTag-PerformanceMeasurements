TString TreeContentFlag = "";

// ttbar MC
TString EOSPathttbar = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL18/MC/TTToHadronic_TuneCP5_13TeV-powheg-pythia8_RunIISummer19UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1/";
const int nttbarRanges = 1;
TString ttbarRange[nttbarRanges] = { "hadronic" };
int nttbarTrees[nttbarRanges] = { 2546 };
double CrossSectionttbar[nttbarRanges] = { 313.9 };
double GeneratedEventsttbar[nttbarRanges] = { 124421000 };

// QCDMu MC
TString EOSPathQCDMu = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL18/MC/QCD_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIISummer19UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1/";
const int nMonteCarloPtHatRanges = 14;
TString MonteCarloPtHatRange[nMonteCarloPtHatRanges] = { "Pt-20to30", "Pt-15to20", "Pt-170to300:1of2", "Pt-170to300:2of2", "Pt-600to800", "Pt-120to170", "Pt-800to1000", "Pt-50to80", "Pt-470to600", "Pt-30to50", "Pt-80to120", "Pt-1000toInf", "Pt-300to470:1of2", "Pt-300to470:2of2" };
int nMonteCarloTrees[nMonteCarloPtHatRanges] = { 436, 88, 670, 670, 391, 408, 421, 283, 441, 429, 568, 406, 987, 987 };
double CrossSection[nMonteCarloPtHatRanges] = { 5.58528E8*0.0053, 1.27319E9*0.003, 117989.*0.07335, 117989.*0.07335, 187.109*0.13412, 469797.*0.05362, 32.3486*0.14552, 1.92225E7*0.02276, 645.528*0.12242, 1.39803E8*0.01182, 2.758420E6*0.03844, 10.4305*0.15544, 7820.25*0.10196, 7820.25*0.10196 };
double GeneratedEvents[nMonteCarloPtHatRanges] = { 28355569, 3934872, 36762384, 36762384, 17081824, 21129958, 17187614, 20087151, 20681561, 30932276, 25473562, 14675000, 49234866, 49234866 };

// JetHT Data
TString EOSPathJetHT = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL18/Data/JetHT_Run2018-12Nov2019_UL2018/";
const int nJetRunRanges = 6;
TString JetRunRangeName[nJetRunRanges] = { "Run2018D:1of3", "Run2018D:2of3", "Run2018D:3of3", "Run2018C", "Run2018A", "Run2018B" };
int nJetTrees[nJetRunRanges] = { 1266, 1266, 1266, 262, 554, 262 };

// BTagMu Data
TString EOSPathBTagMu = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL18/Data/BTagMu_Run2018-12Nov2019_UL2018/";
const int nBTagMuRanges = 6;
TString BTagMuRangeName[nBTagMuRanges] = { "Run2018B", "Run2018C", "Run2018A", "Run2018D:1of3", "Run2018D:2of3", "Run2018D:3of3" };
int nBTagMuTrees[nBTagMuRanges] = { 262, 262, 554, 1258, 1258, 1258 };

// QCD MC:
TString EOSPathQCD = "root://eoscms.cern.ch//eos/cms/store/group/phys_btag/performance/UL18/MC/QCD_TuneCP5_13TeV_pythia8_RunIISummer19UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1/";
const int nMCInclusivePtHatRanges = 22;
TString MCInclusivePtHatRange[nMCInclusivePtHatRanges] = { "Pt-80to120", "Pt-3200toInf", "Pt-1400to1800", "Pt-600to800:1of3", "Pt-600to800:2of3", "Pt-600to800:3of3", "Pt-50to80", "Pt-170to300", "Pt-120to170:1of2", "Pt-120to170:2of2", "Pt-300to470:1of3", "Pt-300to470:2of3", "Pt-300to470:3of3", "Pt-470to600:1of2", "Pt-470to600:2of2", "Pt-15to30", "Pt-1000to1400", "Pt-2400to3200", "Pt-30to50", "Pt-800to1000:1of2", "Pt-800to1000:2of2", "Pt-1800to2400" };
int nMCInclusiveTrees[nMCInclusivePtHatRanges] = { 157, 37, 171, 1465, 1465, 1465, 455, 500, 642, 642, 1089, 1089, 1089, 719, 719, 274, 505, 72, 432, 778, 778, 94 };
double CrossSectionInclusive[nMCInclusivePtHatRanges] = { 2762530., 0.000165445, 0.84265, 186.9, 186.9, 186.9, 1.92043E7, 117276., 471100., 471100., 7823., 7823., 7823., 648.2, 648.2, 1.83741E9, 9.4183, 0.00682981, 1.40932E8, 32.293, 32.293, 0.114943 };
double GeneratedEventsInclusive[nMCInclusivePtHatRanges] = { 6937200, 722400, 5892300, 67514500, 67514500, 67514500, 19983800, 28322900, 20994700, 20994700, 57616400, 57616400, 57616400, 27343500, 27343500, 19995200, 19397100, 1992800, 19976600, 39090300, 39090300, 2990400 };

// https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/SimGeneral/MixingModule/python/mix_2018_25ns_UltraLegacy_PoissonOOTPU_cfi.py
float PileUpScenario[] = { 8.89374611122e-07, 1.1777062868e-05, 3.99725585118e-05, 0.000129888015252, 0.000265224848687,
                           0.000313088635109, 0.000353781668514, 0.000508787237162, 0.000873670065767, 0.00147166880932,
                           0.00228230649018, 0.00330375581273, 0.00466047608406, 0.00624959203029, 0.00810375867901,
                           0.010306521821, 0.0129512453978, 0.0160303925502, 0.0192913204592, 0.0223108613632,
                           0.0249798930986, 0.0273973789867, 0.0294402350483, 0.031029854302, 0.0324583524255,
                           0.0338264469857, 0.0351267479019, 0.0360320204259, 0.0367489568401, 0.0374133183052,
                           0.0380352633799, 0.0386200967002, 0.039124376968, 0.0394201612616, 0.0394673457109,
                           0.0391705388069, 0.0384758587461, 0.0372984548399, 0.0356497876549, 0.0334655175178,
                           0.030823567063, 0.0278340752408, 0.0246009685048, 0.0212676009273, 0.0180250593982,
                           0.0149129830776, 0.0120582333486, 0.00953400069415, 0.00738546929512, 0.00563442079939,
                           0.00422052915668, 0.00312446316347, 0.00228717533955, 0.00164064894334, 0.00118425084792,
                           0.000847785826565, 0.000603466454784, 0.000419347268964, 0.000291768785963, 0.000199761337863,
                           0.000136624574661, 9.46855200945e-05, 6.80243180179e-05, 4.94806013765e-05, 3.53122628249e-05,
                           2.556765786e-05, 1.75845711623e-05, 1.23828210848e-05, 9.31669724108e-06, 6.0713272037e-06,
                           3.95387384933e-06, 2.02760874107e-06, 1.22535149516e-06, 9.79612472109e-07, 7.61730246474e-07,
                           4.2748847738e-07, 2.41170461205e-07, 1.38701083552e-07, 3.37678010922e-08, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0 };

