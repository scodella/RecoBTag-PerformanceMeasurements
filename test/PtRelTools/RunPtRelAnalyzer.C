#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
R__LOAD_LIBRARY(PtRelAnalyzer_C.so) 
//R__LOAD_LIBRARY(BTagCalibrationStandalone_cpp.so)

TString TemplateVariable, PUWeighting, KinWeighting, Selection;
  PtRelAnalyzer *ptRelAna;

void ExecEventCounter(TString DataType) {

  cout << "ExecEventCounter " << DataType << endl; 

  if (DataType.Contains("QCDMu"))
    for (int ph = 0; ph<nMonteCarloPtHatRanges; ph++)
      ptRelAna->EventCounter("QCDMu", ph);

  if (DataType.Contains("QCDX"))
    for (int ph = 0; ph<nMCInclusivePtHatRanges; ph++)
      ptRelAna->EventCounter("QCD", ph);

}

void ExecComputeBTaggingWorkingPoints(TString AlgorithmName) {

  cout << "ExecComputeBTaggingWorkingPoints " << AlgorithmName << endl; 

  bool UseQCD = (AlgorithmName.Contains("UseTTBar")) ? false : true;
  AlgorithmName.ReplaceAll("UseTTBar", "");

  bool RemovePileUpJets = (AlgorithmName.Contains("RemovePU")) ? true : false;
  AlgorithmName.ReplaceAll("RemovePU", "");

  bool ApplyPileUpReweighting = (AlgorithmName.Contains("ReweightPU")) ? true : false;
  AlgorithmName.ReplaceAll("ReweightPU", "");

  ptRelAna->ComputeBTaggingWorkingPoints(AlgorithmName, UseQCD, RemovePileUpJets, ApplyPileUpReweighting);

}

void ExecComputePileUpWeights(TString Operation, int DataRange) {

  cout << "ExecComputePileUpWeights " << Operation << " " << DataRange << endl; 

  if (Operation.Contains("Loop") || Operation.Contains("All"))
    ptRelAna->ComputePileUpWeights("", "Loop" + Operation, DataRange);

  if (Operation.Contains("Merge") || Operation.Contains("All"))
    ptRelAna->ComputePileUpWeights("", "Merge" + Operation, DataRange);

  if (Operation.Contains("TestPV")) {
 
    ptRelAna->ComputePileUpWeights("TestPV", "QCDMu");

  } else if (Operation.Contains("_PV") || Operation.Contains("_PSV")) {

    TString DataType = Operation; DataType.Replace(0, DataType.First(":"), "");
    ptRelAna->ComputePileUpWeights(Operation.ReplaceAll(DataType, ""), DataType);

  } else if (Operation.Contains("PS")) {
    
    ptRelAna->ComputePileUpWeights("PS", Operation);
    ptRelAna->ComputePileUpWeights("PSm05", Operation);
    ptRelAna->ComputePileUpWeights("PSp05", Operation);

  }
      
}

void ExecProduceHistograms(TString DataType, int DataRange) {

  cout << "ExecProduceHistograms " << DataType << " " << DataRange << endl; 

  int nDataRanges;
  if (DataType=="BTagMu") nDataRanges = nBTagMuRanges;
  if (DataType=="QCDMu") nDataRanges = nMonteCarloPtHatRanges;
  
  if (TemplateVariable!="System8") {

    ptRelAna->BookHistograms();
  
    for (int dr = 0; dr<nDataRanges; dr++)
      if (DataRange==-999 || DataRange==dr) 
	ptRelAna->FillHistograms(DataType, dr);
      
    if (nDataRanges>1 && DataRange<0)
      ptRelAna->MergeHistograms(DataType);
    
  } else {

    ptRelAna->BookSystem8Histograms();
    
    for (int dr = 0; dr<nDataRanges; dr++)
      if (DataRange==-999 || DataRange==dr) 
	ptRelAna->FillSystem8Histograms(DataType, dr);

    if (nDataRanges>1 && DataRange<0)
      ptRelAna->MergeSystem8Histograms(DataType);
    
  }
  
}

void ExecProduceLightHistograms(TString DataType, int DataRange) {

  cout << "ExecProduceLightHistograms " << DataType << " " << DataRange << endl; 

  ptRelAna->BookLightHistograms();

   int nDataRanges;
  if (DataType=="JetHT") nDataRanges = nJetRunRanges;
  if (DataType=="QCD") nDataRanges = nMCInclusivePtHatRanges;

  for (int dr = 0; dr<nDataRanges; dr++) 
    if (DataRange==-999 || DataRange==dr)
      ptRelAna->FillLightHistograms(DataType, dr);

  if (nDataRanges>1 && DataRange<0)
    ptRelAna->MergeLightHistograms(DataType, 0, nDataRanges);
  
}

void ExecBuildTemplates(TString TrackTemplates) {

  cout << "ExecBuildTemplates " << TrackTemplates << endl; 
 
  if (TemplateVariable=="PtRel") {

    bool AddTrackTemplates = (TrackTemplates.Contains("All")) ? true : false;
    ptRelAna->BuildTemplates(AddTrackTemplates);

  } else if (TemplateVariable=="System8") ptRelAna->BuildSystem8Templates();
  
}

void ExecComputeKinematicWeights(TString DataType) {

  cout << "ExecComputeKinematicWeights " << DataType << endl; 
  
  for (int is = 0; is<nSystematics; is++) {

    TString ThisSystematic = SystematicName[is]; ThisSystematic.ReplaceAll("_", "");
    if (KinWeighting.Contains(ThisSystematic)) {

      if (TemplateVariable=="PtRel") {

	TString AddTrackTemplates = (DataType.Contains("All")) ? "All" : "";
	if (DataType.Contains("QCDMu")) ptRelAna->ComputeKinematicWeights(SystematicName[is], "QCDMu", AddTrackTemplates); 
	if (DataType.Contains("JetHT")) ptRelAna->ComputeKinematicWeights(SystematicName[is], "JetHT",             "All"); 
	if (DataType.Contains( "QCDX")) ptRelAna->ComputeKinematicWeights(SystematicName[is],   "QCD",             "All"); 

      } 
 
      if (TemplateVariable=="System8") ptRelAna->ComputeKinematicWeights(SystematicName[is], "QCDMu", ""); 

    }

  }

}

void ExecCompareDataToMC(TString EtaBin, int Rebinning) {

  cout << "ExecCompareDataToMC " << EtaBin << " " << Rebinning << endl; 

  for (int nb = 0; nb<nPtRelEtaBins; nb++) {
    if (EtaBin.Contains(PtRelEtaBin[nb])) {
      
      TString LightTemplates = "All";
      if (TemplateVariable=="System8") LightTemplates = "";

      ptRelAna->CompareDataToMC("jetPt",  PtRelEtaBin[nb], "All", Rebinning,   LightTemplates);
      ptRelAna->CompareDataToMC("jetEta", PtRelEtaBin[nb], "All", 2*Rebinning, LightTemplates);
      ptRelAna->CompareDataToMC("PV",     PtRelEtaBin[nb], "All", Rebinning,   LightTemplates);
      ptRelAna->CompareDataToMC("muonPt", PtRelEtaBin[nb], "All", Rebinning,   LightTemplates);
      
    }
  }

}

void ExecComputePtRelScaleFactors(TString Parameters = "CSVv2:_Central:") {

  cout << "ExecComputePtRelScaleFactors " << Parameters << endl; 
 
  TString TaggerList = Parameters;
  TaggerList.Remove(TaggerList.First(":")); 
  TString SystFit = Parameters;
  SystFit.ReplaceAll(TaggerList + ":", "");
  SystFit.Remove(SystFit.First(":")); 
  TString FitOption = Parameters;
  FitOption.ReplaceAll(TaggerList + ":" + SystFit + ":", "");

  TString PlotOption = "png";
  TString PtBinList = "All";
  TString EtaBinFlag = "anyEta";
  TString TemplateFlag = (FitOption.Contains("_LightTemplatesRatio")) ? "All" : "All";

  if (SystFit=="_All") {

    for (int is = 0; is<nFitSystematics; is++) {

      //if (FitSystematicName[is].Contains("_Central_")) continue;

      if (!FitSystematicName[is].Contains("_Central_")) PlotOption = "";
      else PlotOption = "png";

      ptRelAna->ComputePtRelScaleFactors(TaggerList, FitSystematicName[is], FitOption, PlotOption, PtBinList, EtaBinFlag, TemplateFlag);

    } 

  } else {
 
    ptRelAna->ComputePtRelScaleFactors(TaggerList, SystFit, FitOption, PlotOption, PtBinList, EtaBinFlag, TemplateFlag);

  }

}

void ExecPlotBTagPerformance(TString Dependence) {

  cout << "ExecPlotBTagPerformance " << Dependence << endl; 

  TString DrawedSystematics = "NONE";

  TString EtaList[2] = {"anyEta", "-"};

  if (Dependence=="ICHEP2016") {

    TString ConfigurationList[2] = {"_PSICHEP2016_KinPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[2] = {"_Central", "-"};

    for (int tg = 0; tg<nTaggers; tg++) 
      ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 670, DrawedSystematics);//, "ScaleFactors");

  } else if (Dependence=="FitOption") {

    TString ConfigurationList[2] = {"_KinPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[6] = {"_Central_TCF", /*"_Central_LightTemplatesRatio", "_Central_cJets", "_Central_bTempRatioCorr",*/
				 "_Central_LightTemplatesRatio_cJets_TCF" , "_Central_LightTemplatesRatio_bTempRatioCorr_TCF", 
				   "_Central_cJets_bTempRatioCorr_TCF", "_Central_LightTemplatesRatio_cJets_bTempRatioCorr_TCF", "-"};

    for (int tg = 0; tg<nTaggers; tg++) 
      if (TaggerName[tg]=="DeepCSVM")
	ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 670, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="TrgConf") {

    TString ConfigurationList[3] = {"_KinPtBinsCentral_LowPtAwayTrgConf", "_KinPtBinsCentral_LowPtAway", "-"};
    TString SystematicList[2] = {"_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "-"};

    for (int tg = 0; tg<nTaggers; tg++) 
      if (TaggerName[tg].Contains("DeepCSV"))
	ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 670, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="Run2017") {

    TString ConfigurationList[4] = {"_PSRun2017BMoriond18_KinPtBinsCentral_LowPtAwayTrgConf", "_PSRun2017CDEMoriond18_KinPtBinsCentral_LowPtAwayTrgConf", "_PSRun2017EFMoriond18_KinPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[2] = {"_Central_LightTemplatesRatio_bTempRatioCorr_cJets", "-"};

    for (int tg = 0; tg<nTaggers; tg++) 
      if (TaggerName[tg].Contains("DeepCSV"))
	ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 670, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="bJets") {

    TString ConfigurationList[2] = {"_KinPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[4] = {"_Central_LightTemplatesRatio_cJets", "_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "_Central_LightTemplatesRatio__cJets_bJets", "-"};

    for (int tg = 0; tg<nTaggers; tg++) 
      if (TaggerName[tg].Contains("DeepCSV"))
	ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 670, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="Moriond18") {

    TString ConfigurationList[2] = {"_KinPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[2] = {"_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "-"};

    for (int tg = 0; tg<nTaggers; tg++) 
      if (TaggerName[tg].Contains("DeepCSV"))
	ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 670, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="Prompt18") {

    TString ConfigurationList[2] = {"_KinPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[4] = {"_Central_LightTemplatesRatio", "_Central_LightTemplatesRatio_bTempRatioCorr", "_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "-"};

    for (int tg = 0; tg<nTaggers; tg++) 
      if (TaggerName[tg].Contains("Deep"))
	ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 670, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="KinEta") {

    TString ConfigurationList[4] = {"_KinPtBinsCentral_LowPtAwayTrgConf", "_KinPtInEtaBinsCentral_LowPtAwayTrgConf", "_KinEtaAfterPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[2] = {"_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "-"};

    for (int tg = 0; tg<nTaggers; tg++) 
      if (TaggerName[tg].Contains("Deep"))
	ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 670, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="2018Reco") {

    TString ConfigurationList[3] = {"_PSRun2018Prompt18_KinEtaAfterPtBinsCentral_LowPtAwayTrgConf", "_PSRun20182018Ultimate_KinEtaAfterPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[2] = {"_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "-"};

    for (int tg = 0; tg<nTaggers; tg++) 
      if (TaggerName[tg].Contains("Deep"))
	ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 1070, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="Nonoalgo") {

    TString ConfigurationList[3] = {"_PSRun2018ABC2018Ultimate_KinEtaAfterPtBinsCentral_LowPtAwayTrgConf", "_PSRun2018D2018Ultimate_KinEtaAfterPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[3] = {"_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "-"};

    for (int tg = 0; tg<nTaggers; tg++)
      if (TaggerName[tg].Contains("Deep"))
        ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 1070, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="UL17") {

    TString ConfigurationList[2] = {"_PSRun2017UL17_KinEtaAfterPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[2] = {"_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "-"};

    for (int tg = 0; tg<nTaggers; tg++)
      if (TaggerName[tg].Contains("Deep"))
        ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 1070, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="UL18") {
                                                                                      
    TString ConfigurationList[2] = {"_PSRun2018UL18_KinEtaAfterPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[4] = {"_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "_Central_LightTemplatesRatio_bTempRatioCorr", "_Central_LightTemplatesRatio_cJets", "-"};

    for (int tg = 0; tg<nTaggers; tg++)
      if (TaggerName[tg].Contains("Deep"))
        ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 1070, DrawedSystematics, "ScaleFactors");

  } else if (Dependence=="Efficiency2018") {

    TString ConfigurationList[3] = {"_PSRun2018Prompt18_KinEtaAfterPtBinsCentral_LowPtAwayTrgConf", "_PSRun2018UL18_KinEtaAfterPtBinsCentral_LowPtAwayTrgConf", "-"};
    TString SystematicList[2] = {"_Central_LightTemplatesRatio_cJets_bTempRatioCorr", "-"};

    for (int tg = 0; tg<nTaggers; tg++)
      if (TaggerName[tg].Contains("Deep"))
        ptRelAna->PlotBTagPerformance(Dependence, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 1070, DrawedSystematics, "Efficiency");

  }
  
}

void ExecAnalyzeSystematics(TString FitFlag = "") {

  cout << "ExecAnalyzeSystematics " << FitFlag << endl; 

  TString EtaList[2] = {"anyEta", "-"};
  TString ConfigurationList[2] = {KinWeighting + "_LowPtAwayTrgConf", "-"};
  
  for (int sfs = 1; sfs<nScaleFactorSystematics; sfs++) {
    
    TString SystematicList[4] = {"_Central" + FitFlag, FitSystematicName[2*sfs-1] + FitFlag, FitSystematicName[2*sfs] + FitFlag, "-"};
    
    TString DepOn = ScaleFactorSystematicName[sfs]; DepOn.ReplaceAll("_", "");
    
    for (int tg = 1; tg<nTaggers; tg++) 
      ptRelAna->PlotBTagPerformance(DepOn, TaggerName[tg], EtaList, ConfigurationList, SystematicList, 1070, "NONE");
  
  }

  TString SystematicList[2] = {"_Central" + FitFlag, "-"};
  
  for (int tg = 1; tg<nTaggers; tg++)
    ptRelAna->PlotBTagPerformance("Final", TaggerName[tg], EtaList, ConfigurationList, SystematicList, 1070, "_Central" + FitFlag, "ScaleFactors");

}

void ExecStoreScaleFactors(TString FitFlag = "") {

  cout << "ExecStoreScaleFactors " << FitFlag << endl; 

  //string BTagger = "DeepCSV"; TString BTaggerName = BTagger;
  string BTagger = "DeepJet"; TString BTaggerName = BTagger;

  TString OP = "All";

  BTagCalibration calib(BTagger);

  const int nOperatingPoints = 3;
  string OperatingPoint[nOperatingPoints] = {"L", "M", "T"};
  
  string MeasName = "ptrel";

  float MinDisc = 0., MaxDisc = 1.;
  if (BTagger=="JP") MaxDisc = 5.; 

  for (int op = 0; op<nOperatingPoints; op++) 
    if (OP=="All" || OP.Contains(OperatingPoint[op])) {
		  
      BTagEntry::OperatingPoint wp;
      if (OperatingPoint[op]=="L")  wp = BTagEntry::OP_LOOSE;
      if (OperatingPoint[op]=="M") wp = BTagEntry::OP_MEDIUM;
      if (OperatingPoint[op]=="T")  wp = BTagEntry::OP_TIGHT;
      
      ptRelAna->ComputeScaleFactorSystematics(BTagger + OperatingPoint[op], "anyEta", "_Central" + FitFlag);
    
      for (int ThisFlavour = 5; ThisFlavour<=5; ThisFlavour++) {
	
	BTagEntry::JetFlavor jfl;
	if (ThisFlavour==4) jfl = BTagEntry::FLAV_C;
	if (ThisFlavour==5) jfl = BTagEntry::FLAV_B;
	  
	int SysFactor = 6 - ThisFlavour;

	for (int fpt = 0; fpt<nFitPtBins; fpt++) {
	      
	  float MinPt = FitPtEdge[fpt];
	  float MaxPt = FitPtEdge[fpt+1];
					   
	  for (int uu = 0; uu<=nScaleFactorSystematics+1; uu++) 
	    for (int vs = 0; vs<2; vs++) {
	      
	      string SysFlag = "central";
	      
	      float ScaledSF = ScaleFactorValue[fpt][0];
	     	      
	      if (uu==0 && vs==1) continue;

	      if (uu>0) {
	     
		float ScalingSF = 0.;
		
		if (vs==0) SysFlag = "up";
		else if (vs==1) SysFlag = "down";
		
		if (uu==1) ScalingSF = SysFactor*TotalScaleFactorSystematic[fpt][0];
		else {
		  
		  SysFlag += "_"; SysFlag += ScaleFactorSystematicStoredName[uu-2];
		  ScalingSF = SysFactor*ScaleFactorSystematic[fpt][0][uu-2];
		
		}

		if (vs==0) ScaledSF += ScalingSF;
		else if (vs==1) ScaledSF -= ScalingSF;

	      }
	      
	      TString ThisFunction; ThisFunction += ScaledSF;
		  
	      BTagEntry::Parameters params(wp, MeasName, SysFlag, jfl, -2.4, 2.4, MinPt, MaxPt, MinDisc, MaxDisc);
	     
	      const TF1 SFFun("SFFun", ThisFunction, MinPt, MaxPt);
	      BTagEntry e1(&SFFun, params);  
	      
	      calib.addEntry(e1);
	      
	    }

	}

      }

    }
  
  ofstream outFile; outFile.open(BTaggerName + ".csv");
  calib.makeCSV(outFile);
  outFile.close();

}

void RunPtRelAnalyzer() {

  cout << "RunPtRelAnalyzer" << endl; 

  TemplateVariable = gSystem->Getenv("TEMPLATEVARIABLE");
  PUWeighting = gSystem->Getenv("PUWEIGHTING");
  KinWeighting = gSystem->Getenv("KINWEIGHTING");
  Selection = gSystem->Getenv("SELECTION");

  ptRelAna = new PtRelAnalyzer(TemplateVariable, PUWeighting, KinWeighting, Selection); 

  TString Macro = gSystem->Getenv("MACRONAME");
  TString DatasetIndex = gSystem->Getenv("DATARANGEINDEX");
  int DataRange = DatasetIndex.Atoi() - 1; 

  if (Macro.Contains("EventCounter"))
    ExecEventCounter(Macro.ReplaceAll("EventCounter", "")); 

  if (Macro.Contains("ComputeBTaggingWorkingPoints"))
    ExecComputeBTaggingWorkingPoints(Macro.ReplaceAll("ComputeBTaggingWorkingPoints", "")); 

  if (Macro.Contains("ComputePileUpWeights"))
    ExecComputePileUpWeights(Macro.ReplaceAll("ComputePileUpWeights", ""), DataRange);

  if (Macro.Contains("ProduceHistograms"))
    ExecProduceHistograms(Macro.ReplaceAll("ProduceHistograms", ""), DataRange);
  if (Macro.Contains("ProduceLightHistograms"))
    ExecProduceLightHistograms(Macro.ReplaceAll("ProduceLightHistograms", ""), DataRange);
  if (Macro.Contains("BuildTemplates")) ExecBuildTemplates(Macro);

  if (Macro.Contains("ComputeKinematicWeights")) ExecComputeKinematicWeights(Macro);
  if (Macro.Contains("CompareDataToMC")) 
    ExecCompareDataToMC(Macro, DataRange);

  if (Macro.Contains("ComputePtRelScaleFactors"))
    ExecComputePtRelScaleFactors(Macro.ReplaceAll("ComputePtRelScaleFactors", ""));
  if (Macro.Contains("PlotBTagPerformance"))
    ExecPlotBTagPerformance(Macro.ReplaceAll("PlotBTagPerformance", ""));
  if (Macro.Contains("AnalyzeSystematics"))
    ExecAnalyzeSystematics(Macro.ReplaceAll("AnalyzeSystematics", ""));

  if (Macro.Contains("StoreScaleFactors"))
    ExecStoreScaleFactors(Macro.ReplaceAll("StoreScaleFactors", ""));
    
}
 
