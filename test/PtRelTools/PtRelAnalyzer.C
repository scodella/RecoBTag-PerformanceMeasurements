#include "PtRelAnalyzer.h"
#include <TROOT.h>
#include <TStyle.h>

#include "BTagAnalyzerUtilities.C"

PtRelAnalyzer::PtRelAnalyzer(TString TemplateFlag, TString PileUpReweighting, TString KinematicWeighting, TString SelectionFlag) {

  PUWeighting = PileUpReweighting;

  KinWeighting = KinematicWeighting;

  Selection = SelectionFlag;
  if (Selection=="") Selection = "_";
  Selection.ReplaceAll("_", "_" + TriggerTreatment);

  
  TemplateVariable = TemplateFlag;

  if (TemplateVariable=="PtRel") {

    nBinsForTemp = 100;
    LowerEdgeForTemp = 0.;
    UpperEdgeForTemp = 5.;

  } else if (TemplateVariable=="IP3D") {
    
    nBinsForTemp = 100.;
    LowerEdgeForTemp = -10.;
    UpperEdgeForTemp = 0.;

  } else if (TemplateVariable=="IP3DSig") {

    nBinsForTemp = 120.;
    LowerEdgeForTemp = -6.;
    UpperEdgeForTemp = 6.;

  } else if (TemplateVariable=="SoftMu") {
    
    nBinsForTemp = 140.;
    LowerEdgeForTemp = -0.2;
    UpperEdgeForTemp = 1.2;

  } else if (TemplateVariable=="System8") {

    nBinsForTemp = 70;
    LowerEdgeForTemp = 0.;
    UpperEdgeForTemp = 7.;
    
  }

  //GetBTemplateCorrections();

  //GetDeltaRCorrections();

  //GetMuPtWeight();

  //ReadJetEnergyUncertainty("AK5PFchs");

  //ReadReshapeFunction();
  
}

PtRelAnalyzer::~PtRelAnalyzer() {

}

TString PtRelAnalyzer::HistogramName(TString VariableName, TString PtBin, int EtaBin, int TriggerIdx, int SystematicIdx, int TaggerIdx, int FlavourIdx) {

  TString ThisHistogramName = VariableName + "_" + PtBin + "_" + PtRelEtaBin[EtaBin];

  if (TriggerIdx>=0) ThisHistogramName += TriggerName[TriggerIdx];
  if (SystematicIdx>=0) ThisHistogramName += SystematicName[SystematicIdx];
  if (TaggerIdx>=0) { ThisHistogramName += "_"; ThisHistogramName += TaggerName[TaggerIdx]; }

  if (FlavourIdx>=0) {

    if (FlavourIdx<10) ThisHistogramName += "_Tag";
    else ThisHistogramName += "_Untag";
    
    if (FlavourIdx==0 || FlavourIdx==10)      ThisHistogramName += "";
    else if (FlavourIdx==1 || FlavourIdx==11) ThisHistogramName += "_b";
    else if (FlavourIdx==2 || FlavourIdx==12) ThisHistogramName += "_c";
    else if (FlavourIdx==3 || FlavourIdx==13) ThisHistogramName += "_lg";
    else if (FlavourIdx==4 || FlavourIdx==14) ThisHistogramName += "_trk";
    else cout << "PtRelAnalyzer::HistogramName: bad choice for paratemeter FlavourIdx " << FlavourIdx << endl;

    if (FlavourIdx==4 || FlavourIdx==14) {
      
      ThisHistogramName.ReplaceAll("_Data", "_Jet");
      ThisHistogramName.ReplaceAll( "_QCD", "_Inc");
      
    }
    
  }

  if ((VariableName.Contains("_Jet") || VariableName.Contains("_Inc")) && SystematicName[SystematicIdx]=="_TCHPM" && TriggerIdx>=0)
    ThisHistogramName.ReplaceAll("_TCHPM", "_JBPM4");

  return ThisHistogramName;
  
}

TString PtRelAnalyzer::HistogramName(TString VariableName, int PtBin, int EtaBin, int TriggerIdx, int SystematicIdx, int TaggerIdx, int FlavourIdx) {

  return HistogramName(VariableName, PtRelPtBin[PtBin], EtaBin, TriggerIdx, SystematicIdx, TaggerIdx, FlavourIdx);

}

void PtRelAnalyzer::BookHistograms() {

  cout << "Booking histograms" << endl;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++) {
     
      TString ThisHistoName;

      for (int is = 0; is<nSystematics; is++) {
	
	ThisHistoName = HistogramName("Observed_DataType", ptb, nb, -1, is, -1, -1);  
	Observed[ptb][nb][is] = new TH2D(ThisHistoName, ThisHistoName, nTriggers, 0., nTriggers, nTriggers, 0., nTriggers);

      }

      for (int tr = 0; tr<nTriggers; tr++)  {
	
	for (int is = 0; is<nSystematics; is++) {

	  ThisHistoName = HistogramName("muonPt_DataType", ptb, nb, tr, is, -1, -1);	  
	  MuPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, 0., 60.); 
	  
	  ThisHistoName = HistogramName("muonDR_DataType", ptb, nb, tr, is, -1, -1);
	  MuDRForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 50, 0., 0.5); 
	  
	  ThisHistoName = HistogramName("jetPt_DataType", ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 600, 0., 600.); 
	  
	  ThisHistoName = HistogramName("jetEta_DataType", ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.);
	  
	  ThisHistoName = HistogramName("nGoodPV_DataType", ptb, nb, tr, is, -1, -1);
	  PVMultiplicity[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 35, 0., 35.); 
	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	    	   
	    for (int fl = 0; fl<2; fl++) {
 
	      ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, fl);
	      PtRelTagForWeighting[fl*nTaggers+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 

	      ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, 10 + fl);
	      PtRelUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
	      
	      ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg,  2 + fl);
	      PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
	    
	      ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, 12 + fl);
	      PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
	    	    
	    }
	  
	  }
	  
	}
	
      }

    }
      
  bool DoErrors = false;
  
  if (DoErrors)
    for (int is = 0; is<nSystematics; is++) 
      for (int tg = 0; tg<nTaggers; tg++) 
	for (int ptb = 0; ptb<nPtRelPtBins; ptb++) 
	  for (int nb = 0; nb<nPtRelEtaBins; nb++)	    
	    for (int tr = 0; tr<nTriggers; tr++) 
	      for (int dt = 0; dt<2; dt++) {
	      
		PtRelTagForWeighting[dt*nTaggers+tg][tr][is][ptb][nb]->Sumw2();
		PtRelUntagForWeighting[dt*nTaggers+tg][tr][is][ptb][nb]->Sumw2();
		PtRelLightTagForWeighting[dt*nTaggers+tg][tr][is][ptb][nb]->Sumw2();
		PtRelLightUntagForWeighting[dt*nTaggers+tg][tr][is][ptb][nb]->Sumw2();
	      
	      } 

}

void PtRelAnalyzer::ResetHistograms(TString DataType) {

  cout << "Resetting histograms" << endl;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int etab = 0; etab<nPtRelEtaBins; etab++) 
      for (int is = 0; is<nSystematics; is++) {

	TString ThisHistoName = HistogramName("Observed_" + DataType, ptb, etab, -1, is, -1, -1);
	Observed[ptb][etab][is]->SetNameTitle(ThisHistoName, ThisHistoName);
	Observed[ptb][etab][is]->Reset();
	
	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  ThisHistoName = HistogramName("muonPt_" + DataType, ptb, etab, tr, is, -1, -1);
	  MuPtForWeighting[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  MuPtForWeighting[tr][is][ptb][etab]->Reset();
	  
	  ThisHistoName = HistogramName("muonDR_" + DataType, ptb, etab, tr, is, -1, -1);
	  MuDRForWeighting[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  MuDRForWeighting[tr][is][ptb][etab]->Reset();
	  
	  ThisHistoName = HistogramName("jetPt_" + DataType, ptb, etab, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetPtForWeighting[tr][is][ptb][etab]->Reset();
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, ptb, etab, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetEtaForWeighting[tr][is][ptb][etab]->Reset();
	  
	  ThisHistoName = HistogramName("nGoodPV_" + DataType, ptb, etab, tr, is, -1, -1);
	  PVMultiplicity[tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	  PVMultiplicity[tr][is][ptb][etab]->Reset();
	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	  	    
	    for (int fl = 0; fl<2; fl++) {
         
	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, fl);
	      PtRelTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	      PtRelTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Reset();
         
	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 10 + fl);
	      PtRelUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	      PtRelUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Reset();
         
	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg,  2 + fl);
	      PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	      PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Reset(); 
         
	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 12 + fl);
	      PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->SetNameTitle(ThisHistoName, ThisHistoName);
	      PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Reset();
	      
	    }
	    
	  }
	  
	}
	
      }
  
}

void PtRelAnalyzer::FillHistograms(TString DataType, int DataRange) {
  
  ResetHistograms(DataType);
  
  cout << "Filling histograms" << endl;

  GetKinematicWeights(DataType);
  
  GetPileUpWeights(DataType);
  
  GetBTemplateCorrections();
  
  const int nAwayTaggers = 16;
  TString AwayTaggerName[nAwayTaggers] = {"TCHPM", "TCHPL", "TCHPT", "JBPT", "JBPT1", "JBPM", "JBPM1", "JBPM2", "JBPM3", "JBPM4", "JBPM5", "JBPM6", "JBPM7", "JBPM8", "JBPM9", "JBPL"};
  
  int ph = DataRange;
  
  double PtHatWeight = 1.;
  if (DataType=="QCD") PtHatWeight = CrossSection[ph]/GeneratedEvents[ph];
  
  int nEventsFromTrees = 0;
  
  TString DataRangeName;
  if (DataType=="Data") DataRangeName = BTagMuRangeName[ph];
  else if (DataType=="QCD") DataRangeName = MonteCarloPtHatRange[ph];
  
  int FirstTree = 1;
  int nTrees = nBTagMuTrees[ph]; 
  if (DataType=="QCD") nTrees = nMonteCarloTrees[ph];
  
  for (int tf = FirstTree; tf<=nTrees; tf++) {
    
    TString FileNumber = "_"; FileNumber += tf; FileNumber += ".root";
    
    TString FileName;	  
    if (DataType=="Data") FileName = EOSPathData + DataRangeName + "/JetTree_data" + TreeContentFlag + FileNumber;
    //if (DataType=="QCD") FileName = EOSPathQCDMu + DataRangeName + "/JetTree_mc" + TreeContentFlag + FileNumber;

    if (DataType=="QCD") FileName = EOSPathQCDMu + "/JetTree_mc" + TreeContentFlag + FileNumber;
    if (DataType=="QCD") FileName.ReplaceAll("Pt-", DataRangeName);
    if (FileName.Contains("15to20")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150822_131817/0000/JetTree"); 
    if (FileName.Contains("20to30")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150822_132009/0000/JetTree"); 
    if (FileName.Contains("30to50")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150821_132612/0000/JetTree"); 
    if (FileName.Contains("50to80")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150821_143819/0000/JetTree"); 
    if (FileName.Contains("80to120")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150822_132157/0000/JetTree"); 
    if (FileName.Contains("120to170")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150822_132251/0000/JetTree"); 
    if (FileName.Contains("170to300")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150822_132349/0000/JetTree"); 
    if (FileName.Contains("300to470")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150822_132445/0000/JetTree"); 
    if (FileName.Contains("470to600")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150822_132622/0000/JetTree"); 
    if (FileName.Contains("600to800")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150822_132821/0000/JetTree"); 
    if (FileName.Contains("800to1000")) FileName.ReplaceAll("/JetTree", "/QCD_" + DataRangeName + "_50ns_jec/150822_132940/0000/JetTree"); 
    if (DataType=="QCD") FileName.ReplaceAll("8/QCD_Pt-", "8/QCD_Pt");
    if (DataType=="QCD") FileName.ReplaceAll("8//QCD_Pt-", "8/QCD_Pt");
    if (DataType=="Data") FileName.ReplaceAll("/JetTree", "/150904_181757/0000/JetTree"); 
    if (DataType=="Data" && tf==93) continue;

    //if (FileName.Contains("15to20") || FileName.Contains("50to80") || FileName.Contains("80to120") || FileName.Contains("120to170") || FileName.Contains("600to800")) FileName.ReplaceAll("/JetTree", "/0000/JetTree"); 
    
    TFile *ThisTree = TFile::Open(FileName);
    
    TTree *tchain = GetChain(ThisTree, false);
    
    if (DataType=="Data") cout << "  Getting " << DataType << "   sample: BTagMu    " << DataRangeName << " File " << tf << "/" << nTrees << endl;
    if (DataType=="QCD") cout << "  Getting " << DataType << "    sample: MuPt5 " << DataRangeName << " File " << tf << "/" << nTrees << endl;
    
    Int_t nentries = (Int_t)tchain->GetEntries();
	    
    nEventsFromTrees += nentries;
    
    for (Int_t i = 0; i<nentries; i++) {
      
      tchain->GetEntry(i);
      
      int iMu; 
      int jMu = GetPFMuonJet(&iMu);
      
      if (jMu>=0 && nPV>0 && Jet_pt[jMu]<PtRelPtEdge[nPtRelPtBins] && fabs(Jet_eta[jMu])<PtRelEtaEdge[nPtRelEtaBins-1]) {
	
	int ptBin = -1;
	for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
	  if (Jet_pt[jMu]>PtRelPtEdge[ptb]) ptBin = ptb;
	
	int tr = -1;
	
	bool PassTrigger[nTriggers];
	
	for (int trg = 0; trg<nTriggers; trg++) {
	  
	  PassTrigger[trg] = false;
	  
	  if (Jet_pt[jMu]>MinPtJetTrigger[trg] && Jet_pt[jMu]<MaxPtJetTrigger[trg]) {
	    
	    bool FiredTrigger = false;
	    if (DataType=="Data" || !Selection.Contains("TrgEmul")) {
	      if (PassTriggerBit(TriggerName[trg])) FiredTrigger = true;
	    } else {
	      if (PassTriggerEmulation(trg, jMu)) FiredTrigger = true;
	    }
	    
	    if (FiredTrigger) {
	      
	      PassTrigger[trg] = true;
	      tr = trg;
	      
	    }
	    
	  }
	  
	}
	
	if (DataType=="QCD") 
	  if (!PassPtHat("MuEnrichedPt5", jMu)) tr = -1;
				     
	if (tr>=0) {

	  bool JetAwayPass[nSystematics], EventJetAwayPass = false;
	  for (int is = 0; is<nSystematics; is++) {
	    
	    TString ThisAwayTaggerName = "JBPM4";//"TCHPM";
	    //if (Jet_pt[jMu]<50.) ThisAwayTaggerName = "JBPM7"; 
	    if (Jet_pt[jMu]>=140.) ThisAwayTaggerName = "JBPM5"; 
	    for (int at = 0; at<nAwayTaggers; at++) 
	      if (SystematicName[is].Contains(AwayTaggerName[at]))
		ThisAwayTaggerName = AwayTaggerName[at];
	    
	    int aJet = GetAwayJet(ThisAwayTaggerName, jMu, 1.5, true);
	    
	    JetAwayPass[is] = false;
	    
	    if (aJet>=0) {

	      float PtAwayJetCut = PtAwayTrigger[tr];

	      if (TriggerName[tr]=="_DiJet20" && Jet_pt[jMu]>=50.) PtAwayJetCut = 30.;

	      if (Jet_pt[aJet]>PtAwayJetCut) {

		JetAwayPass[is] = true;
		EventJetAwayPass = true;
		
	      }
	    
	    }

	  }
	  
	  if (EventJetAwayPass) {
	    
	    int SavePtBin = ptBin;
		    
	    for (int is = 0; is<nSystematics; is++) {
	      
	      ptBin = SavePtBin;
	      
	      float TrackPtCut =  MuonPtCut[ptBin];
	      if (SystematicName[is]=="_MuPt6") TrackPtCut = 6.;
	      if (SystematicName[is]=="_MuPt8") TrackPtCut = 8.;
	      
	      float TrackDRCut = 999., TrackMinDRCut = 0.;
	      GetDRCuts(SystematicName[is], Jet_pt[jMu], &TrackDRCut, &TrackMinDRCut);
	      
	      double JetMuonDR = DeltaR(Jet_eta[jMu], Jet_phi[jMu], PFMuon_eta[iMu], PFMuon_phi[iMu]);
	      
	      if (PFMuon_pt[iMu]>TrackPtCut && JetMuonDR<TrackDRCut  && JetMuonDR>=TrackMinDRCut && JetAwayPass[is]) {
		
		double ThisJetWeight = 1.;

		if (DataType=="QCD") {

		  int ipt = Jet_pt[jMu]; 
		  if (ipt>=PtRelPtEdge[nPtRelPtBins]) ipt = PtRelPtEdge[nPtRelPtBins] - 1;
		  int imu = PFMuon_pt[iMu]; if (imu>=60) imu = 59;
		  int ieta = 0; if (fabs(Jet_eta[jMu])>1.2) ieta = 1;

		  int EventPileUp = nPUtrue;
		  if (PUWeighting.Contains("PV")) EventPileUp = nPV;
		  if (EventPileUp>=nMaxPU) EventPileUp = nMaxPU - 1;

		  ThisJetWeight = PtHatWeight*KinematicWeight[tr][is][ipt][imu][ieta]*PileUpWeight[tr][EventPileUp][is];

		}

		if (DataType=="QCD" && SystematicName[is].Contains("_GluonSplitting")) {
		  if (fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_GluonSplittingB")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 5)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==4 && SystematicName[is].Contains("_GluonSplittingC")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 4)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==21 && SystematicName[is].Contains("_GluonSplittingG")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 21)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } 
		}
		
		if (DataType=="QCD" && fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_BTemplates")) {
		  int iB = GetBHadron(jMu);
		  if (iB>=0) {
		    float EnergyFraction = BHadron_pT[iB]/Jet_genpt[jMu]; int efbin = EnergyFraction/0.02;
		    if (efbin>=0 && efbin<100) {
		      if (SystematicName[is].Contains("Minus")) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][0];
		      if (SystematicName[is].Contains("Plus") ) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][1];
		    }
		  }
		}
		/*
		if (DataType=="QCD" && SystematicName[is].Contains("_MuPtWeight")) {
		  
		  int ThisMuonPtBin = PFMuon_pt[iMu]; if (ThisMuonPtBin>=60) ThisMuonPtBin = 59;
		  ThisJetWeight *= MuPtWeight[ptBin][ThisMuonPtBin];
		  
		  }*/
			
		/*if (SystematicName[is].Contains("_JEU") && dt==1) {
		  
		  float ThisJEU = GetJEU(1, Jet_eta[jMu], Jet_pt[jMu]);
		  if (SystematicName[is].Contains("Twice")) ThisJEU *= 2.;
		  if (SystematicName[is].Contains("JEU2")) ThisJEU = sqrt(ThisJEU*ThisJEU + 0.02*0.02);
		  
		  int JEUSign = 1.; if (SystematicName[is].Contains("Down")) JEUSign = -1.; 
		  float ScaledJetPt = Jet_pt[jMu]*(1. + JEUSign*ThisJEU);
		  
		  ptBin = -1;
		  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
		    if (ScaledJetPt>PtRelPtEdge[ptb]) ptBin = ptb;
		  
		}
		
		if (ptBin==-1) return;*/
			
		bool JetEtaBin[nPtRelEtaBins]; JetEtaBin[0] = true;
		for (int jeta = 1; jeta<nPtRelEtaBins; jeta++) {

		  JetEtaBin[jeta] = false;
		  if (fabs(Jet_eta[jMu])>PtRelEtaEdge[jeta-1] && fabs(Jet_eta[jMu])<PtRelEtaEdge[jeta]) JetEtaBin[jeta] = true;
		  
		}
		
		for (int nb = 0; nb<nPtRelEtaBins; nb++) 
		  if (JetEtaBin[nb]) {
		    
		    double Discriminator = PFMuon_ptrel[iMu];
		    		    
		    if (TemplateVariable=="IP3D")
		      Discriminator = log(fabs(PFMuon_IP[iMu]));
		    
		    if (TemplateVariable=="LinearIP3D") 
		      Discriminator = fabs(PFMuon_IP[iMu]);
		    
		    if (TemplateVariable=="IP3DSig")
		      Discriminator = log(fabs(PFMuon_IP[iMu]));
		    
		    if (TemplateVariable=="JetProb")	      
		      Discriminator = Jet_Proba[jMu];		    
		    
		    if (TemplateVariable=="SoftMu")
		      Discriminator = Jet_SoftMu[jMu];	
		    
		    if (Discriminator==0.) continue;

		    if (DataType=="Data")
		      for (int trg1 = 0; trg1<nTriggers; trg1++) 
			for (int trg2 = 0; trg2<nTriggers; trg2++) 
			  if (PassTrigger[trg1] || PassTrigger[trg2]) 
			    Observed[ptBin][nb][is]->Fill(trg1, trg2); 	     
		    
		    PVMultiplicity[tr][is][ptBin][nb]->Fill(nPV, ThisJetWeight);      
		    MuDRForWeighting[tr][is][ptBin][nb]->Fill(JetMuonDR, ThisJetWeight);    
		    MuPtForWeighting[tr][is][ptBin][nb]->Fill(PFMuon_pt[iMu], ThisJetWeight);
		    JetPtForWeighting[tr][is][ptBin][nb]->Fill(Jet_pt[jMu], ThisJetWeight); 
		    JetEtaForWeighting[tr][is][ptBin][nb]->Fill(Jet_eta[jMu], ThisJetWeight);
		    
		    for (int tg = 0; tg<nTaggers; tg++) 
		      if (IsTaggedJet(jMu, TaggerName[tg])) {
			
			if (DataType=="Data") 
			  PtRelTagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else if (fabs(Jet_flavour[jMu])==5)	  
			  PtRelTagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
		        else if (fabs(Jet_flavour[jMu])==4) 
			  PtRelLightTagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else 	  
			  PtRelLightTagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
				
		      } else {
						
			if (DataType=="Data")
			  PtRelUntagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else if (fabs(Jet_flavour[jMu])==5)
			  PtRelUntagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			
			else if (fabs(Jet_flavour[jMu])==4)
			  PtRelLightUntagForWeighting[tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			else 	  
			  PtRelLightUntagForWeighting[nTaggers+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			
		      }
			    
		  } // Loop on eta bins
		
	      } // Muon pt and DR cuts
	      
	    } // Loop on systematics
	    
	  } // Has jet away
	  
	} // Pass trigger
	
      } // Muon-jet selection
      
    } // Loop on events
    
    ThisTree->Close();
	    
  } // Loop on files
  
  cout << "  nEventsFromTrees " << nEventsFromTrees << endl;
  
  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataRangeName + PUWeighting + KinWeighting + Selection + ".root";

  if (DataType=="Data") OutputFileName.ReplaceAll(PUWeighting + KinWeighting, "");
  
  SaveHistograms(OutputFileName, DataType);

}

void PtRelAnalyzer::SaveHistograms(TString OutputFileName, TString DataType) {

  cout << "Saving histograms" << endl;

  TFile *OutFile = new TFile(OutputFileName, "recreate");

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int etab = 0; etab<nPtRelEtaBins; etab++) 
      for (int is = 0; is<nSystematics; is++) {

	if (DataType=="Data") 
	  Observed[ptb][etab][is]->Write();

	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  MuDRForWeighting[tr][is][ptb][etab]->Write();
	  MuPtForWeighting[tr][is][ptb][etab]->Write();
	  JetPtForWeighting[tr][is][ptb][etab]->Write();
	  JetEtaForWeighting[tr][is][ptb][etab]->Write();
	  PVMultiplicity[tr][is][ptb][etab]->Write();
	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	      
	    if (DataType=="Data") {

	      PtRelTagForWeighting[tg][tr][is][ptb][etab]->Write();
	      PtRelUntagForWeighting[tg][tr][is][ptb][etab]->Write();
	    
	    } else if (DataType=="QCD") {
		
	      PtRelTagForWeighting[nTaggers+tg][tr][is][ptb][etab]->Write();
	      PtRelUntagForWeighting[nTaggers+tg][tr][is][ptb][etab]->Write();

	      for (int fl = 0; fl<2; fl++) {
		  
		PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Write(); 
		PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Write();
		  
	      }
	      
	    }
	    
	  }
	  
	}

      }
  cout << OutputFileName << endl;
  OutFile->Close();

  cout << "Exiting" << endl;
  
}
 
void PtRelAnalyzer::MergeHistograms(TString DataType, int FirstBin, int LastBin) {

  ResetHistograms(DataType);

  cout << "Merging histograms: " << DataType << endl;

  if (LastBin==100) {

    if (DataType=="Data") LastBin = nBTagMuRanges - 1;
    else if (DataType=="QCD") LastBin = nMonteCarloPtHatRanges - 1;

  }
  
  for (int ph = FirstBin; ph<=LastBin; ph++) {
    
    TString DataRangeName;
    if (DataType=="Data") DataRangeName = BTagMuRangeName[ph];
    else if (DataType=="QCD") DataRangeName = MonteCarloPtHatRange[ph];

    cout << "  Merging " << DataRangeName << endl;
  
    TString InputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataRangeName + PUWeighting + KinWeighting + Selection + ".root";
    if (DataType=="Data") InputFileName.ReplaceAll(PUWeighting + KinWeighting, "");
    TFile *HistogramFile = TFile::Open(InputFileName); 

    for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {

      cout << "    Merging " << PtRelPtBin[ptb] << endl;

      for (int etab = 0; etab<nPtRelEtaBins; etab++) 
	for (int is = 0; is<nSystematics; is++) {
	  
	  TString ThisHistoName;

	  if (DataType=="Data") {
	    
	    ThisHistoName = HistogramName("Observed_" + DataType, ptb, etab, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) HistogramFile->Get(ThisHistoName);
	    Observed[ptb][etab][is]->Add(ThisObserved);
	  
	  }

	  for (int tr = 0; tr<nTriggers; tr++) {
	    
	    ThisHistoName = HistogramName("muonPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisMuPt = (TH1D*) HistogramFile->Get(ThisHistoName);
	    MuPtForWeighting[tr][is][ptb][etab]->Add(ThisMuPt);
	    
	    ThisHistoName = HistogramName("muonDR_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisMuDR = (TH1D*) HistogramFile->Get(ThisHistoName);
	    MuDRForWeighting[tr][is][ptb][etab]->Add(ThisMuDR);
	    
	    ThisHistoName = HistogramName("jetPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetPt = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetPtForWeighting[tr][is][ptb][etab]->Add(ThisJetPt);
	    
	    ThisHistoName = HistogramName("jetEta_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetEta = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetEtaForWeighting[tr][is][ptb][etab]->Add(ThisJetEta);
	    
	    ThisHistoName = HistogramName("nGoodPV_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisPVMult = (TH1D*) HistogramFile->Get(ThisHistoName);
	    PVMultiplicity[tr][is][ptb][etab]->Add(ThisPVMult);
	    
	    for (int tg = 0; tg<nTaggers; tg++) {
	      
	      if (DataType=="Data") {
		
		ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg,  0);
		TH1D *ThisTag = (TH1D*) HistogramFile->Get(ThisHistoName);
		PtRelTagForWeighting[tg][tr][is][ptb][etab]->Add(ThisTag);
		
		ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 10);
		TH1D *ThisUntag = (TH1D*) HistogramFile->Get(ThisHistoName);
		PtRelUntagForWeighting[tg][tr][is][ptb][etab]->Add(ThisUntag);
		
	      } else if (DataType=="QCD") {
		
		
		ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg,  1);
		TH1D *ThisTagB = (TH1D*) HistogramFile->Get(ThisHistoName);
		PtRelTagForWeighting[nTaggers+tg][tr][is][ptb][etab]->Add(ThisTagB);
		
		ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 11);
		TH1D *ThisUntagB = (TH1D*) HistogramFile->Get(ThisHistoName);
		PtRelUntagForWeighting[nTaggers+tg][tr][is][ptb][etab]->Add(ThisUntagB);
		
		for (int fl = 0; fl<2; fl++) {
		  
		
		  ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg,  2 + fl);
		  TH1D *ThisTagL = (TH1D*) HistogramFile->Get(ThisHistoName);
		  PtRelLightTagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Add(ThisTagL);
		
		  ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 12 + fl);
		  TH1D *ThisUntagL = (TH1D*) HistogramFile->Get(ThisHistoName); 
		  PtRelLightUntagForWeighting[fl*nTaggers+tg][tr][is][ptb][etab]->Add(ThisUntagL);
		  
		}
		
	      }
	      
	    }
	    
	  }
	  
	}

    }

  }

  TString LastBinFlag = "";
  if (DataType=="QCD") {

    if (MonteCarloPtHatRange[LastBin]=="Pt-300470") LastBinFlag = "PtHat470";

  }
  
  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataType + PUWeighting + KinWeighting + Selection + LastBinFlag + ".root";

  if (DataType=="Data") OutputFileName.ReplaceAll(PUWeighting + KinWeighting, "");
  
  SaveHistograms(OutputFileName, DataType);
  
}

void PtRelAnalyzer::BookLightHistograms() {
  
  cout << "Booking light histograms" << endl;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      for (int tr = 0; tr<nTriggers; tr++) 
	for (int is = 0; is<nSystematics; is++)	{
	  
	  TString ThisHistoName = HistogramName("jetPt_DataType", ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 600, 0., 600.); 
	  
	  ThisHistoName = HistogramName("jetEta_DataType", ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.);
  
	  for (int tg = 0; tg<nTaggers; tg++) {
	      
	    ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, 4);
	    PtRelLightTagForSystematic[tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
	    
	    ThisHistoName = HistogramName(TemplateVariable + "_DataType", ptb, nb, tr, is, tg, 14);
	    PtRelLightUntagForSystematic[tg][tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
 
	    bool DoErrors = false;
  
	    if (DoErrors) {
	      
	      PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->Sumw2();
	      PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->Sumw2();
	      
	    }
	    
	  } 

        } 
  
}

void PtRelAnalyzer::ResetLightHistograms(TString DataType) {

  cout << "Resetting light histograms" << endl; 

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++) 
      for (int is = 0; is<nSystematics; is++)	
	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  TString ThisHistoName = HistogramName("jetPt_" + DataType, ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetPtForWeighting[tr][is][ptb][nb]->Reset();
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetEtaForWeighting[tr][is][ptb][nb]->Reset();
	  
	  for (int tg = 0; tg<nTaggers; tg++) {
	      
	    ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, nb, tr, is, tg, 4);
	    PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);;
	    PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->Reset();
	    
	    ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, nb, tr, is, tg, 14);
	    PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);;
	    PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->Reset();
	    
	  }

        } 
  
}

void PtRelAnalyzer::FillLightHistograms(TString DataType, int DataRange) {

  ResetLightHistograms(DataType);

  cout << "Filling light histograms" << endl;
  
  GetKinematicWeights(DataType);
  GetPileUpWeights(DataType);
  
  double PtHatWeight = 1.;
  if (DataType=="Inc") PtHatWeight = CrossSectionInclusive[DataRange]/GeneratedEventsInclusive[DataRange];
             
  int nEventsFromTrees = 0;
      
  TString DataRangeName;
  if (DataType=="Jet") DataRangeName = JetRunRangeName[DataRange];
  else if (DataType=="Inc") DataRangeName = MCInclusivePtHatRange[DataRange];

  int nTrees = nJetTrees[DataRange]; int FirstTree = 1;
  if (DataType=="Inc") nTrees = nMCInclusiveTrees[DataRange];
 
  for (int tf = FirstTree; tf<nTrees+FirstTree; tf++) {
	
    TString FileNumber = "_"; FileNumber += tf; FileNumber += ".root";
    
    TString FileName;
    if (DataType=="Jet") FileName = EOSPathJet + JetRunRangeName [DataRange] + "/JetTree_data" + TreeContentFlag + FileNumber;
    if (DataType=="Inc") FileName = EOSPathQCD + MCInclusivePtHatRange[DataRange] + "/JetTree_mc" + TreeContentFlag + FileNumber;
	
    if (DataType=="Jet")
      if (tf==38 || tf==67 || tf==93) continue;
	
    if (DataType=="Inc") 
      if (DataRangeName.Contains("300to470") && tf==346) continue;
 
    TFile *ThisTree = TFile::Open(FileName);
	
    TTree *tchain = GetChain(ThisTree, true);
    
    if (DataType=="Jet") cout << "  Getting data   sample: Jet       " << JetRunRangeName[DataRange] << " " << tf << "/" << nTrees << endl; 
    if (DataType=="Inc") cout << "  Getting QCD    sample: Inclusive " << MCInclusivePtHatRange[DataRange] << " " << tf << "/" << nTrees << endl;
	
    Int_t nentries = (Int_t)tchain->GetEntries();

    nEventsFromTrees += nentries;

    for (Int_t i = 0; i<nentries; i++) {

      tchain->GetEntry(i);

      for (int ijet = 0; ijet<nJet; ijet++) {
	    
	if (Jet_pt[ijet]>PtRelPtEdge[0] && !HasPFMuon(ijet, false) && !HasTaggedJet(ijet) && 
	    nPV>0 && Jet_pt[ijet]<PtRelPtEdge[nPtRelPtBins] && fabs(Jet_eta[ijet])<PtRelEtaEdge[nPtRelEtaBins-1]) {
	      
	  int ptBin = -1;
	  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
            if (Jet_pt[ijet]>PtRelPtEdge[ptb]) ptBin = ptb;
	   
	  int IdxAwayJet = GetAwayJet("NONE", ijet, 1.5, false);
	      
          if (IdxAwayJet<0) continue;
	      
	  int tr = -1;

	  for (int trg = 0; trg<nTriggers; trg++)
            if (Jet_pt[IdxAwayJet]>PtAwayTrigger[trg] && Jet_pt[ijet]>MinPtJetTrigger[trg] && Jet_pt[ijet]<MaxPtJetTrigger[trg]) {
		  
       	      bool PassTrigger = PassTriggerBit(JetTriggerName[trg]);
		  
	      if (trg==1 && Jet_pt[ijet]>50. && Jet_pt[ijet]< 80. && i%2!=1) PassTrigger = false;
              if (trg==1 && Jet_pt[ijet]>80. && Jet_pt[ijet]<140. && i%3!=1) PassTrigger = false;
	      if (trg==2 && Jet_pt[ijet]>80. && Jet_pt[ijet]<140. && i%3!=2) PassTrigger = false;
		
	      if (PassTrigger) tr = trg;
		  
	   }
		  
	   if (DataType=="Inc") 
             if (!PassPtHat("Inclusive", ijet)) tr = -1;
	      
	   if (tr>=0) {
	     
	     int SavePtBin = ptBin;
		
	     for (int is = 0; is<nSystematics; is++) {
		  
	       ptBin = SavePtBin;
		  
               double ThisJetWeight = 1.;
	       int ipt = Jet_pt[ijet]; if (ipt>=600) ipt = 599;
	       int ieta = 0; if (fabs(Jet_eta[ijet])>1.2) ieta = 1;
	       int EventPileUp = nPUtrue;
               if (PUWeighting.Contains("PV") || DataType=="Jet") EventPileUp = nPV;
	       if (EventPileUp>=nMaxPU) EventPileUp = nMaxPU - 1;
	       ThisJetWeight = PtHatWeight*KinematicWeight[tr][is][ipt][0][ieta]*PileUpWeight[tr][EventPileUp][is];
		  
	       float TrackPtCut = MuonPtCut[ptBin];
	       if (SystematicName[is]=="_MuPt6") TrackPtCut = 6.;
	       if (SystematicName[is]=="_MuPt8") TrackPtCut = 8.;
		  
	       /*if (SystematicName[is].Contains("_JEU") && dt==1) {
		    
		    float ThisJEU = GetJEU(1, Jet_eta[ijet], Jet_pt[ijet]);
		    if (SystematicName[is].Contains("Twice")) ThisJEU *= 2.;
		    if (SystematicName[is].Contains("JEU2")) ThisJEU = sqrt(ThisJEU*ThisJEU + 0.02*0.02);
		    
		    int JEUSign = 1.; if (SystematicName[is].Contains("Down")) JEUSign = -1.; 
		    float ScaledJetPt = Jet_pt[ijet]*(1. + JEUSign*ThisJEU);
		    
		    ptBin = -1;
		    for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
		      if (ScaledJetPt>PtRelPtEdge[ptb]) ptBin = ptb;
		    
		  }
		  
		  if (ptBin==-1) return;
		  */
	       float TrackDRCut = 999., TrackMinDRCut = 0.;
	       GetDRCuts(SystematicName[is], Jet_pt[ijet], &TrackDRCut, &TrackMinDRCut);
		  
	       bool JetEtaBin[nPtRelEtaBins], FilledJet[nPtRelEtaBins]; 
	       JetEtaBin[0] = true; FilledJet[0] = false; 
	       for (int jeta = 1; jeta<nPtRelEtaBins; jeta++) {
		 JetEtaBin[jeta] = false; FilledJet[jeta] = false; 
                 if (fabs(Jet_eta[ijet])>PtRelEtaEdge[jeta-1] && fabs(Jet_eta[ijet])<PtRelEtaEdge[jeta]) JetEtaBin[jeta] = true;
               }

	       int nLightTracks = 0;//Jet_nLastTrkInc[ijet] - Jet_nFirstTrkInc[ijet];
		  
	       for (int tk = Jet_nFirstTrkInc[ijet]; tk<Jet_nLastTrkInc[ijet]; tk++) {
		    
		 double JetTrkDR = DeltaR(Jet_eta[ijet], Jet_phi[ijet], TrkInc_eta[tk], TrkInc_phi[tk]);
		    
		 if (TrkInc_pt[tk]>TrackPtCut && JetTrkDR<TrackDRCut && JetTrkDR>=TrackMinDRCut) nLightTracks++;
		 //if (fabs(Jet_eta[ijet]-TrkInc_eta[tk])<TrackDRCut && DeltaPhi(Jet_phi[ijet], TrkInc_phi[tk])<TrackDRCut) nLightTracks++;
		    
	       }
	     
	       for (int tk = Jet_nFirstTrkInc[ijet]; tk<Jet_nLastTrkInc[ijet]; tk++) {
		    
		 double JetTrkDR = DeltaR(Jet_eta[ijet], Jet_phi[ijet], TrkInc_eta[tk], TrkInc_phi[tk]);
		    
		 if (TrkInc_pt[tk]>TrackPtCut && JetTrkDR<TrackDRCut && JetTrkDR>=TrackMinDRCut) {
		 //if (TrkInc_pt[tk]>TrackPtCut && fabs(Jet_eta[ijet]-TrkInc_eta[tk])<TrackDRCut && DeltaPhi(Jet_phi[ijet], TrkInc_phi[tk])<TrackDRCut) {
		      
		 double ThisTrackWeight = ThisJetWeight;

		 /*if (DataType=="Inc" && ApplyDeltaRCorrections) {
		 int idrpt = 0, idreta = 0, idrdr = 0;
			float drjetpt = Jet_pt[ijet], drjeteta = fabs(Jet_eta[ijet]);
			if (drjetpt<30) idrpt = 0;
			else if (drjetpt<50) idrpt = 1;
			else if (drjetpt<80) idrpt = 2;
			else if (drjetpt<120) idrpt = 3;
			else if (drjetpt<160) idrpt = 4;
			else if (drjetpt<320) idrpt = 5;
			else if (drjetpt<500) idrpt = 6;
			else if (drjetpt<800) idrpt = 7;
			if (drjeteta<0.4) idreta = 0;
			else if (drjeteta<0.8) idreta = 1;
			else if (drjeteta<1.2) idreta = 2;
			else if (drjeteta<1.8) idreta = 3;
			else if (drjeteta<12.4) idreta = 4;
			idrdr = JetTrkDR/0.02;
			int iDEta = fabs(Jet_eta[ijet]-TrkInc_eta[tk])/0.01;
			//if (iDEta>=0 && iDEta<50)
			//ThisTrackWeight *= DeltaRCorrections[iDEta][0][idrpt][idreta];
			int ThisDPhi = fabs(Jet_phi[ijet]-TrkInc_phi[tk]);
			if (ThisDPhi>=3.14159) ThisDPhi = 2*3.14159 - ThisDPhi;
			int iDPhi = ThisDPhi/0.01;
			//if (iDPhi>=0 && iDPhi<50)
			//ThisTrackWeight *= DeltaRCorrections[iDPhi][0][idrpt][idreta];
			if (idrdr<50)
			  ThisTrackWeight *= DeltaRCorrections[idrdr][0][idrpt][idreta];
		      }*/
		      
		 double Discriminator = TrkInc_ptrel[tk];
		      
		 if (TemplateVariable=="IP2D")
		   Discriminator = log(fabs(TrkInc_IP[tk]));
		      
		 if (TemplateVariable=="IP2DSig")
		   Discriminator = log(fabs(TrkInc_IPsig[tk]));
		      
		 if (TemplateVariable=="IP3D")
	           Discriminator = log(fabs(TrkInc_IP[tk]));
		      
		 if (TemplateVariable=="LinearIP3D")
	           Discriminator = TrkInc_IP[tk];
		      
		 if (TemplateVariable=="IP3DSig")
	           Discriminator = log(fabs(TrkInc_IPsig[tk]));
		      		      
		 for (int nb = 0; nb<nPtRelEtaBins; nb++) 
	           if (JetEtaBin[nb]) {

                     if (!FilledJet[nb]) {			  
		       JetPtForWeighting[tr][is][ptBin][nb]->Fill(Jet_pt[ijet], ThisJetWeight);
                       JetEtaForWeighting[tr][is][ptBin][nb]->Fill(Jet_eta[ijet], ThisJetWeight);
                       FilledJet[nb] = true;
		
                     }

	             PtRelLightUntagForSystematic[nTaggers-1][tr][is][ptBin][nb]->Fill(Discriminator, ThisTrackWeight/nLightTracks);
			  
	           } // Loop on eta bins
		      
		 } // Muon pt and DR cuts
			
	       } // Loop on tracks
		      
	     } // Loop on systematics
		
	   } // Pass trigger
	    
	 } // Light jet selection
	    
       } // Loop on jets
	  
     } // Loop on events
	
     ThisTree->Close();
	
   } // Loop on files
      
   cout << "  nEventsFromTrees " << nEventsFromTrees << endl;
    
   for (int tr = 0; tr<nTriggers; tr++)
     for (int ptb = 0; ptb<nPtRelPtBins; ptb++) 
       for (int nb = 0; nb<nPtRelEtaBins; nb++)
	 for (int is = 0; is<nSystematics; is++) 
	   for (int tg = 0; tg<nTaggers; tg++) {

             PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->Add(PtRelLightUntagForSystematic[nTaggers-1][tr][is][ptb][nb]);

             if (tg<nTaggers-1)
	       PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->Add(PtRelLightUntagForSystematic[nTaggers-1][tr][is][ptb][nb]);

	   }

  cout << "  nEventsFromTrees " << nEventsFromTrees << endl;
  
  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_" + DataRangeName + PUWeighting + KinWeighting + Selection + ".root";
  
  SaveLightHistograms(OutputFileName);
  
}

void PtRelAnalyzer::SaveLightHistograms(TString TemplateFileName) {

  cout << "Saving light histograms" << endl;

  TFile *OutFile = new TFile(TemplateFileName, "recreate");

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      for (int tr = 0; tr<nTriggers; tr++) 
	for (int is = 0; is<nSystematics; is++)	{

	  JetPtForWeighting[tr][is][ptb][nb]->Write();
	  JetEtaForWeighting[tr][is][ptb][nb]->Write();

	  for (int tg = 0; tg<nTaggers; tg++) {
	      
	    PtRelLightTagForSystematic[tg][tr][is][ptb][nb]->Write();
	    PtRelLightUntagForSystematic[tg][tr][is][ptb][nb]->Write();
	    
	  }

        }
  
  OutFile->Close();

  cout << "Exiting" << endl;

}

void PtRelAnalyzer::MergeLightHistograms(TString DataType, int FirstBin, int LastBin) {

  ResetLightHistograms(DataType);

  cout << "Merging light histograms: " << DataType << endl;

  if (LastBin==100) {

    if (DataType=="Jet") LastBin = nJetRunRanges - 1;
    else if (DataType=="Inc") LastBin = nMCInclusivePtHatRanges - 1;

  }
  
  for (int ph = FirstBin; ph<=LastBin; ph++) {
    
    TString DataRangeName;
    if (DataType=="Jet") DataRangeName = JetRunRangeName[ph];
    else if (DataType=="Inc") DataRangeName = MCInclusivePtHatRange[ph];

    cout << "  Merging " << DataRangeName << endl;
  
    TString InputFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_" + DataRangeName + PUWeighting + KinWeighting + Selection + ".root";
    TFile *HistogramFile = TFile::Open(InputFileName); 

    for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {

      cout << "    Merging " << PtRelPtBin[ptb] << endl;

      for (int etab = 0; etab<nPtRelEtaBins; etab++) 
	for (int is = 0; is<nSystematics; is++) {
	  
	  TString ThisHistoName;

	  for (int tr = 0; tr<nTriggers; tr++) {
	    
	    ThisHistoName = HistogramName("jetPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetPt = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetPtForWeighting[tr][is][ptb][etab]->Add(ThisJetPt);
	    
	    ThisHistoName = HistogramName("jetEta_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetEta = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetEtaForWeighting[tr][is][ptb][etab]->Add(ThisJetEta);
	    
	    for (int tg = 0; tg<nTaggers; tg++) {

	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 4);
              TH1D *ThisTag = (TH1D*) HistogramFile->Get(ThisHistoName);
	      PtRelLightTagForSystematic[tg][tr][is][ptb][etab]->Add(ThisTag);

	      ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, ptb, etab, tr, is, tg, 14);
              TH1D *ThisUntag = (TH1D*) HistogramFile->Get(ThisHistoName);
	      PtRelLightUntagForSystematic[tg][tr][is][ptb][etab]->Add(ThisUntag);
	      
	    }
	    
	  }
	  
	}

    }

  }

  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_" + DataType + PUWeighting + KinWeighting + Selection + ".root";
  
  SaveLightHistograms(OutputFileName);

}

void PtRelAnalyzer::BookSystem8Histograms() {
  
  cout << "Booking System8 histograms" << endl;

  TString ThisHistoName;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++) {
      
      for (int is = 0; is<nSystematics; is++) {
	
	ThisHistoName = HistogramName("Observed_DataType", ptb, nb, -1, is, -1, -1);  
	Observed[ptb][nb][is] = new TH2D(ThisHistoName, ThisHistoName, nTriggers, 0., nTriggers, nTriggers, 0., nTriggers);
	
      }
      
      for (int tr = 0; tr<nTriggers; tr++) {
	
	for (int is = 0; is<nSystematics; is++)	{
	  
	  ThisHistoName = HistogramName("muonPt_DataType", ptb, nb, tr, is, -1, -1);	  
	  MuPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, 0., 60.); 
	  
	  ThisHistoName = HistogramName("jetPt_DataType", ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 600, 0., 600.); 
	  
	  ThisHistoName = HistogramName("jetEta_DataType", ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.);
	  
	  ThisHistoName = HistogramName("nGoodPV_DataType", ptb, nb, tr, is, -1, -1);
	  PVMultiplicity[tr][is][ptb][nb] = new TH1D(ThisHistoName, ThisHistoName, 35, 0., 35.); 
		  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {

	      if (tg==0)
		ThisHistoName = HistogramName("ntag_pT", ptb, nb, tr, is, -1, -1);
	      else ThisHistoName = HistogramName("ntag_pT", ptb, nb, tr, is, tg-1, -1);
	      if (st==1) ThisHistoName.ReplaceAll("ntag", "ptag");
	      if (tg==0) ThisHistoName.ReplaceAll("tag", "");
	      
	      System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName,         ThisHistoName,         nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
	      System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName + "_b",  ThisHistoName + "_b",  nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
	      System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb] = new TH1D(ThisHistoName + "_cl", ThisHistoName + "_cl", nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);

	    } 
	  
        }

      }
      
    }

}

void PtRelAnalyzer::ResetSystem8Histograms(TString DataType) {

  cout << "Resetting System8 histograms" << endl; 

  TString ThisHistoName;

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      for (int is = 0; is<nSystematics; is++)	{

	ThisHistoName = HistogramName("Observed_" + DataType, ptb, nb, -1, is, -1, -1);
	Observed[ptb][nb][is]->SetNameTitle(ThisHistoName, ThisHistoName);
	Observed[ptb][nb][is]->Reset();
	
	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  ThisHistoName = HistogramName("muonPt_" + DataType, ptb, nb, tr, is, -1, -1);	
	  MuPtForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  MuPtForWeighting[tr][is][ptb][nb]->Reset();
	  	  
	  ThisHistoName = HistogramName("jetPt_" + DataType, ptb, nb, tr, is, -1, -1);
	  JetPtForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetPtForWeighting[tr][is][ptb][nb]->Reset();
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, ptb, nb, tr, is, -1, -1);
	  JetEtaForWeighting[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  JetEtaForWeighting[tr][is][ptb][nb]->Reset();
	  
	  ThisHistoName = HistogramName("nGoodPV_" + DataType, ptb, nb, tr, is, -1, -1);
	  PVMultiplicity[tr][is][ptb][nb]->SetNameTitle(ThisHistoName, ThisHistoName);
	  PVMultiplicity[tr][is][ptb][nb]->Reset();
		  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {

	      System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Reset();
	      System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Reset();
	      System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Reset();

	    } 
	  
        } 

      }
  
}
 
void PtRelAnalyzer::FillSystem8Histograms(TString DataType, int DataRange) {

  ResetSystem8Histograms(DataType);

  cout << "Filling System8 histograms" << endl;

  GetKinematicWeights(DataType);
  
  GetPileUpWeights(DataType);
  
  GetBTemplateCorrections();
  
  const int nAwayTaggers = 16;
  TString AwayTaggerName[nAwayTaggers] = {"TCHPM", "TCHPL", "TCHPT", "JBPT", "JBPT1", "JBPM", "JBPM1", "JBPM2", "JBPM3", "JBPM4", "JBPM5", "JBPM6", "JBPM7", "JBPM8", "JBPM9", "JBPL"};
  
  int ph = DataRange;
  
  double PtHatWeight = 1.;
  if (DataType=="QCD") PtHatWeight = CrossSection[ph]/GeneratedEvents[ph];
  
  int nEventsFromTrees = 0;
  
  TString DataRangeName;
  if (DataType=="Data") DataRangeName = BTagMuRangeName[ph];
  else if (DataType=="QCD") DataRangeName = MonteCarloPtHatRange[ph];
  
  int FirstTree = 1;
  int nTrees = nBTagMuTrees[ph]; 
  if (DataType=="QCD") nTrees = nMonteCarloTrees[ph];
  
  for (int tf = FirstTree; tf<=nTrees; tf++) {
    
    TString FileNumber = "_"; FileNumber += tf; FileNumber += ".root";
    
    TString FileName;	  
    if (DataType=="Data") FileName = EOSPathData + DataRangeName + "/JetTree_data" + TreeContentFlag + FileNumber;
    if (DataType=="QCD") FileName = EOSPathQCDMu + DataRangeName + "/JetTree_mc" + TreeContentFlag + FileNumber;
        
    if (FileName.Contains("15to20") || FileName.Contains("50to80") || FileName.Contains("80to120") || FileName.Contains("120to170") || FileName.Contains("600to800")) FileName.ReplaceAll("/JetTree", "/0000/JetTree"); 
    
    TFile *ThisTree = TFile::Open(FileName);
    
    TTree *tchain = GetChain(ThisTree, false);
    
    if (DataType=="Data") cout << "  Getting " << DataType << "   sample: BTagMu    " << DataRangeName << " File " << tf << "/" << nTrees << endl;
    if (DataType=="QCD") cout << "  Getting " << DataType << "    sample: MuPt5 " << DataRangeName << " File " << tf << "/" << nTrees << endl;
    
    Int_t nentries = (Int_t)tchain->GetEntries();
	    
    nEventsFromTrees += nentries;
    
    for (Int_t i = 0; i<nentries; i++) {
      
      tchain->GetEntry(i);
      
      int iMu; 
      int jMu = GetPFMuonJet(&iMu);
      
      if (jMu>=0 && nPV>0 && Jet_pt[jMu]<PtRelPtEdge[nPtRelPtBins] && fabs(Jet_eta[jMu])<PtRelEtaEdge[nPtRelEtaBins-1]) {
	
	int ptBin = -1;
	for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
	  if (Jet_pt[jMu]>PtRelPtEdge[ptb]) ptBin = ptb;
	
	int tr = -1;
	
	bool PassTrigger[nTriggers];
	
	for (int trg = 0; trg<nTriggers; trg++) {
	  
	  PassTrigger[trg] = false;
	  
	  if (Jet_pt[jMu]>MinPtJetTrigger[trg] && Jet_pt[jMu]<MaxPtJetTrigger[trg]) {
	    
	    bool FiredTrigger = false;
	    if (DataType=="Data" || !Selection.Contains("TrgEmul")) {
	      if (PassTriggerBit(TriggerName[trg])) FiredTrigger = true;
	    } else {
	      if (PassTriggerEmulation(trg, jMu)) FiredTrigger = true;
	    }
	    
	    if (FiredTrigger) {
	      
	      PassTrigger[trg] = true;
	      tr = trg;
	      
	    }
	    
	  }
	  
	}
	
	if (DataType=="QCD") 
	  if (!PassPtHat("MuEnrichedPt5", jMu)) tr = -1;
	
	if (tr>=0) {

	  bool EventJetAwayPass = false;

	  int aJet = GetAwayJet("NONE", jMu, 0., false);
	  
	  if (aJet>=0) {
	    
	    float PtAwayJetCut = PtAwayTrigger[tr];
	    if (TriggerName[tr]=="_DiJet20" && Jet_pt[jMu]>=50.) PtAwayJetCut = 30.;

	    if (Jet_pt[aJet]>PtAwayJetCut)
	      EventJetAwayPass = true;
	    
	  }
	  
	  if (EventJetAwayPass) {
	    
	    int SavePtBin = ptBin;
	    
	    for (int is = 0; is<nSystematics; is++) {
	      
	      ptBin = SavePtBin;
	      
	      float TrackPtCut =  MuonPtCut[ptBin];
	      if (SystematicName[is]=="_MuPt6") TrackPtCut = 6.;
	      if (SystematicName[is]=="_MuPt8") TrackPtCut = 8.;
	      
	      float TrackDRCut = 999., TrackMinDRCut = 0.;
	      GetDRCuts(SystematicName[is], Jet_pt[jMu], &TrackDRCut, &TrackMinDRCut);
	      
	      double JetMuonDR = DeltaR(Jet_eta[jMu], Jet_phi[jMu], PFMuon_eta[iMu], PFMuon_phi[iMu]);
	      
	      if (PFMuon_pt[iMu]>TrackPtCut && JetMuonDR<TrackDRCut  && JetMuonDR>=TrackMinDRCut) {
		
		double ThisJetWeight = 1.;

		if (DataType=="QCD") {
		  
		  int ipt = Jet_pt[jMu]; 
		  if (ipt>=PtRelPtEdge[nPtRelPtBins]) ipt = PtRelPtEdge[nPtRelPtBins] - 1;
		  int imu = PFMuon_pt[iMu]; if (imu>=60) imu = 59;
		  int ieta = 0; if (fabs(Jet_eta[jMu])>1.2) ieta = 1;
		  
		  int EventPileUp = nPUtrue;
		  if (PUWeighting.Contains("PV")) EventPileUp = nPV;
		  if (EventPileUp>=nMaxPU) EventPileUp = nMaxPU - 1;
		  
		  ThisJetWeight = PtHatWeight*KinematicWeight[tr][is][ipt][imu][ieta]*PileUpWeight[tr][EventPileUp][is];
		  
		}
		
		if (DataType=="QCD" && SystematicName[is].Contains("_GluonSplitting")) {
		  if (fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_GluonSplittingB")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 5)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==4 && SystematicName[is].Contains("_GluonSplittingC")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 4)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } else if (fabs(Jet_flavour[jMu])==21 && SystematicName[is].Contains("_GluonSplittingG")) {
		    if (IsFromGluonSplittingFromHadron(jMu, 21)) {
		      if (SystematicName[is].Contains("Down")) ThisJetWeight *= 0.5;
		      if (SystematicName[is].Contains("Up")) ThisJetWeight *= 1.5;
		    }
		  } 
		}
		
		if (DataType=="QCD" && fabs(Jet_flavour[jMu])==5 && SystematicName[is].Contains("_BTemplates")) {
		  int iB = GetBHadron(jMu);
		  if (iB>=0) {
		    float EnergyFraction = BHadron_pT[iB]/Jet_genpt[jMu]; int efbin = EnergyFraction/0.02;
		    if (efbin>=0 && efbin<100) {
		      if (SystematicName[is].Contains("Minus")) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][0];
		      if (SystematicName[is].Contains("Plus") ) ThisJetWeight *= BTemplateCorrections[efbin][ptBin][1];
		    }
		  }
		}
		/*
		  if (DataType=="QCD" && SystematicName[is].Contains("_MuPtWeight")) {
		  
		  int ThisMuonPtBin = PFMuon_pt[iMu]; if (ThisMuonPtBin>=60) ThisMuonPtBin = 59;
		  ThisJetWeight *= MuPtWeight[ptBin][ThisMuonPtBin];
		  
		  }*/
			
		/*if (SystematicName[is].Contains("_JEU") && dt==1) {
		  
		  float ThisJEU = GetJEU(1, Jet_eta[jMu], Jet_pt[jMu]);
		  if (SystematicName[is].Contains("Twice")) ThisJEU *= 2.;
		  if (SystematicName[is].Contains("JEU2")) ThisJEU = sqrt(ThisJEU*ThisJEU + 0.02*0.02);
		  
		  int JEUSign = 1.; if (SystematicName[is].Contains("Down")) JEUSign = -1.; 
		  float ScaledJetPt = Jet_pt[jMu]*(1. + JEUSign*ThisJEU);
		  
		  ptBin = -1;
		  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)  
		    if (ScaledJetPt>PtRelPtEdge[ptb]) ptBin = ptb;
		  
		}
		
		if (ptBin==-1) return;*/
		
		double Discriminator = PFMuon_ptrel[iMu];
		
		if (Discriminator==0.) continue;
	       
		TString ThisAwayTaggerName = "JBPM4";//"TCHPM";
		//if (Jet_pt[jMu]<50.) ThisAwayTaggerName = "JBPM7"; 
		if (Jet_pt[jMu]>=140.) ThisAwayTaggerName = "JBPM5"; 
		for (int at = 0; at<nAwayTaggers; at++) 
		  if (SystematicName[is].Contains(AwayTaggerName[at]))
		    ThisAwayTaggerName = AwayTaggerName[at];
		
		bool IsPSet = IsTaggedJet(aJet, ThisAwayTaggerName);

		bool JetEtaBin[nPtRelEtaBins]; JetEtaBin[0] = true;
		for (int jeta = 1; jeta<nPtRelEtaBins; jeta++) {

		  JetEtaBin[jeta] = false;
		  if (fabs(Jet_eta[jMu])>PtRelEtaEdge[jeta-1] && fabs(Jet_eta[jMu])<PtRelEtaEdge[jeta]) JetEtaBin[jeta] = true;
		  
		}
		
		for (int nb = 0; nb<nPtRelEtaBins; nb++) 
		  if (JetEtaBin[nb]) {
		    
		    if (DataType=="Data")
		      for (int trg1 = 0; trg1<nTriggers; trg1++) 
			for (int trg2 = 0; trg2<nTriggers; trg2++) 
			  if (PassTrigger[trg1] || PassTrigger[trg2]) 
			    Observed[ptBin][nb][is]->Fill(trg1, trg2); 	     
		    
		    PVMultiplicity[tr][is][ptBin][nb]->Fill(nPV, ThisJetWeight);      
		    MuPtForWeighting[tr][is][ptBin][nb]->Fill(PFMuon_pt[iMu], ThisJetWeight);
		    JetPtForWeighting[tr][is][ptBin][nb]->Fill(Jet_pt[jMu], ThisJetWeight); 
		    JetEtaForWeighting[tr][is][ptBin][nb]->Fill(Jet_eta[jMu], ThisJetWeight);
		    
		    for (int st = 0; st<2; st++)
		      if (st==0 || IsPSet) 
			for (int tg = 0; tg<nTaggers+1; tg++) 
			  if (tg==0 || IsTaggedJet(jMu, TaggerName[tg-1])) {
			    
			    if (DataType=="Data")
			      System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			    else {
			      
			      if (fabs(Jet_flavour[jMu])==5) 
				System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			      
			      else if (fabs(Jet_flavour[jMu])!=5) 
				System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptBin][nb]->Fill(Discriminator, ThisJetWeight);
			      
			    }
			    
			  }
		    
		  } // Loop on eta bins
		
	      } // Muon pt and DR cuts
	      
	    } // Loop on systematics
	    
	  } // Has jet away
	  
	} // Pass trigger
	
      } // Muon-jet selection
      
    } // Loop on events
    
    ThisTree->Close();
	    
  } // Loop on files
  
  cout << "  nEventsFromTrees " << nEventsFromTrees << endl;

  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataRangeName + PUWeighting + KinWeighting + Selection + ".root";

  if (DataType=="Data") OutputFileName.ReplaceAll(PUWeighting + KinWeighting, "");
  
  SaveSystem8Histograms(OutputFileName, DataType);
  
}

void PtRelAnalyzer::SaveSystem8Histograms(TString OutputFileName, TString DataType) {

  cout << "Saving System8 histograms" << endl;

  TFile *OutFile = new TFile(OutputFileName, "recreate");

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      for (int is = 0; is<nSystematics; is++)	{

	if (DataType=="Data") 
	  Observed[ptb][nb][is]->Write();

	for (int tr = 0; tr<nTriggers; tr++) {
	  
	  MuPtForWeighting[tr][is][ptb][nb]->Write();
	  JetPtForWeighting[tr][is][ptb][nb]->Write();
	  JetEtaForWeighting[tr][is][ptb][nb]->Write();
	  PVMultiplicity[tr][is][ptb][nb]->Write();
	  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {
	      
	      if (DataType=="Data") System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Write();
	      else {

		System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Write();
		System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptb][nb]->Write();
		
	      }
	      
	    } 
	  
        } 
	
      }
  
  OutFile->Close();

  cout << "Exiting" << endl;
  
}

void PtRelAnalyzer::MergeSystem8Histograms(TString DataType, int FirstBin, int LastBin) {

  ResetSystem8Histograms(DataType);

  cout << "Merging System8 histograms: " << DataType << endl;
  
  if (LastBin==100) {
    
    if (DataType=="Data") LastBin = nBTagMuRanges - 1;
    else if (DataType=="QCD") LastBin = nMonteCarloPtHatRanges - 1;
    
  }
  
  for (int ph = FirstBin; ph<=LastBin; ph++) {
    
    TString DataRangeName;
    if (DataType=="Data") DataRangeName = BTagMuRangeName[ph];
    else if (DataType=="QCD") DataRangeName = MonteCarloPtHatRange[ph];

    cout << "  Merging " << DataRangeName << endl;
    
    TString InputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataRangeName + PUWeighting + KinWeighting + Selection + ".root";
    if (DataType=="Data") InputFileName.ReplaceAll(PUWeighting + KinWeighting, "");
    TFile *HistogramFile = TFile::Open(InputFileName); 

    TString ThisHistoName;

    for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {

      cout << "    Merging " << PtRelPtBin[ptb] << endl;

      for (int etab = 0; etab<nPtRelEtaBins; etab++) 
	for (int is = 0; is<nSystematics; is++) {

	  if (DataType=="Data") {
	    
	    ThisHistoName = HistogramName("Observed_" + DataType, ptb, etab, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) HistogramFile->Get(ThisHistoName);
	    Observed[ptb][etab][is]->Add(ThisObserved);
	  
	  }

	  for (int tr = 0; tr<nTriggers; tr++) {
	    
	    ThisHistoName = HistogramName("muonPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisMuPt = (TH1D*) HistogramFile->Get(ThisHistoName);
	    MuPtForWeighting[tr][is][ptb][etab]->Add(ThisMuPt);
	    
	    ThisHistoName = HistogramName("jetPt_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetPt = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetPtForWeighting[tr][is][ptb][etab]->Add(ThisJetPt);
	    
	    ThisHistoName = HistogramName("jetEta_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisJetEta = (TH1D*) HistogramFile->Get(ThisHistoName);
	    JetEtaForWeighting[tr][is][ptb][etab]->Add(ThisJetEta);
	    
	    ThisHistoName = HistogramName("nGoodPV_" + DataType, ptb, etab, tr, is, -1, -1);
	    TH1D *ThisPVMult = (TH1D*) HistogramFile->Get(ThisHistoName);
	    PVMultiplicity[tr][is][ptb][etab]->Add(ThisPVMult);
		  
	    for (int st = 0; st<2; st++)
	      for (int tg = 0; tg<nTaggers+1; tg++) {
		
		if (tg==0)
		  ThisHistoName = HistogramName("ntag_pT", ptb, etab, tr, is, -1, -1);
		else ThisHistoName = HistogramName("ntag_pT", ptb, etab, tr, is, tg-1, -1);
		if (st==1) ThisHistoName.ReplaceAll("ntag", "ptag");
		if (tg==0) ThisHistoName.ReplaceAll("tag", "");
		
		if (DataType=="Data") {

		  TH1D *ThisData = (TH1D*) HistogramFile->Get(ThisHistoName);
		  System8ForDataWeighting[st*nTaggers+st+tg][tr][is][ptb][etab]->Add(ThisData);
		
		} else {
		  
		  TH1D *ThisBJet = (TH1D*) HistogramFile->Get(ThisHistoName + "_b");
		  System8ForBJetWeighting[st*nTaggers+st+tg][tr][is][ptb][etab]->Add(ThisBJet);
		  
		  TH1D *ThisLJet = (TH1D*) HistogramFile->Get(ThisHistoName + "_cl");
		  System8ForLJetWeighting[st*nTaggers+st+tg][tr][is][ptb][etab]->Add(ThisLJet);
		  
		}
		
	      }
	    
	  }
      
	}
      
    }

  }
  
  TString OutputFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_" + DataType + PUWeighting + KinWeighting + Selection + ".root";
  if (DataType=="Data") OutputFileName.ReplaceAll(PUWeighting + KinWeighting, "");

  SaveSystem8Histograms(OutputFileName, DataType);

}

void PtRelAnalyzer::BuildSystem8Templates() {

  BookSystem8Templates();

  cout << "Building System8 templates" << endl;
 
  TString DataFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_Data" + Selection + ".root"; DataFileName.ReplaceAll("TrgEmul", "");
  TFile *DataFile = TFile::Open(DataFileName); 
  
  TString QCDFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_QCD" + PUWeighting + KinWeighting + Selection + ".root";
  TFile *QCDFile = TFile::Open(QCDFileName); 

  TString ThisHistoName;
  
  for (int fpt = 0; fpt<nFitPtBins; fpt++) {
    
    cout << "    Merging " << FitPtBin[fpt] << endl;

    float LowFitPtBinEdge = FitPtEdge[fpt];
    float HighFitPtBinEdge = FitPtEdge[fpt+1];

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      float MiddlePtBin = (PtRelPtEdge[bpt] + PtRelPtEdge[bpt+1])/2.;
      
      if (MiddlePtBin>LowFitPtBinEdge && MiddlePtBin<HighFitPtBinEdge) {
	    
	for (int nb = 0; nb<nPtRelEtaBins; nb++) { 
	      
	  for (int is = 0; is<nSystematics; is++) {

	    ThisHistoName = HistogramName("Observed_Data", bpt, nb, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) DataFile->Get(ThisHistoName);
	    
	    for (int tr = 0; tr<nTriggers; tr++)
	      if (AllowedTrigger[bpt][tr]) {

		float LuminosityWeight =  1.;
		
		ThisHistoName = HistogramName("muonPt_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetEta = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataPVMult = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetEta = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDPVMult = (TH1D*) QCDFile->Get(ThisHistoName);
		
		//if (dt==1)
		//if (ThisMuPt->Integral(0, 61)>0.)
		//  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		
		int tr1 = -1, tr2 = -1;
		for (int trg = 0; trg<nTriggers; trg++)
		  if (trg!=tr)
		    if (AllowedTrigger[bpt][trg]) {
		      
		      if (tr1==-1) tr1 = trg;
		      else if (tr2==-1) tr2 = trg;
		      else 
			cout << "PtRelAnalyzer::BuildSystem8Templates Warning:" 
			     << " too many triggers merged in jet pt bin " << PtRelPtBin[bpt] << endl; 
		      
		    }
		
		if (tr1==-1) {
		  
		  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		  
		} else if (tr2==-1) {
		  
		  int T1 = tr, T2 = tr1;
		  if (tr1<tr) { T1 = tr1; T2 = tr; } 
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT2 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) - nT1;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else {
		    
		    cout << "PtRelAnalyzer::BuildSystem8Templates Warning:" 
			 << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << endl; 
		    
		  }
		  
		} else {
		  
		  int T1 = tr, T2 = tr1, T3 = tr2;
		  
		  if (tr<tr1 && tr<tr2 && tr1>tr2) { T1 = tr;  T2 = tr2; T3 = tr1; }
		  if (tr>tr1 && tr<tr2 && tr1<tr2) { T1 = tr1; T2 = tr;  T3 = tr2; }
		  //if (tr<tr1 && tr>tr2 && tr1<tr2) { }
		  if (tr>tr1 && tr>tr2 && tr1<tr2) { T1 = tr1; T2 = tr2; T3 = tr;  }
		  //if (tr>tr1 && tr<tr2 && tr1>tr2) { }
		  if (tr<tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr;  T3 = tr1; }
		  if (tr>tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr1; T3 = tr;  }
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  float nT2 = (ThisObserved->GetBinContent(T1+1, T2+1) - nT1)*(1-(ThisObserved->GetBinContent(T3+1, T3+1)/ThisObserved->GetBinContent(T2+1, T2+1))*(TriggerLuminosity[T2]/TriggerLuminosity[T3]));
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T3, is, -1, -1);
		  TH1D *ThisDataMuPt3 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT3 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) + ThisDataMuPt3->Integral(0, 61) - nT1 - nT2;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else if (tr==T3) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T3, is, -1, -1);
		    TH1D *ThisQCDMuPt3 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT3/ThisQCDMuPt3->Integral(0, 61);
		    
		  } else {
		    
		    cout << "PtRelAnalyzer::BuildSystem8Templates Warning:" 
			   << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << endl; 
		    
		  }
		  
		} 
		
		MuonPt[0][is][fpt][nb]->Add(ThisDataMuPt);
		JetPt[0][is][fpt][nb]->Add(ThisDataJetPt);
		JetEta[0][is][fpt][nb]->Add(ThisDataJetEta);
		PVEvent[0][is][fpt][nb]->Add(ThisDataPVMult);

		if (LuminosityWeight<0. || LuminosityWeight>1000000000.) {

		  cout << "PtRelAnalyzer::BuildSystem8Templates: LuminosityWeight " << LuminosityWeight
		       << " for " << TriggerName[tr] << " in " << PtRelPtBin[bpt] << endl;

		  LuminosityWeight = 0.;

		}

		ThisQCDMuPt->Scale(LuminosityWeight);
		MuonPt[1][is][fpt][nb]->Add(ThisQCDMuPt);

		ThisQCDJetPt->Scale(LuminosityWeight);
		JetPt[1][is][fpt][nb]->Add(ThisQCDJetPt);

		ThisQCDJetEta->Scale(LuminosityWeight);
		JetEta[1][is][fpt][nb]->Add(ThisQCDJetEta);

		ThisQCDPVMult->Scale(LuminosityWeight);
		PVEvent[1][is][fpt][nb]->Add(ThisQCDPVMult);
		
		for (int st = 0; st<2; st++)
		  for (int tg = 0; tg<nTaggers+1; tg++) {
		    
		    if (tg==0)
		      ThisHistoName = HistogramName("ntag_pT", bpt, nb, tr, is, -1, -1);
		    else ThisHistoName = HistogramName("ntag_pT", bpt, nb, tr, is, tg-1, -1);
		    if (st==1) ThisHistoName.ReplaceAll("ntag", "ptag");
		    if (tg==0) ThisHistoName.ReplaceAll("tag", "");
		
		    TH1D *ThisDataSystem8 = (TH1D*) DataFile->Get(ThisHistoName);
		    System8[st*nTaggers+st+tg][is][fpt][nb][0]->Add(ThisDataSystem8);
 
		    TH1D *ThisBJetSystem8 = (TH1D*) QCDFile->Get(ThisHistoName + "_b");
		    ThisBJetSystem8->Scale(LuminosityWeight);
		    System8[st*nTaggers+st+tg][is][fpt][nb][1]->Add(ThisBJetSystem8);
 
		    TH1D *ThisLJetSystem8 = (TH1D*) QCDFile->Get(ThisHistoName + "_cl");
		    ThisLJetSystem8->Scale(LuminosityWeight);
		    System8[st*nTaggers+st+tg][is][fpt][nb][2]->Add(ThisLJetSystem8);
 
		  }		
		
	      }
	    
	  }
	  
	}
	
      }
      
    }
    
  }

  TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + PUWeighting + KinWeighting + Selection + Production + ".root";
  
  SaveSystem8Templates(TemplateFileName);

}
  
void PtRelAnalyzer::BookSystem8Templates() {
  
  cout << "Booking System8 templates" << endl;

  TString ThisHistoName;
  
  for (int fpt = 0; fpt<nFitPtBins; fpt++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)       
      for (int is = 0; is<nSystematics; is++) 
	for (int dt = 0; dt<2; dt++) {

	  TString DataType = "Data";
	  if (dt==1) DataType = "QCD";
	  
	  ThisHistoName =  ThisHistoName = HistogramName("jetPt_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  JetPt[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 600, 0., 600.); 
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  JetEta[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.); 

	  ThisHistoName = HistogramName("muonPt_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  MuonPt[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, 0., 60.); 
	  
	  ThisHistoName = HistogramName("PV_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  PVEvent[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 35, 0., 35.);
	  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {

	      if (tg==0)
		ThisHistoName = HistogramName("ntag_pT", FitPtBin[fpt], nb, -1, is, -1, -1);
	      else ThisHistoName = HistogramName("ntag_pT", FitPtBin[fpt], nb, -1, is, tg-1, -1);
	      if (st==1) ThisHistoName.ReplaceAll("ntag", "ptag");
	      if (tg==0) ThisHistoName.ReplaceAll("tag", "");
	      
	      if (dt==0) System8[st*nTaggers+st+tg][is][fpt][nb][0] = new TH1D(ThisHistoName,         ThisHistoName,         nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
	      else {
		System8[st*nTaggers+st+tg][is][fpt][nb][1] = new TH1D(ThisHistoName + "_b",  ThisHistoName + "_b",  nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
		System8[st*nTaggers+st+tg][is][fpt][nb][2] = new TH1D(ThisHistoName + "_cl", ThisHistoName + "_cl", nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp);
	      }

	    }
	  
	}
  
}

void PtRelAnalyzer::SaveSystem8Templates(TString TemplateFileName) {
  
  cout << "Saving System8 templates" << endl;

  TFile *OutFile = new TFile(TemplateFileName, "recreate");

  TString ThisHistoName;
  
  for (int fpt = 0; fpt<nFitPtBins; fpt++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)       
      for (int is = 0; is<nSystematics; is++) 
	for (int dt = 0; dt<2; dt++) {
	  	  
	  JetPt[dt][is][fpt][nb]->Write();
	  JetEta[dt][is][fpt][nb]->Write();
	  MuonPt[dt][is][fpt][nb]->Write();
	  PVEvent[dt][is][fpt][nb]->Write();
	  
	  for (int st = 0; st<2; st++)
	    for (int tg = 0; tg<nTaggers+1; tg++) {
	      
	      if (dt==0) System8[st*nTaggers+st+tg][is][fpt][nb][0]->Write();
	      else {
		System8[st*nTaggers+st+tg][is][fpt][nb][1]->Write();
		System8[st*nTaggers+st+tg][is][fpt][nb][2]->Write();
	      }

	    }
	  
	}
  
  OutFile->Close();
  
  cout << "Exiting" << endl;
  
}

int PtRelAnalyzer::GetAwayJet(TString ThisAwayTaggerName, int jMu, float AwayDRCut, bool IsUnique) {
  
  int aJet = -1;

  int nAwayJets = 0;
  float PtLeadingAwayJet = -1;
    
  for (int ijet = 0; ijet<nJet; ijet++) 
    if (Jet_pt[ijet]>PtRelPtEdge[0] && ijet!=jMu && fabs(Jet_eta[ijet])<2.4) {
    
      if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], Jet_eta[jMu], Jet_phi[jMu])>AwayDRCut) {
	
	bool pass = IsTaggedJet(ijet, ThisAwayTaggerName);
		
	if (pass) {
	
	  if (Jet_pt[ijet]>PtLeadingAwayJet) {
	    
	    aJet = ijet;
	    PtLeadingAwayJet = Jet_pt[ijet];

	  }

	  nAwayJets++;

	}
      
      }

    }

  if (IsUnique && nAwayJets>1) aJet = -1;

  return aJet;
   
}

void PtRelAnalyzer::GetDRCuts(TString ThisSystematicName, float MuonJetPt, float *TrackDRCut, float *TrackMinDRCut) {

  *TrackDRCut    = 999.; 
  *TrackMinDRCut = 0.;
  
  if (!ThisSystematicName.Contains("_MuDR")) {

    if (TemplateVariable=="PtRel") {

      if (MuonJetPt<30) *TrackDRCut    = 0.20; 
      else if (MuonJetPt<80) *TrackDRCut    = 0.15; 
      else *TrackDRCut    = 0.12;
    
    } else if (TemplateVariable=="System8") {

      *TrackDRCut    = 0.40;
    
    } else {

      *TrackDRCut    = 999.; 
      *TrackMinDRCut = 0.;
      
    } 

  } else if (ThisSystematicName=="_MuDRMinus") {
    
    if (TemplateVariable=="PtRel") {

      if (MuonJetPt<30) *TrackDRCut    = 0.15; 
      else if (MuonJetPt<80) *TrackDRCut    = 0.12; 
      else *TrackDRCut    = 0.09; 

    } else if (TemplateVariable=="System8") {

      *TrackDRCut    = 0.30;
    
    } else {
    
      if (MuonJetPt<30) *TrackDRCut    = 0.20; 
      else if (MuonJetPt<80) *TrackDRCut    = 0.15; 
      else *TrackDRCut    = 0.12; 
      
    }
    
  } else if (ThisSystematicName=="_MuDRPlus") {
    
    *TrackDRCut    = 999.; 
    *TrackMinDRCut = 0.;
    
  } else {
    
    if (ThisSystematicName=="_MuDR003") *TrackDRCut = 0.03;
    else if (ThisSystematicName=="_MuDR004") *TrackDRCut = 0.04;
    else if (ThisSystematicName=="_MuDR005") *TrackDRCut = 0.05;
    else if (ThisSystematicName=="_MuDR006") *TrackDRCut = 0.06;
    else if (ThisSystematicName=="_MuDR007") *TrackDRCut = 0.07;
    else if (ThisSystematicName=="_MuDR008") *TrackDRCut = 0.08;
    else if (ThisSystematicName=="_MuDR009") *TrackDRCut = 0.09;
    else if (ThisSystematicName=="_MuDR010") *TrackDRCut = 0.10;
    else if (ThisSystematicName=="_MuDR012") *TrackDRCut = 0.12;
    else if (ThisSystematicName=="_MuDR015") *TrackDRCut = 0.15;
    else if (ThisSystematicName=="_MuDR020") *TrackDRCut = 0.20;
    else if (ThisSystematicName=="_MuDR050") *TrackDRCut = 999.;
    else if (ThisSystematicName=="_MuDRp08") *TrackMinDRCut = 0.08;
    else if (ThisSystematicName=="_MuDRp10") *TrackMinDRCut = 0.10; 
    else if (ThisSystematicName=="_MuDRp15") *TrackMinDRCut = 0.15; 
    else std::cout << "Wrong choise for MuDRCut" << std::endl;
    
  }
  
}
 
void PtRelAnalyzer::GetPileUpWeights(TString DataType) { 

  if (DataType=="Data") {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int pu = 0; pu<nMaxPU; pu++)
	for (int is = 0; is<nSystematics; is++)
	  PileUpWeight[tr][pu][is] = 1.;
      
  } else if (PUWeighting.Contains("PUWeighting")) {
    /*
    for (int tr = 0; tr<nTriggers; tr++)
      if(TriggerCode.Contains(TriggerName[tr])) 
	for (int is = 0; is<nSystematics; is++) {
	
	  TString WeightsFileName = "./PileUp/LowPt_DiJet20_" + DataForPileUpWeights + ".txt";

	  WeightsFileName.ReplaceAll("_DiJet20", TriggerName[tr]);
	  
	  if (Type==0) WeightsFileName.ReplaceAll("LowPt", "Jet");
	  if (Type==1) WeightsFileName.ReplaceAll("LowPt", "LightQCD");
	  if (TriggerName[tr].Contains("_Jet")) WeightsFileName.ReplaceAll("LowPt", "HighPt");
	  
	  if (is==1) WeightsFileName.ReplaceAll(".txt", "_m05.txt");
	  if (is==2) WeightsFileName.ReplaceAll(".txt", "_p05.txt");
	  
	  ifstream WeightsFile; WeightsFile.open(WeightsFileName);
	  
	  while (WeightsFile) {
	  
	    int ThisPU; float ThisWeight;
	    WeightsFile >> ThisPU >> ThisWeight;
	    
	    PileUpWeight[tr][ThisPU][is] = ThisWeight;
	    
	  }
	
	}
    */
  } else if (PUWeighting.Contains("_PV")) {
    
    TString PVFileFlag = PUWeighting; PVFileFlag.ReplaceAll("PV", "");

    for (int tr = 0; tr<nTriggers; tr++) {
      
      ifstream PVWeightFile; PVWeightFile.open("./Weights/PileUp/PVMultTriggered" + TriggerName[tr] + "_" + DataType + PVFileFlag + ".txt");

      for (int npv = 0; npv<nMaxPU; npv++) {

	int pvm; float weightpv;
	PVWeightFile >> pvm >> weightpv;

	for (int is = 0; is<nSystematics; is++)  
	  PileUpWeight[tr][pvm][is] = weightpv;

      }
      
    }

  } else {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int pu = 0; pu<nMaxPU; pu++)
	for (int is = 0; is<nSystematics; is++)
	  PileUpWeight[tr][pu][is] = 1.;
    
  }
  
}

void PtRelAnalyzer::GetKinematicWeights(TString DataType) {

  cout << "Getting kinematic weights" << endl;

  if (DataType=="Data") {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<600; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 1; nb<nPtRelEtaBins; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb-1] = 1.;
    
  } else if (KinWeighting.Contains("_KinPtEtaBins")) {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<600; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 0; nb<nPtRelEtaBins-1; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb] = 1.;

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {

      //for (int nb = 0; nb<nPtRelEtaBins; nb++) {
      
        TString CorrectionFileName = "./Weights/KinematicWeights/" + TemplateVariable + KinWeighting + "_Pt_" + PtRelPtBin[bpt] + "_anyEta_" + DataType + PUWeighting + Selection + ".txt";

	if (DataType=="QCD") CorrectionFileName.ReplaceAll("IP3D_", "IP3DSig_");
	CorrectionFileName.ReplaceAll("SoftMu_", "IP3DSig_");
	ifstream PtWeightsFile; PtWeightsFile.open(CorrectionFileName);
	
	int nPtSteps = PtRelPtEdge[bpt+1] - PtRelPtEdge[bpt];
	for (int stpt = 0; stpt<nPtSteps; stpt++) {
	  
	  float PtStep, ThisPtWeight; 
	  PtWeightsFile >> PtStep >> ThisPtWeight;	 

	  for (int tr = 0; tr<nTriggers; tr++)
	    for (int is = 0; is<nSystematics; is++) 
	      for (int imu = 0; imu<60; imu++) {

		KinematicWeight[tr][is][int(PtStep)][imu][0] = ThisPtWeight;
		KinematicWeight[tr][is][int(PtStep)][imu][1] = ThisPtWeight;
		
	      }  
	  
	}
	      
	//ifstream EtaWeightsFile; EtaWeightsFile.open("./Weights/KinematicWeights/" + KinWeighting + "_Eta_" + PtRelPtBin[PtBin] + "_anyEta" + PUWeighting + Selection + ".txt");
    
	//int const nEtaSteps = 48;

	//}
      
    }

  } else if (KinWeighting.Contains("_KinPtEta2Bins")) {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<600; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 0; nb<nPtRelEtaBins-1; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb] = 1.;

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      for (int nb = 1; nb<nPtRelEtaBins; nb++) {
	
        TString CorrectionFileName = "./Weights/KinematicWeights/" + TemplateVariable + KinWeighting + "_Pt_" + PtRelPtBin[bpt] + "_" + PtRelEtaBin[nb] + "_" + DataType + PUWeighting + Selection + ".txt";
	CorrectionFileName.ReplaceAll("2Bins", "Bins");
	cout << CorrectionFileName << endl;
	CorrectionFileName.ReplaceAll("IP3D_", "IP3DSig_");
	CorrectionFileName.ReplaceAll("SoftMu_", "IP3DSig_");

	ifstream PtWeightsFile; PtWeightsFile.open(CorrectionFileName);
	
	int nPtSteps = PtRelPtEdge[bpt+1] - PtRelPtEdge[bpt];
	for (int stpt = 0; stpt<nPtSteps; stpt++) {
	  
	  float PtStep, ThisPtWeight; 
	  PtWeightsFile >> PtStep >> ThisPtWeight;	 
	  cout << " " << PtStep << " " << ThisPtWeight << endl;
	  for (int tr = 0; tr<nTriggers; tr++)
	    for (int is = 0; is<nSystematics; is++) 
	      for (int imu = 0; imu<60; imu++) {
		
		KinematicWeight[tr][is][int(PtStep)][imu][nb-1] = ThisPtWeight;
		
	      }  
	  
	}
	      
	//ifstream EtaWeightsFile; EtaWeightsFile.open("./Weights/KinematicWeights/" + KinWeighting + "_Eta_" + PtRelPtBin[PtBin] + "_anyEta" + PUWeighting + Selection + ".txt");
	
	//int const nEtaSteps = 48;
	
      }
      
    }

  } else if (KinWeighting.Contains("_KinWeights") || KinWeighting.Contains("_FullKinWeights")) {
    /*
    for (int tr = 0; tr<nTriggers; tr++)
      for (int npt = 0; npt<nKinWeightsPtBins; npt++) 
	if( AllowedTrigger[npt].Contains(TriggerName[tr])) 
	  for (int is = 0; is<nSystematics; is++) {
	  
	    TString ThisAwayTagger = AwayTaggerName[0];
	    for (int at = 0; at<nAwayTaggers; at++) 
	      if (SystematicName[is].Contains(AwayTaggerName[at]))
		ThisAwayTagger = AwayTaggerName[at];

	    TString WeightsFileName = "./KinematicWeights/KinWeights-" + ThisAwayTagger; 
	    if (KinWeighting.Contains("_FullKinW")) WeightsFileName.ReplaceAll("Weights/Kin", "Weights/FullKin");
	    WeightsFileName += "-";
	    WeightsFileName += TriggerName[tr];
	    WeightsFileName.ReplaceAll("-_", "-"); 
	    WeightsFileName += "-";
	    WeightsFileName += KinWeightsPtBin[npt];
	  
	    if (Type==0) WeightsFileName += "_LightTemplatesJet";
	    if (Type==1) WeightsFileName += "_LightTemplatesQCD";
	  
	    WeightsFileName += ".txt"; 
	  
	    ifstream WeightsFile; WeightsFile.open(WeightsFileName);
	  
	    while (WeightsFile) {
	    
	      int ipt, nb, imu; double ThisWeight;
	      WeightsFile >> ipt >> nb >> imu >> ThisWeight;
	      
	      KinematicWeight[tr][is][ipt][imu][nb] = ThisWeight;
	      
	    }

	  }
    */
  } else {
    
    for (int tr = 0; tr<nTriggers; tr++)
      for (int is = 0; is<nSystematics; is++) 
	for (int ipt = 0; ipt<600; ipt++) 
	  for (int imu = 0; imu<60; imu++)
	    for (int nb = 0; nb<nPtRelEtaBins-1; nb++) 
	      KinematicWeight[tr][is][ipt][imu][nb] = 1.;
    
  }
  
}

void PtRelAnalyzer::GetBTemplateCorrections() {

  for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {

    for (int ib = 0; ib<100; ib++) {

      BTemplateCorrections[ib][ptb][0] = 1.;
      BTemplateCorrections[ib][ptb][1] = 1.;

    }

    TString ThisPtRelPtBin = PtRelPtBin[ptb];
    if (ThisPtRelPtBin=="Pt120140" || ThisPtRelPtBin=="Pt140160") ThisPtRelPtBin = "Pt120160";
    if (ThisPtRelPtBin=="Pt260300" || ThisPtRelPtBin=="Pt300320") ThisPtRelPtBin = "Pt260320";
    if (ThisPtRelPtBin=="Pt500600" || ThisPtRelPtBin=="Pt600")    ThisPtRelPtBin = "Pt500";

    ifstream MnusCorrectionsFile; MnusCorrectionsFile.open("./Weights/BTemplatesCorrections/EnergyFraction_" + ThisPtRelPtBin + "_m5.txt");

    while (MnusCorrectionsFile) {

      float xBin, efcorr;
      MnusCorrectionsFile >> xBin >> efcorr;
      
      if (efcorr>4.) efcorr = 1.;

      int ib = xBin/0.02;
      BTemplateCorrections[ib][ptb][0] = efcorr;

    }

    ifstream PlusCorrectionsFile; PlusCorrectionsFile.open("./Weights/BTemplatesCorrections/EnergyFraction_" + ThisPtRelPtBin + "_p5.txt");

    while (PlusCorrectionsFile) {

      float xBin, efcorr;
      PlusCorrectionsFile >> xBin >> efcorr;
      
      if (efcorr>4.) efcorr = 1.;

      int ib = xBin/0.02;
      BTemplateCorrections[ib][ptb][1] = efcorr;

    }

  }

}

void PtRelAnalyzer::BookTemplates() {

  cout << "Booking templates" << endl;

  TString ThisHistoName;

  for (int fpt = 0; fpt<nFitPtBins; fpt++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++)       
      for (int is = 0; is<nSystematics; is++) 
	for (int dt = 0; dt<4; dt++) {

	  TString DataType = "Data";
	  if (dt==1) DataType = "QCD";
	  if (dt==2) DataType = "Jet";
	  if (dt==3) DataType = "Inc";
	  
	  ThisHistoName =  ThisHistoName = HistogramName("jetPt_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  JetPt[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 600, 0., 600.); 
	  
	  ThisHistoName = HistogramName("jetEta_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	  JetEta[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, -3., 3.); 

	  if (dt<2) {

	    ThisHistoName = HistogramName("muonPt_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	    MuonPt[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 60, 0., 60.); 
	    
	    ThisHistoName = HistogramName("muonDR_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	    MuonDR[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 50, 0., 0.5);
	  
	    ThisHistoName = HistogramName("PV_" + DataType, FitPtBin[fpt], nb, -1, is, -1, -1);
	    PVEvent[dt][is][fpt][nb] = new TH1D(ThisHistoName, ThisHistoName, 35, 0., 35.);

	    for (int tg = 0; tg<nTaggers; tg++)
	      for (int tp = 0; tp<2; tp++)
		for (int fl = 0; fl<4; fl++) {
		  
		  if (dt==0 && fl!=0 && fl!=3) continue;
  
		  int lf = 0;
		  if (dt==0 && fl==3) lf = 4;
		  if (dt==1 && fl==0) lf = 1;
		  if (dt==1 && fl==1) lf = 2;
		  if (dt==1 && fl==2) lf = 3;
		  if (dt==1 && fl==3) lf = 4;
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_" + DataType, FitPtBin[fpt], nb, -1, is, tg, 10*tp + lf);
		  PtRel[dt*nTaggers+tg][is][fpt][nb][4*tp+fl] = new TH1D(ThisHistoName, ThisHistoName, nBinsForTemp, LowerEdgeForTemp, UpperEdgeForTemp); 
		  
		}

	  }

	}
  
}

void PtRelAnalyzer::BuildTemplates(bool AddLightTemplates) {
  
  BookTemplates();

  cout << "Building templates" << endl;
  
  TFile *DataFile, *QCDFile, *JetFile, *IncFile;

  TString DataFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_Data" + Selection + ".root"; DataFileName.ReplaceAll("TrgEmul", "");
  //DataFileName.ReplaceAll("JEC", "");
  DataFile = TFile::Open(DataFileName); 
  
  TString QCDFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_QCD" + PUWeighting + KinWeighting + Selection + ".root";
  QCDFile = TFile::Open(QCDFileName); 

  if (AddLightTemplates) {

    TString JetFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_Jet" + PUWeighting + KinWeighting + Selection + ".root"; JetFileName.ReplaceAll("TrgEmul", "");
    JetFile = TFile::Open(JetFileName); 
    
    TString IncFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_Inc" + PUWeighting + KinWeighting + Selection + ".root"; IncFileName.ReplaceAll("TrgEmul", "");
    IncFile = TFile::Open(IncFileName); 

  }

  TString ThisHistoName;

  for (int fpt = 0; fpt<nFitPtBins; fpt++) {
    
    cout << "    Merging " << FitPtBin[fpt] << endl;

    float LowFitPtBinEdge = FitPtEdge[fpt];
    float HighFitPtBinEdge = FitPtEdge[fpt+1];

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      float MiddlePtBin = (PtRelPtEdge[bpt] + PtRelPtEdge[bpt+1])/2.;
      
      if (MiddlePtBin>LowFitPtBinEdge && MiddlePtBin<HighFitPtBinEdge) {
	    
	for (int nb = 0; nb<nPtRelEtaBins; nb++) { 
	      
	  for (int is = 0; is<nSystematics; is++) {

	    ThisHistoName = HistogramName("Observed_Data", bpt, nb, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) DataFile->Get(ThisHistoName);
	    
	    for (int tr = 0; tr<nTriggers; tr++)
	      if (AllowedTrigger[bpt][tr]) {
		
		float LuminosityWeight =  1.;
		
		ThisHistoName = HistogramName("muonPt_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonDR_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuDR = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetEta = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataPVMult = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonDR_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuDR = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetEta = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDPVMult = (TH1D*) QCDFile->Get(ThisHistoName);
		
		//if (dt==1)
		//if (ThisMuPt->Integral(0, 61)>0.)
		//  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		
		int tr1 = -1, tr2 = -1;
		for (int trg = 0; trg<nTriggers; trg++)
		  if (trg!=tr)
		    if (AllowedTrigger[bpt][trg]) {
		      
		      if (tr1==-1) tr1 = trg;
		      else if (tr2==-1) tr2 = trg;
		      else 
			cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			     << " too many triggers merged in jet pt bin " << PtRelPtBin[bpt] << endl; 
		      
		    }
		
		if (tr1==-1) {
		  
		  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		  
		} else if (tr2==-1) {
		  
		  int T1 = tr, T2 = tr1;
		  if (tr1<tr) { T1 = tr1; T2 = tr; } 
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT2 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) - nT1;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else {
		    
		    cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			 << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << endl; 
		    
		  }
		  
		} else {
		  
		  int T1 = tr, T2 = tr1, T3 = tr2;
		  
		  if (tr<tr1 && tr<tr2 && tr1>tr2) { T1 = tr;  T2 = tr2; T3 = tr1; }
		  if (tr>tr1 && tr<tr2 && tr1<tr2) { T1 = tr1; T2 = tr;  T3 = tr2; }
		  //if (tr<tr1 && tr>tr2 && tr1<tr2) { }
		  if (tr>tr1 && tr>tr2 && tr1<tr2) { T1 = tr1; T2 = tr2; T3 = tr;  }
		  //if (tr>tr1 && tr<tr2 && tr1>tr2) { }
		  if (tr<tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr;  T3 = tr1; }
		  if (tr>tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr1; T3 = tr;  }
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  float nT2 = (ThisObserved->GetBinContent(T1+1, T2+1) - nT1)*(1-(ThisObserved->GetBinContent(T3+1, T3+1)/ThisObserved->GetBinContent(T2+1, T2+1))*(TriggerLuminosity[T2]/TriggerLuminosity[T3]));
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T3, is, -1, -1);
		  TH1D *ThisDataMuPt3 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT3 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) + ThisDataMuPt3->Integral(0, 61) - nT1 - nT2;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else if (tr==T3) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T3, is, -1, -1);
		    TH1D *ThisQCDMuPt3 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT3/ThisQCDMuPt3->Integral(0, 61);
		    
		  } else {
		    
		    cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			   << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << endl; 
		    
		  }
		  
		} 
		
		MuonPt[0][is][fpt][nb]->Add(ThisDataMuPt);
		MuonDR[0][is][fpt][nb]->Add(ThisDataMuDR);
		JetPt[0][is][fpt][nb]->Add(ThisDataJetPt);
		JetEta[0][is][fpt][nb]->Add(ThisDataJetEta);
		PVEvent[0][is][fpt][nb]->Add(ThisDataPVMult);

		if (LuminosityWeight<0. || LuminosityWeight>1000000000.) {

		  cout << "PtRelAnalyzer::BuildTemplates: LuminosityWeight " << LuminosityWeight
		       << " for " << TriggerName[tr] << " in " << PtRelPtBin[bpt] << endl;

		  LuminosityWeight = 0.;

		}

		ThisQCDMuPt->Scale(LuminosityWeight);
		MuonPt[1][is][fpt][nb]->Add(ThisQCDMuPt);

		ThisQCDMuDR->Scale(LuminosityWeight);
		MuonDR[1][is][fpt][nb]->Add(ThisQCDMuDR);

		ThisQCDJetPt->Scale(LuminosityWeight);
		JetPt[1][is][fpt][nb]->Add(ThisQCDJetPt);

		ThisQCDJetEta->Scale(LuminosityWeight);
		JetEta[1][is][fpt][nb]->Add(ThisQCDJetEta);

		ThisQCDPVMult->Scale(LuminosityWeight);
		PVEvent[1][is][fpt][nb]->Add(ThisQCDPVMult);
 
		if (AddLightTemplates) {
		 
		  int nScaleBins = ThisQCDJetPt->GetNbinsX();

		  ThisHistoName = HistogramName("jetPt_Jet", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisJetJetPt = (TH1D*) JetFile->Get(ThisHistoName);

                  float JetLuminosityWeight = ThisQCDJetPt->Integral(0, nScaleBins)/ThisJetJetPt->Integral(0, nScaleBins);

                  ThisJetJetPt->Scale(JetLuminosityWeight);
		  JetPt[2][is][fpt][nb]->Add(ThisJetJetPt);
		  
		  ThisHistoName = HistogramName("jetEta_Jet", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisJetJetEta = (TH1D*) JetFile->Get(ThisHistoName);

                  ThisJetJetEta->Scale(JetLuminosityWeight);
		  JetEta[2][is][fpt][nb]->Add(ThisJetJetEta);
		  
		  ThisHistoName = HistogramName("jetPt_Inc", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisIncJetPt = (TH1D*) IncFile->Get(ThisHistoName);
		  
                  float IncLuminosityWeight = ThisQCDJetPt->Integral(0, nScaleBins)/ThisIncJetPt->Integral(0, nScaleBins);

                  ThisIncJetPt->Scale(IncLuminosityWeight);
		  JetPt[3][is][fpt][nb]->Add(ThisIncJetPt);

		  ThisHistoName = HistogramName("jetEta_Inc", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisIncJetEta = (TH1D*) IncFile->Get(ThisHistoName);

                  ThisIncJetEta->Scale(IncLuminosityWeight);
		  JetEta[3][is][fpt][nb]->Add(ThisIncJetEta);
		  
		}

		for (int tg = 0; tg<nTaggers; tg++) {
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_Data", bpt, nb, tr, is, tg,  0);
		  TH1D *ThisDataPtRelTag = (TH1D*) DataFile->Get(ThisHistoName);
		  PtRel[tg][is][fpt][nb][0]->Add(ThisDataPtRelTag);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_Data", bpt, nb, tr, is, tg, 10);
		  TH1D *ThisDataPtRelUntag = (TH1D*) DataFile->Get(ThisHistoName);
		  PtRel[tg][is][fpt][nb][4]->Add(ThisDataPtRelUntag);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg,  1);
		  TH1D *ThisQCDPtRelTagB = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagB->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][0]->Add(ThisQCDPtRelTagB);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg, 11);
		  TH1D *ThisQCDPtRelUntagB = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagB->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][4]->Add(ThisQCDPtRelUntagB);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg,  2);
		  TH1D *ThisQCDPtRelTagC = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagC->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][1]->Add(ThisQCDPtRelTagC);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg, 12);
		  TH1D *ThisQCDPtRelUntagC = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagC->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][5]->Add(ThisQCDPtRelUntagC);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg,  3);
		  TH1D *ThisQCDPtRelTagLG = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagLG->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][2]->Add(ThisQCDPtRelTagLG);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg, 13);
		  TH1D *ThisQCDPtRelUntagLG = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagLG->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][6]->Add(ThisQCDPtRelUntagLG);

		  if (AddLightTemplates) {

		    int nScaleBins = ThisQCDPtRelTagLG->GetNbinsX();
		    float TagLuminosityWeight   = ThisQCDPtRelTagLG->Integral(0, nScaleBins);
		    float UntagLuminosityWeight = ThisQCDPtRelUntagLG->Integral(0, nScaleBins);

		    ThisHistoName = HistogramName(TemplateVariable + "_Jet", bpt, nb, tr, is, tg,  4);
		    TH1D *ThisJetPtRelTagTrk = (TH1D*) JetFile->Get(ThisHistoName);
                    float JetTagLuminosityWeight = TagLuminosityWeight/ThisJetPtRelTagTrk->Integral(0, nScaleBins);
		    ThisJetPtRelTagTrk->Scale(JetTagLuminosityWeight);
		    PtRel[tg][is][fpt][nb][3]->Add(ThisJetPtRelTagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_Jet", bpt, nb, tr, is, tg, 14);
		    TH1D *ThisJetPtRelUntagTrk = (TH1D*) JetFile->Get(ThisHistoName);
                    float JetUntagLuminosityWeight = UntagLuminosityWeight/ThisJetPtRelUntagTrk->Integral(0, nScaleBins);
		    ThisJetPtRelUntagTrk->Scale(JetUntagLuminosityWeight);
		    PtRel[tg][is][fpt][nb][7]->Add(ThisJetPtRelUntagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_Inc", bpt, nb, tr, is, tg,  4);
		    TH1D *ThisIncPtRelTagTrk = (TH1D*) IncFile->Get(ThisHistoName);
                    float IncTagLuminosityWeight = TagLuminosityWeight/ThisIncPtRelTagTrk->Integral(0, nScaleBins);
		    ThisIncPtRelTagTrk->Scale(IncTagLuminosityWeight);
		    PtRel[nTaggers+tg][is][fpt][nb][3]->Add(ThisIncPtRelTagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_Inc", bpt, nb, tr, is, tg, 14);
		    TH1D *ThisIncPtRelUntagTrk = (TH1D*) IncFile->Get(ThisHistoName);
                    float IncUntagLuminosityWeight = UntagLuminosityWeight/ThisIncPtRelUntagTrk->Integral(0, nScaleBins);
		    ThisIncPtRelUntagTrk->Scale(IncUntagLuminosityWeight);
		    PtRel[nTaggers+tg][is][fpt][nb][7]->Add(ThisIncPtRelUntagTrk);

		  }

		}

	      }
		      
	    }
		    
	  }
		  
      }
		
    }
    
  }
    
  TString LightTemplates = "";
  if (AddLightTemplates) LightTemplates = "All";	
  TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + LightTemplates + PUWeighting + KinWeighting + Selection + Production + ".root";
  
  SaveTemplates(TemplateFileName, AddLightTemplates);
  
}

void PtRelAnalyzer::BuildTemplatesEtaBins(bool AddLightTemplates) {
  
  BookTemplates();

  cout << "Building templates" << endl;
  
  TFile *DataFile, *QCDFile, *JetFile, *IncFile;

  TString DataFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_Data" + Selection + ".root"; DataFileName.ReplaceAll("TrgEmul", "");
  DataFile = TFile::Open(DataFileName); 
  
  TString QCDFileName = "./Templates/Histograms/" + TemplateVariable + "_Histograms_QCD" + PUWeighting + KinWeighting + Selection + ".root";
  QCDFile = TFile::Open(QCDFileName); 

  if (AddLightTemplates) {

    TString JetFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_Jet" + PUWeighting + KinWeighting + Selection + ".root"; JetFileName.ReplaceAll("TrgEmul", "");
    JetFile = TFile::Open(JetFileName); 
    
    TString IncFileName = "./Templates/Histograms/" + TemplateVariable + "_LightHistograms_Inc" + PUWeighting + KinWeighting + Selection + ".root"; IncFileName.ReplaceAll("TrgEmul", "");
    IncFile = TFile::Open(IncFileName); 

  }

  TString ThisHistoName;

  for (int fpt = 0; fpt<nFitPtBins; fpt++) {
    
    cout << "    Merging " << FitPtBin[fpt] << endl;

    float LowFitPtBinEdge = FitPtEdge[fpt];
    float HighFitPtBinEdge = FitPtEdge[fpt+1];

    for (int bpt = 0; bpt<nPtRelPtBins; bpt++) {
      
      float MiddlePtBin = (PtRelPtEdge[bpt] + PtRelPtEdge[bpt+1])/2.;
      
      if (MiddlePtBin>LowFitPtBinEdge && MiddlePtBin<HighFitPtBinEdge) {
	    
	for (int nb = 1; nb<nPtRelEtaBins; nb++) { 
	      
	  for (int is = 0; is<nSystematics; is++) {

	    ThisHistoName = HistogramName("Observed_Data", bpt, nb, -1, is, -1, -1);
	    TH2D *ThisObserved = (TH2D*) DataFile->Get(ThisHistoName);
	    
	    for (int tr = 0; tr<nTriggers; tr++)
	      if (AllowedTrigger[bpt][tr]) {
		
		float LuminosityWeight =  1.;
		
		ThisHistoName = HistogramName("muonPt_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonDR_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataMuDR = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetPt = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataJetEta = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_Data", bpt, nb, tr, is, -1, -1);
		TH1D *ThisDataPVMult = (TH1D*) DataFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("muonDR_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDMuDR = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetPt_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetPt = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("jetEta_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDJetEta = (TH1D*) QCDFile->Get(ThisHistoName);
		
		ThisHistoName = HistogramName("nGoodPV_QCD", bpt, nb, tr, is, -1, -1);
		TH1D *ThisQCDPVMult = (TH1D*) QCDFile->Get(ThisHistoName);
		
		//if (dt==1)
		//if (ThisMuPt->Integral(0, 61)>0.)
		//  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		
		int tr1 = -1, tr2 = -1;
		for (int trg = 0; trg<nTriggers; trg++)
		  if (trg!=tr)
		    if (AllowedTrigger[bpt][trg]) {
		      
		      if (tr1==-1) tr1 = trg;
		      else if (tr2==-1) tr2 = trg;
		      else 
			cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			     << " too many triggers merged in jet pt bin " << PtRelPtBin[bpt] << endl; 
		      
		    }
		
		if (tr1==-1) {
		  
		  LuminosityWeight = ThisDataMuPt->Integral(0, 61)/ThisQCDMuPt->Integral(0, 61);
		  
		} else if (tr2==-1) {
		  
		  int T1 = tr, T2 = tr1;
		  if (tr1<tr) { T1 = tr1; T2 = tr; } 
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT2 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) - nT1;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else {
		    
		    cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			 << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << endl; 
		    
		  }
		  
		} else {
		  
		  int T1 = tr, T2 = tr1, T3 = tr2;
		  
		  if (tr<tr1 && tr<tr2 && tr1>tr2) { T1 = tr;  T2 = tr2; T3 = tr1; }
		  if (tr>tr1 && tr<tr2 && tr1<tr2) { T1 = tr1; T2 = tr;  T3 = tr2; }
		  //if (tr<tr1 && tr>tr2 && tr1<tr2) { }
		  if (tr>tr1 && tr>tr2 && tr1<tr2) { T1 = tr1; T2 = tr2; T3 = tr;  }
		  //if (tr>tr1 && tr<tr2 && tr1>tr2) { }
		  if (tr<tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr;  T3 = tr1; }
		  if (tr>tr1 && tr>tr2 && tr1>tr2) { T1 = tr2; T2 = tr1; T3 = tr;  }
		  
		  float nT1 = ThisObserved->GetBinContent(T1+1, T1+1) - ThisObserved->GetBinContent(T2+1, T2+1)*TriggerLuminosity[T1]/TriggerLuminosity[T2];
		  
		  float nT2 = (ThisObserved->GetBinContent(T1+1, T2+1) - nT1)*(1-(ThisObserved->GetBinContent(T3+1, T3+1)/ThisObserved->GetBinContent(T2+1, T2+1))*(TriggerLuminosity[T2]/TriggerLuminosity[T3]));
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T1, is, -1, -1);
		  TH1D *ThisDataMuPt1 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T2, is, -1, -1);
		  TH1D *ThisDataMuPt2 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  ThisHistoName = HistogramName("muonPt_Data", bpt, nb, T3, is, -1, -1);
		  TH1D *ThisDataMuPt3 = (TH1D*) DataFile->Get(ThisHistoName);
		  
		  float nT3 = ThisDataMuPt1->Integral(0, 61) + ThisDataMuPt2->Integral(0, 61) + ThisDataMuPt3->Integral(0, 61) - nT1 - nT2;
		  
		  if (tr==T1) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T1, is, -1, -1);
		    TH1D *ThisQCDMuPt1 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT1/ThisQCDMuPt1->Integral(0, 61);
		    
		  } else if (tr==T2) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T2, is, -1, -1);
		    TH1D *ThisQCDMuPt2 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT2/ThisQCDMuPt2->Integral(0, 61);
		    
		  } else if (tr==T3) {
		    
		    ThisHistoName = HistogramName("muonPt_QCD", bpt, nb, T3, is, -1, -1);
		    TH1D *ThisQCDMuPt3 = (TH1D*) QCDFile->Get(ThisHistoName);
		    
		    LuminosityWeight = nT3/ThisQCDMuPt3->Integral(0, 61);
		    
		  } else {
		    
		    cout << "PtRelAnalyzer::BuildTemplates Warning:" 
			   << " error weighting trigger " << TriggerName[tr]
			 << " in jet pt bin " << PtRelPtBin[bpt] << endl; 
		    
		  }
		  
		} 
		
		MuonPt[0][is][fpt][nb]->Add(ThisDataMuPt);
		MuonDR[0][is][fpt][nb]->Add(ThisDataMuDR);
		JetPt[0][is][fpt][nb]->Add(ThisDataJetPt);
		JetEta[0][is][fpt][nb]->Add(ThisDataJetEta);
		PVEvent[0][is][fpt][nb]->Add(ThisDataPVMult);

		if (LuminosityWeight<0. || LuminosityWeight>1000000000.) {

		  cout << "PtRelAnalyzer::BuildTemplates: LuminosityWeight " << LuminosityWeight
		       << " for " << TriggerName[tr] << " in " << PtRelPtBin[bpt] << endl;

		  LuminosityWeight = 0.;

		}

		ThisQCDMuPt->Scale(LuminosityWeight);
		MuonPt[1][is][fpt][nb]->Add(ThisQCDMuPt);

		ThisQCDMuDR->Scale(LuminosityWeight);
		MuonDR[1][is][fpt][nb]->Add(ThisQCDMuDR);

		ThisQCDJetPt->Scale(LuminosityWeight);
		JetPt[1][is][fpt][nb]->Add(ThisQCDJetPt);

		ThisQCDJetEta->Scale(LuminosityWeight);
		JetEta[1][is][fpt][nb]->Add(ThisQCDJetEta);

		ThisQCDPVMult->Scale(LuminosityWeight);
		PVEvent[1][is][fpt][nb]->Add(ThisQCDPVMult);
 
		if (AddLightTemplates) {
		 
		  int nScaleBins = ThisQCDJetPt->GetNbinsX();

		  ThisHistoName = HistogramName("jetPt_Jet", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisJetJetPt = (TH1D*) JetFile->Get(ThisHistoName);

                  float JetLuminosityWeight = ThisQCDJetPt->Integral(0, nScaleBins)/ThisJetJetPt->Integral(0, nScaleBins);

                  ThisJetJetPt->Scale(JetLuminosityWeight);
		  JetPt[2][is][fpt][nb]->Add(ThisJetJetPt);
		  
		  ThisHistoName = HistogramName("jetEta_Jet", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisJetJetEta = (TH1D*) JetFile->Get(ThisHistoName);

                  ThisJetJetEta->Scale(JetLuminosityWeight);
		  JetEta[2][is][fpt][nb]->Add(ThisJetJetEta);
		  
		  ThisHistoName = HistogramName("jetPt_Inc", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisIncJetPt = (TH1D*) IncFile->Get(ThisHistoName);
		  
                  float IncLuminosityWeight = ThisQCDJetPt->Integral(0, nScaleBins)/ThisIncJetPt->Integral(0, nScaleBins);

                  ThisIncJetPt->Scale(IncLuminosityWeight);
		  JetPt[3][is][fpt][nb]->Add(ThisIncJetPt);

		  ThisHistoName = HistogramName("jetEta_Inc", bpt, nb, tr, is, -1, -1);
		  TH1D *ThisIncJetEta = (TH1D*) IncFile->Get(ThisHistoName);

                  ThisIncJetEta->Scale(IncLuminosityWeight);
		  JetEta[3][is][fpt][nb]->Add(ThisIncJetEta);
		  
		}

		for (int tg = 0; tg<nTaggers; tg++) {
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_Data", bpt, nb, tr, is, tg,  0);
		  TH1D *ThisDataPtRelTag = (TH1D*) DataFile->Get(ThisHistoName);
		  PtRel[tg][is][fpt][nb][0]->Add(ThisDataPtRelTag);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_Data", bpt, nb, tr, is, tg, 10);
		  TH1D *ThisDataPtRelUntag = (TH1D*) DataFile->Get(ThisHistoName);
		  PtRel[tg][is][fpt][nb][4]->Add(ThisDataPtRelUntag);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg,  1);
		  TH1D *ThisQCDPtRelTagB = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagB->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][0]->Add(ThisQCDPtRelTagB);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg, 11);
		  TH1D *ThisQCDPtRelUntagB = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagB->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][4]->Add(ThisQCDPtRelUntagB);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg,  2);
		  TH1D *ThisQCDPtRelTagC = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagC->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][1]->Add(ThisQCDPtRelTagC);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg, 12);
		  TH1D *ThisQCDPtRelUntagC = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagC->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][5]->Add(ThisQCDPtRelUntagC);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg,  3);
		  TH1D *ThisQCDPtRelTagLG = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelTagLG->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][2]->Add(ThisQCDPtRelTagLG);
		  
		  ThisHistoName = HistogramName(TemplateVariable + "_QCD", bpt, nb, tr, is, tg, 13);
		  TH1D *ThisQCDPtRelUntagLG = (TH1D*) QCDFile->Get(ThisHistoName);
		  ThisQCDPtRelUntagLG->Scale(LuminosityWeight);
		  PtRel[nTaggers+tg][is][fpt][nb][6]->Add(ThisQCDPtRelUntagLG);

		  if (AddLightTemplates) {

		    int nScaleBins = ThisQCDPtRelTagLG->GetNbinsX();
		    float TagLuminosityWeight   = ThisQCDPtRelTagLG->Integral(0, nScaleBins);
		    float UntagLuminosityWeight = ThisQCDPtRelUntagLG->Integral(0, nScaleBins);

		    ThisHistoName = HistogramName(TemplateVariable + "_Jet", bpt, nb, tr, is, tg,  4);
		    TH1D *ThisJetPtRelTagTrk = (TH1D*) JetFile->Get(ThisHistoName);
                    float JetTagLuminosityWeight = TagLuminosityWeight/ThisJetPtRelTagTrk->Integral(0, nScaleBins);
		    ThisJetPtRelTagTrk->Scale(JetTagLuminosityWeight);
		    PtRel[tg][is][fpt][nb][3]->Add(ThisJetPtRelTagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_Jet", bpt, nb, tr, is, tg, 14);
		    TH1D *ThisJetPtRelUntagTrk = (TH1D*) JetFile->Get(ThisHistoName);
                    float JetUntagLuminosityWeight = UntagLuminosityWeight/ThisJetPtRelUntagTrk->Integral(0, nScaleBins);
		    ThisJetPtRelUntagTrk->Scale(JetUntagLuminosityWeight);
		    PtRel[tg][is][fpt][nb][7]->Add(ThisJetPtRelUntagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_Inc", bpt, nb, tr, is, tg,  4);
		    TH1D *ThisIncPtRelTagTrk = (TH1D*) IncFile->Get(ThisHistoName);
                    float IncTagLuminosityWeight = TagLuminosityWeight/ThisIncPtRelTagTrk->Integral(0, nScaleBins);
		    ThisIncPtRelTagTrk->Scale(IncTagLuminosityWeight);
		    PtRel[nTaggers+tg][is][fpt][nb][3]->Add(ThisIncPtRelTagTrk);
		  
		    ThisHistoName = HistogramName(TemplateVariable + "_Inc", bpt, nb, tr, is, tg, 14);
		    TH1D *ThisIncPtRelUntagTrk = (TH1D*) IncFile->Get(ThisHistoName);
                    float IncUntagLuminosityWeight = UntagLuminosityWeight/ThisIncPtRelUntagTrk->Integral(0, nScaleBins);
		    ThisIncPtRelUntagTrk->Scale(IncUntagLuminosityWeight);
		    PtRel[nTaggers+tg][is][fpt][nb][7]->Add(ThisIncPtRelUntagTrk);

		  }

		}

	      }
		      
	  }
	  
	}

      }

    }
	    
    for (int nb = 1; nb<nPtRelEtaBins; nb++) { 
      
      for (int is = 0; is<nSystematics; is++) {
	
	for (int dt = 0; dt<2; dt++) {
	  
	  MuonPt[dt][is][fpt][0]->Add(MuonPt[dt][is][fpt][nb]);
	  MuonDR[dt][is][fpt][0]->Add(MuonDR[dt][is][fpt][nb]);
	  JetPt[dt][is][fpt][0]->Add(JetPt[dt][is][fpt][nb]);
	  JetEta[dt][is][fpt][0]->Add(JetEta[dt][is][fpt][nb]);
	  PVEvent[dt][is][fpt][0]->Add(PVEvent[dt][is][fpt][nb]);
	  
	}
	
	if (AddLightTemplates) {
	  
	  JetPt[2][is][fpt][0]->Add(JetPt[2][is][fpt][nb]);
	  JetEta[2][is][fpt][0]->Add(JetEta[2][is][fpt][nb]);
	  
	  JetPt[3][is][fpt][0]->Add(JetPt[3][is][fpt][nb]);
	  JetEta[3][is][fpt][0]->Add(JetEta[3][is][fpt][nb]);
	  
	}
	
	for (int tg = 0; tg<nTaggers; tg++) {
	  
	  PtRel[tg][is][fpt][0][0]->Add(PtRel[tg][is][fpt][nb][0]);
	  PtRel[tg][is][fpt][0][4]->Add(PtRel[tg][is][fpt][nb][4]);
	  
	  PtRel[nTaggers+tg][is][fpt][0][0]->Add(PtRel[nTaggers+tg][is][fpt][nb][0]);
	  PtRel[nTaggers+tg][is][fpt][0][4]->Add(PtRel[nTaggers+tg][is][fpt][nb][4]);
	  
	  PtRel[nTaggers+tg][is][fpt][0][1]->Add(PtRel[nTaggers+tg][is][fpt][nb][1]);
	  PtRel[nTaggers+tg][is][fpt][0][5]->Add(PtRel[nTaggers+tg][is][fpt][nb][5]);
	  
	  PtRel[nTaggers+tg][is][fpt][0][2]->Add(PtRel[nTaggers+tg][is][fpt][nb][2]);
	  PtRel[nTaggers+tg][is][fpt][0][6]->Add(PtRel[nTaggers+tg][is][fpt][nb][6]);
  
	  if (AddLightTemplates) {
	    
	    PtRel[tg][is][fpt][0][3]->Add(PtRel[tg][is][fpt][nb][3]);
	    PtRel[tg][is][fpt][0][7]->Add(PtRel[tg][is][fpt][nb][7]);
	    
	    PtRel[nTaggers+tg][is][fpt][0][3]->Add(PtRel[nTaggers+tg][is][fpt][nb][3]);
	    PtRel[nTaggers+tg][is][fpt][0][7]->Add(PtRel[nTaggers+tg][is][fpt][nb][7]);
	    
	  }
	  
	}
	
      }
      
    }
      
  }

  TString LightTemplates = "";
  if (AddLightTemplates) LightTemplates = "All";	
  TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + LightTemplates + PUWeighting + KinWeighting + "EtaBins" + Selection + Production + ".root";
  
  SaveTemplates(TemplateFileName, AddLightTemplates);
  
}

void PtRelAnalyzer::SaveTemplates(TString TemplateFileName, bool AddLightTemplates) {

  cout << "Saving templates" << endl;

  TFile *OutFile = new TFile(TemplateFileName, "recreate");
  
  for (int fpt = 0; fpt<nFitPtBins; fpt++)
    for (int nb = 0; nb<nPtRelEtaBins; nb++) 
      for (int is = 0; is<nSystematics; is++) 
	for (int dt = 0; dt<2; dt++) {
	  
	  JetPt[dt][is][fpt][nb]->Write();
	  JetEta[dt][is][fpt][nb]->Write();

	  if (AddLightTemplates) {

	    JetPt[dt+2][is][fpt][nb]->Write();
	    JetEta[dt+2][is][fpt][nb]->Write();
	    
	  }

	  MuonPt[dt][is][fpt][nb]->Write();
	  MuonDR[dt][is][fpt][nb]->Write();
	  PVEvent[dt][is][fpt][nb]->Write();
	    
	  for (int tg = 0; tg<nTaggers; tg++)
	    for (int tp = 0; tp<2; tp++)
	      for (int fl = 0; fl<4; fl++)
		if ((dt==1 && fl<3) || fl==0 || (fl==3 && AddLightTemplates))
		  PtRel[dt*nTaggers+tg][is][fpt][nb][4*tp+fl]->Write();
	  
	}
  
  OutFile->Close();

  cout << "Exiting" << endl;
  
}

void PtRelAnalyzer::ComputeKinematicWeights(TString SystematicFlag, TString DataType, TString LightTemplates) {

  cout << "Producing kinematic weights" << endl;
  
  if (KinWeighting.Contains("_KinPtEtaBins")) {
	
    TString TemplateFileName = "./Templates/" + TemplateVariable + "_Templates" + LightTemplates + PUWeighting + Selection + "_BaseProduction.root";
    TFile *TemplateFile = TFile::Open(TemplateFileName); 
    cout << TemplateFileName << endl;
    bool GoodSystematicFlag = false;

    for (int is = 0; is<nSystematics; is++) 
      if (SystematicFlag==SystematicName[is]) {
	
	for (int bpt = 0; bpt<nPtRelPtBins; bpt++) 
	  for (int nb = 0; nb<nPtRelEtaBins; nb++) {

            //if (TemplateVariable=="PtRel" && bpt>=9) continue;
            //if (TemplateVariable=="IP3D"  && bpt<=8) continue;

	    cout << "  Weights for " << PtRelPtBin[bpt] << " " << PtRelEtaBin[nb] << endl;
	    
	    TString ThisRefPtName = HistogramName("jetPt_" + DataType, PtRelPtBin[bpt], nb, -1, is, -1, -1);
	    TH1D *ThisRefPt = (TH1D*) TemplateFile->Get(ThisRefPtName);
	    
	    TString ThisRefEtaName = HistogramName("jetEta_" + DataType, PtRelPtBin[bpt], nb, -1, is, -1, -1);
	    TH1D *ThisRefEta = (TH1D*) TemplateFile->Get(ThisRefEtaName);
	    
	    TString ThisDataPtName = HistogramName("jetPt_Data", PtRelPtBin[bpt], nb, -1, is, -1, -1);
	    TH1D *ThisDataPt = (TH1D*) TemplateFile->Get(ThisDataPtName);
	    
	    TString ThisDataEtaName = HistogramName("jetEta_Data", bpt, nb, -1, is, -1, -1);
	    TH1D *ThisDataEta = (TH1D*) TemplateFile->Get(ThisDataEtaName);
	    
	    float DataEvents = ThisDataPt->Integral();
	    float RefEvents = ThisRefPt->Integral();
	    ThisRefPt->Scale(DataEvents/RefEvents);
	    
	    int PtRebin = 2;
	    if (PtRelPtEdge[bpt]>=120.) PtRebin = 4;
	    if (PtRelPtEdge[bpt]>=160.) PtRebin = 10;
	    if (PtRelPtEdge[bpt]>=320.) PtRebin = 20;

	    ThisDataPt->Rebin(PtRebin);
	    ThisRefPt->Rebin(PtRebin);
	    ThisDataPt->Divide(ThisRefPt);
	    cout << "Pt " << DataEvents << " " << RefEvents << endl;
	    WriteKinematicWeights(ThisDataPt, "_Pt", bpt, nb, DataType);
	    
	    DataEvents = ThisDataEta->Integral();
	    RefEvents = ThisRefEta->Integral();
	    ThisRefEta->Scale(DataEvents/RefEvents);
	    
	    int EtaRebin = 6;
	    ThisDataEta->Rebin(EtaRebin);
	    ThisRefEta->Rebin(EtaRebin);
	    ThisDataEta->Divide(ThisRefEta);
	    cout << "Eta " << DataEvents << " " << RefEvents << endl;
	    WriteKinematicWeights(ThisDataEta, "_Eta", bpt, nb, DataType);
	   	    
	    GoodSystematicFlag = true;
	    
	  }
	
      }
	
    if (!GoodSystematicFlag) 
      cout << "PtRelAnalyzer::ComputeKinematicWeights: SystematicFlag " << SystematicFlag << " not supported" << endl;
    
  }
  
}

void PtRelAnalyzer::WriteKinematicWeights(TH1D *&HistoDivide, TString Variable, int PtBin, int EtaBin, TString DataType) {

  cout << "Writing kinematic weights" << endl;

  ofstream KinWeightsFile; 
  KinWeightsFile.open("./Weights/KinematicWeights/" + TemplateVariable + KinWeighting + Variable + "_" + PtRelPtBin[PtBin] + "_" + PtRelEtaBin[EtaBin] + "_" + DataType + PUWeighting + Selection + ".txt");
  
  int nSteps; float FirstStep, StepSize;
  if (Variable=="_Pt") { 

    nSteps = PtRelPtEdge[PtBin+1] - PtRelPtEdge[PtBin];
    FirstStep = PtRelPtEdge[PtBin];
    StepSize = 1.;

  } else { 

    nSteps = 48;
    FirstStep = -PtRelEtaEdge[nPtRelEtaBins-1];
    StepSize = 0.1;

  }

  for (int st = 0; st<nSteps; st++) {

    float ThisStep = FirstStep + st*StepSize;
    int ThisStepBin = HistoDivide->FindBin(ThisStep);
    
    float ThisWeight = HistoDivide->GetBinContent(ThisStepBin);
    if (ThisWeight<0.) ThisWeight = 1.;

    KinWeightsFile << ThisStep << " " << ThisWeight << endl; 

  }

  KinWeightsFile.close();

}

bool PtRelAnalyzer::PassTriggerEmulation(int TriggerIdx, int MuonJetIdx) {

  if (TriggerName[TriggerIdx].Contains("_Jet")) return true; 
  else if (TriggerName[TriggerIdx].Contains("_DiJet")) {
    
    int tJet = GetAwayJet("NONE", MuonJetIdx, 0., false);

    if (tJet>=0) {    
      
      float PtThreshold = 1000000.;
      if (TriggerName[TriggerIdx].Contains("_DiJet20"))       PtThreshold =  30.;
      else if (TriggerName[TriggerIdx].Contains("_DiJet40"))  PtThreshold =  50.;
      else if (TriggerName[TriggerIdx].Contains("_DiJet70"))  PtThreshold =  80.;
      else if (TriggerName[TriggerIdx].Contains("_DiJet110")) PtThreshold = 120.;
      else 
	cout << "PtRelAnalyzer::PassTriggerEmulation: trigger " << TriggerName[TriggerIdx] << " not supported" << endl;
     
      if (Jet_pt[tJet]>=PtThreshold) return true;
      
    }

  } else
    cout << "PtRelAnalyzer::PassTriggerEmulation: trigger " << TriggerName[TriggerIdx] << " not supported" << endl;

  return false;

}

FitResult PtRelAnalyzer::FitResultFromTable(TString TableName) {

  struct FitResult ThisFitResult;
  
  ifstream Table; Table.open(TableName);
  
  TString EtaBin, PtBin;
  Table >> EtaBin >> PtBin;
  
  TString Eff, MC, Eq, PM; 
  double ThisEff, ThisEffError;
  Table >> Eff >> MC >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisFitResult.EffMC = ThisEff; ThisFitResult.EffMCError = ThisEffError;

  TString Chi2NString; 
  double Chi2N_Untag, Chi2N_Tag;
  Table >> Chi2NString >> Eq >> Chi2N_Untag >> Chi2N_Tag;
    
  ThisFitResult.Chi2NTag = Chi2N_Tag; ThisFitResult.Chi2NUntag = Chi2N_Untag;

  TString bottom, frac;
  Table >> bottom >> frac >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisFitResult.FracUntag = ThisEff; ThisFitResult.FracUntagError = ThisEffError;
  
  TString tag;
  Table >> bottom >> frac >> tag >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisFitResult.FracTag = ThisEff; ThisFitResult.FracTagError = ThisEffError;
  
  TString Data;
  Table >> Eff >> Data >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisFitResult.EffData = ThisEff; ThisFitResult.EffDataError = ThisEffError;
  
  TString OpenPar, Raw1, Raw2, Raw3, Eq2, EffString;
  Table >> OpenPar >> Eff >> Eq >> Raw1 >> Raw2 >> Raw3 >> Eq2 >> EffString;
  
  EffString.ReplaceAll(")", "");
  
  ThisFitResult.EffDataRaw = EffString.Atof();
  
  TString Scale, Factor;
  Table >> Scale >> Factor >> Eq >> ThisEff >> PM >> ThisEffError;
  
  ThisFitResult.SF = ThisEff; ThisFitResult.SFError = ThisEffError;

  return ThisFitResult;
  
}

TString PtRelAnalyzer::FitResultTableName(TString Tagger, TString EtaBin, TString PtBin, TString Systematic, TString PUWeightingFlag, TString KinWeightingFlag, TString SelectionFlag, TString ProductionFlag) {

  return FitResultTableName(Tagger, EtaBin, PtBin, Systematic, PUWeightingFlag, KinWeightingFlag + SelectionFlag, ProductionFlag);

}

TString PtRelAnalyzer::FitResultTableName(TString Tagger, TString EtaBin, TString PtBin, TString Systematic, TString PUWeightingFlag, TString Configuration, TString ProductionFlag) {

    TString TableName = "./Tables/" + TemplateVariable + "Fit_" + Tagger + "_" + EtaBin + "_" + PtBin + Systematic + PUWeightingFlag + Configuration + ProductionFlag + ".txt";
    if (TableName.Contains("JEC")) TableName.ReplaceAll(".txt", "_JEC.txt");
    else TableName.ReplaceAll(".txt", "_DP.txt");
    
    return TableName;

}

TCanvas *PtRelAnalyzer::BTagPerformanceCanvas(TString Type) {

  gStyle->SetTitleBorderSize( 0);
  gStyle->SetTitleFillColor (10);
  gStyle->SetStatFont       (42);
  gStyle->SetTitleFont      (42);

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat("");
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);

  float CanvasHeight = 400.;
  if (Type=="Performance") CanvasHeight = 800.;

  TCanvas* ThisCanvas = new TCanvas("BTagPerformance", "BTagPerformance", 1200., CanvasHeight);
  
  if (Type=="Performance") {

    ThisCanvas->Divide(1,2,0,0);

    TPad *EfficiencyPad = (TPad*) ThisCanvas->GetPad(1);
    EfficiencyPad->SetPad(0.01, 0.51, 0.95, 0.95);

    TPad *ScaleFactorsPad = (TPad*) ThisCanvas->GetPad(2);
    ScaleFactorsPad->SetPad(0.01, 0.07, 0.95, 0.51);

  } else {

    ThisCanvas->Divide(1,1,0,0);

    TPad *ScaleFactorsPad = (TPad*) ThisCanvas->GetPad(1);
    ScaleFactorsPad->SetPad(0.01, 0.07, 0.95, 0.95);

  }
   
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);
  gPad->SetLeftMargin(0.19);
  //gPad->SetGrid();
    
  gPad->SetTickx();
  gPad->SetTicky();

  return ThisCanvas;
  
}

void PtRelAnalyzer::FillBTagPerformanceHistograms(TString Tagger, TString EtaBin, TString Configuration, TString Systematic, int ColorIdx, TString DrawedSystematics) {

  MCEfficiency->Reset();
  DataEfficiency->Reset();
  DataMCSF->Reset();

  float MinEfficiency = 1.;
  float MaxEfficiency = 0.;

  for (int ppt = 2; ppt<MCEfficiency->GetNbinsX(); ppt++) {

    TString TableName = FitResultTableName(Tagger, EtaBin, FitPtBin[ppt-1], Systematic, PUWeighting, Configuration, Production);

    FitResult ThisFitResult;

    ThisFitResult.EffMC = 0.; ThisFitResult.EffData = 0.; ThisFitResult.SF = 0.;

    //cout << TableName <<" " << ThisFitResult.SF<<endl;
    ThisFitResult = FitResultFromTable(TableName);
    //cout << TableName <<" " << ThisFitResult.SF<<endl;
    MCEfficiency->SetBinContent(ppt+1, ThisFitResult.EffMC);
    MCEfficiency->SetBinError(ppt+1, ThisFitResult.EffMCError);

    DataEfficiency->SetBinContent(ppt+1, ThisFitResult.EffData);
    DataEfficiency->SetBinError(ppt+1, ThisFitResult.EffDataError);

    DataMCSF->SetBinContent(ppt+1, ThisFitResult.SF);
    DataMCSF->SetBinError(ppt+1, ThisFitResult.SFError);

    if (ThisFitResult.EffMC<MinEfficiency) MinEfficiency = ThisFitResult.EffMC;
    if (ThisFitResult.EffData<MinEfficiency) MinEfficiency = ThisFitResult.EffData;

    if (ThisFitResult.EffMC>MaxEfficiency) MaxEfficiency = ThisFitResult.EffMC;
    if (ThisFitResult.EffData>MaxEfficiency) MaxEfficiency = ThisFitResult.EffData;

  }

  MCEfficiency->GetYaxis()->SetRangeUser(0.05*(int(MinEfficiency/0.05)-2), 0.05*(int(MaxEfficiency/0.05)+3));

  MCEfficiency->SetMarkerSize(1.);
  MCEfficiency->SetMarkerStyle(25);
  MCEfficiency->SetMarkerColor(ColorIdx);
  
  MCEfficiency->GetXaxis()->SetLabelFont(42);
  MCEfficiency->GetYaxis()->SetLabelFont(42);
  MCEfficiency->GetXaxis()->SetTitleFont(42);
  MCEfficiency->GetYaxis()->SetTitleFont(42);
  MCEfficiency->SetTitle("");
  MCEfficiency->SetXTitle("Jet p_{T} [GeV]");
  MCEfficiency->SetYTitle("b-tag Efficiency #epsilon_{b}");
  MCEfficiency->GetYaxis()->SetTitleSize(0.09);
  MCEfficiency->GetYaxis()->SetTitleOffset(0.5);
  MCEfficiency->GetYaxis()->SetLabelSize(0.08);
  
  DataEfficiency->SetMarkerSize(1.);
  DataEfficiency->SetMarkerStyle(20);
  DataEfficiency->SetMarkerColor(ColorIdx);
  
  DataMCSF->SetMarkerSize(1.);
  DataMCSF->SetMarkerStyle(20);
  DataMCSF->SetMarkerColor(ColorIdx);
  
  DataMCSF->GetXaxis()->SetLabelFont(42);
  DataMCSF->GetYaxis()->SetLabelFont(42);
  DataMCSF->GetXaxis()->SetTitleFont(42);
  DataMCSF->GetYaxis()->SetTitleFont(42);
  DataMCSF->SetTitle("");
  DataMCSF->SetXTitle("Jet p_{T} [GeV]");
  DataMCSF->SetYTitle(Tagger + "Data/Sim. SF_{b}");
  DataMCSF->GetXaxis()->SetTitleOffset(1.3);
  DataMCSF->GetXaxis()->SetTitleSize(0.09);
  DataMCSF->GetXaxis()->SetLabelSize(0.08);
  DataMCSF->GetYaxis()->SetTitle("Data/Sim. SF_{b}");
  DataMCSF->GetYaxis()->SetTitleSize(0.09);
  DataMCSF->GetYaxis()->SetTitleOffset(0.5);
  DataMCSF->GetYaxis()->CenterTitle(true);
  DataMCSF->GetYaxis()->SetLabelSize(0.09); 
  
  DataMCSF->SetMinimum(0.6);
  DataMCSF->SetMaximum(1.4);

  if (DrawedSystematics!="NONE") {

    DataMCSFSystematics->Add(DataMCSF);
    
    ComputeScaleFactorSystematics(Tagger, EtaBin, Configuration, DrawedSystematics);
 
    for (int nb = 0; nb<nPtRelEtaBins; nb++)
      if (EtaBin==PtRelEtaBin[nb]) 
	for (int ppt = 1; ppt<MCEfficiency->GetNbinsX(); ppt++)
	  DataMCSFSystematics->SetBinError(ppt+1, TotalScaleFactorSystematic[ppt-1][nb]);

    DataMCSFSystematics->SetMarkerStyle(20);
    DataMCSFSystematics->SetMarkerColor(2);
    DataMCSFSystematics->SetFillStyle(3005);
    DataMCSFSystematics->SetFillColor(ColorIdx+1);
    
  }
  
}

void PtRelAnalyzer::PlotBTagPerformance(TString PlotFlag, TString Tagger, TString EtaList[], TString ConfigurationList[], TString SystematicList[], float MaxJetPt, TString DrawedSystematics, TString Type) {

  TCanvas *ThisCanvas = BTagPerformanceCanvas(Type);
  
  int nPlotPtBins = 0;
  float PlotPtEdge[1000];
  
  PlotPtEdge[0] = 0.;
  
  for (int fpt = 0; fpt<nFitPtBins+1; fpt++) 
    if (FitPtEdge[fpt]<=MaxJetPt) {
      
      nPlotPtBins++;
      PlotPtEdge[nPlotPtBins] = FitPtEdge[fpt];
      
    }
  
  MCEfficiency = new TH1D("MCEfficiency", "", nPlotPtBins, PlotPtEdge);
  DataEfficiency = new TH1D("DataEfficiency", "", nPlotPtBins, PlotPtEdge);
  DataMCSF = new TH1D("DataMCSF", "", nPlotPtBins, PlotPtEdge);
  DataMCSFSystematics = new TH1D("DataMCSFSystematics", "", nPlotPtBins, PlotPtEdge);
  
  TString Option;
  
  if (ConfigurationList[0]=="") ConfigurationList[0] = KinWeighting + Selection;

  int ColorIndex = 1;

  for (int nb = 0; nb<100 && EtaList[nb]!="-"; nb++)   
    for (int cfg = 0; cfg<100 && ConfigurationList[cfg]!="-"; cfg++)
      for (int is = 0; is<100 && SystematicList[is]!="-"; is++) {
	
	FillBTagPerformanceHistograms(Tagger, EtaList[nb], ConfigurationList[cfg], SystematicList[is], ColorIndex, DrawedSystematics);
	
	ThisCanvas->cd(1);

	gPad->SetTickx();
	gPad->SetTicky();
	
	if (Type!="ScaleFactors") {
	  
	  MCEfficiency->DrawCopy("pe1" + Option);
	  DataEfficiency->DrawCopy("pe1same");
	  
	}
	
	if (Type=="Performance") {

	  ThisCanvas->cd(2);
	
	  gPad->SetTickx();
	  gPad->SetTicky();

	}

	if (Type!="Efficiency") {

	  DataMCSF->DrawCopy("pe1" + Option);
	  
	  if (DrawedSystematics!="NONE")  {

	    DataMCSFSystematics->DrawCopy("pe2same");
	    DataMCSF->DrawCopy("pe1same");

	  }

	}

	Option = "same";
	
	ColorIndex++;

      }

  TString ThisPlotName = "BTag" + Type + "_DepOn" + PlotFlag + "_" + Tagger + "_" + EtaList[0] + PUWeighting + ConfigurationList[0] + Production + ".png";
  
  if (EtaList[1]!="-") ThisPlotName.ReplaceAll("_" + EtaList[0], "");
  if (ConfigurationList[1]!="-") ThisPlotName.ReplaceAll(ConfigurationList[0], "");

  ThisCanvas->Print("./Plots/" + ThisPlotName);

  cout << "http://scodella.web.cern.ch/scodella/Work/CMS/BTagging/ptRel/Run2015B/" << ThisPlotName << endl;
  
}

void PtRelAnalyzer::ComputeScaleFactorSystematics(TString Tagger, int PtBin, int EtaBin, TString Configuration, TString CentralFlag) {

  TString TableNameCentral = FitResultTableName(Tagger, PtRelEtaBin[EtaBin], FitPtBin[PtBin], CentralFlag, PUWeighting, Configuration, Production);
  FitResult CentralFitResult = FitResultFromTable(TableNameCentral);
  
  float CentralScaleFactor = CentralFitResult.SF;
  
  ScaleFactorValue[PtBin][EtaBin] = CentralScaleFactor;

  ScaleFactorSystematic[PtBin][EtaBin][0] = CentralFitResult.SFError;

  TotalScaleFactorSystematic[PtBin][EtaBin] = CentralFitResult.SFError*CentralFitResult.SFError;
  
  for (int sfs = 1; sfs<nScaleFactorSystematics; sfs++) {

    ScaleFactorSystematic[PtBin][EtaBin][sfs] = 0.;
    
    TString VariationSign = "down";

    for (int is = 2*sfs-1; is<2*sfs+1; is++) {

      if (ScaleFactorSystematicName[sfs]=="_MuPt" || ScaleFactorSystematicName[sfs]=="_AwayTagger")
	VariationSign = "up";
      
      TString SystematicTableName = TableNameCentral; 
      SystematicTableName.ReplaceAll("_Central", SystematicName[is]);
      FitResult SystematicFitResult = FitResultFromTable(SystematicTableName);
      
      float SystematicScaleFactor = SystematicFitResult.SF - CentralScaleFactor;
      float SystematicSFError = SystematicFitResult.SFError;
      
      if (SystematicSFError<0.13)
	if (fabs(SystematicScaleFactor)>fabs(ScaleFactorSystematic[PtBin][EtaBin][sfs])) {

	  if (VariationSign=="down") ScaleFactorSystematic[PtBin][EtaBin][sfs] = -1.*SystematicScaleFactor;
	  else ScaleFactorSystematic[PtBin][EtaBin][sfs] = SystematicScaleFactor;

	}

      VariationSign = "up";

    }
    
    if (VerboseSystematics)
      cout << Tagger << " " << FitPtBin[PtBin] << " " << ScaleFactorSystematicName[sfs] << " " << ScaleFactorSystematic[PtBin][EtaBin][sfs] <<endl;
 
    TotalScaleFactorSystematic[PtBin][EtaBin] += ScaleFactorSystematic[PtBin][EtaBin][sfs]*ScaleFactorSystematic[PtBin][EtaBin][sfs];
    
  }
  
  TotalScaleFactorSystematic[PtBin][EtaBin] = sqrt(TotalScaleFactorSystematic[PtBin][EtaBin]);
  
}

void PtRelAnalyzer::ComputeScaleFactorSystematics(TString Tagger, int PtBin, int EtaBin, TString CentralFlag) {

  ComputeScaleFactorSystematics(Tagger, PtBin, EtaBin, KinWeighting + Selection, CentralFlag);

}

void PtRelAnalyzer::ComputeScaleFactorSystematics(TString Taggers, TString EtaBin, TString CentralFlag) {
  
  ComputeScaleFactorSystematics(Taggers, EtaBin, KinWeighting + Selection, CentralFlag);
  
}

void PtRelAnalyzer::ComputeScaleFactorSystematics(TString Taggers, TString EtaBin, TString Configuration, TString CentralFlag) {
  
  for (int tg = 0; tg<nTaggers; tg++)
    if (Taggers=="All" || TaggerName[tg].Contains(Taggers) || Taggers.Contains(TaggerName[tg]))
      for (int fpt = 0; fpt<nFitPtBins; fpt++)
	for (int nb = 0; nb<nPtRelEtaBins; nb++)
	  if (EtaBin==PtRelEtaBin[nb])
	    ComputeScaleFactorSystematics(TaggerName[tg], fpt, nb, Configuration, CentralFlag);
  
}
