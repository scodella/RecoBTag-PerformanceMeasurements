#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

#include "CampaignParameters/Run2015B/BaseParameters.h"
#include "CampaignParameters/Run2015B/TriggerInfo.h"
#include "CampaignParameters/Run2015B/Taggers.h"

// Choose the production 
//#include "CampaignParameters/Run2015B/BaseProduction.h"
//#include "CampaignParameters/Run2015B/RunIProduction.h"
#include "CampaignParameters/Run2015B/Run2015BProduction.h"
//#include "CampaignParameters/Run2015B/FlatProduction.h"
//#include "CampaignParameters/Run2015B/FlatProduction300.h"

#include "CampaignParameters/Run2015B/Systematics.h"
  
float TotalScaleFactorSystematic[nFitPtBins][nPtRelEtaBins];
float ScaleFactorSystematic[nFitPtBins][nPtRelEtaBins][nScaleFactorSystematics];
float ScaleFactorValue[nFitPtBins][nPtRelEtaBins];

TString TemplateVariable;

struct FitResult {double EffData; double EffDataRaw; double EffDataError; double EffMC; double EffMCError; double SF; double SFError; double FracTag; double FracTagError; double FracUntag; double FracUntagError; double Chi2NTag; double Chi2NUntag;};

class PtRelAnalyzer {
  
 public:
  
  PtRelAnalyzer(TString TemplateFlag, TString PileUpReweighting = "", TString KinematicWeighting = "", TString SelectionFlag = "");
  ~PtRelAnalyzer();

  void BookHistograms();
  void FillHistograms(TString DataType, int DataRange);
  void MergeHistograms(TString DataType, int FirstBin = 0, int LastBin = 100);

  void BookLightHistograms();
  void FillLightHistograms(TString DataType, int DataRange);
  void MergeLightHistograms(TString DataType, int FirstBin = 0, int LastBin = 100);

  void BuildTemplates(bool AddLightTemplates = false);
  void BuildTemplatesEtaBins(bool AddLightTemplates = false);

  void BookSystem8Histograms();
  void FillSystem8Histograms(TString DataType, int DataRange);
  void MergeSystem8Histograms(TString DataType, int FirstBin = 0, int LastBin = 100);

  void BuildSystem8Templates();
  
  void ComputeKinematicWeights(TString SystematicFlag = "_Central", TString DataType = "QCD", TString LightTemplates = "");


  void PlotBTagPerformance(TString PlotFlag, TString Tagger, TString EtaList[], TString ConfigurationList[], TString SystematicList[], float MaxJetPt = 1000000., TString DrawedSystematics = "NONE", TString Type = "Performance");

  void ComputeScaleFactorSystematics(TString Taggers = "All", TString EtaBin = "anyEta", TString CentralFlag = "_Central");
  
 private:

  TString PUWeighting, KinWeighting, Selection;

  int nBinsForTemp;
  float LowerEdgeForTemp, UpperEdgeForTemp;

  int GetAwayJet(TString ThisAwayTaggerName, int jMu, float AwayDRCut, bool IsUnique);
  void GetDRCuts(TString ThisSystematicName, float JetPt, float *TrackDRCut, float *TrackMinDRCut);

  double KinematicWeight[nTriggers][nSystematics][600][60][nPtRelEtaBins];
  float PileUpWeight[nTriggers][nMaxPU][nSystematics];
  float BTemplateCorrections[100][nPtRelPtBins][2]; 

  void ResetHistograms(TString DataType);
  void SaveHistograms(TString OutputFileName, TString DataType);

  void ResetLightHistograms(TString DataType);
  void SaveLightHistograms(TString OutputFileName);

  void ResetSystem8Histograms(TString DataType);
  void SaveSystem8Histograms(TString OutputFileName, TString DataType);

  void GetKinematicWeights(TString DataType);
  void GetPileUpWeights(TString DataType);
  void GetBTemplateCorrections();

  TString HistogramName(TString VariableName, TString PtBin, int EtaBin, int TriggerIdx, int SystematicIdx, int TaggerIdx, int FlavourIdx);
  TString HistogramName(TString VariableName, int PtBin, int EtaBin, int TriggerIdx, int SystematicIdx, int TaggerIdx, int FlavourIdx);

  void BookTemplates();
  void SaveTemplates(TString TemplateFileName, bool AddLightTemplates);

  void BookSystem8Templates();
  void SaveSystem8Templates(TString TemplateFileName);

  bool PassTriggerEmulation(int TriggerIdx, int MuonJetIdx);

  // Histograms
  TH1D *MuDRForWeighting[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *MuPtForWeighting[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *JetPtForWeighting[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *JetEtaForWeighting[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PVMultiplicity[nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  
  TH1D *PtRelTagForWeighting[2*nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PtRelUntagForWeighting[2*nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PtRelLightTagForWeighting[2*nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *PtRelLightUntagForWeighting[2*nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
 
  TH1D *PtRelLightTagForSystematic[nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PtRelLightUntagForSystematic[nTaggers][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  
  TH2D *Observed[nPtRelPtBins][nPtRelEtaBins][nSystematics];

  // Templates
  TH1D *JetPt[4][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *JetEta[4][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *MuonPt[2][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *MuonDR[2][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *PVEvent[2][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *PtRel[2*nTaggers][nSystematics][nPtRelPtBins][nPtRelEtaBins][8]; 

  TH1D *System8[2*nTaggers+2][nSystematics][nPtRelPtBins][nPtRelEtaBins][3]; 

  TH1D *System8ForDataWeighting[2*nTaggers+2][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];
  TH1D *System8ForBJetWeighting[2*nTaggers+2][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins]; 
  TH1D *System8ForLJetWeighting[2*nTaggers+2][nTriggers][nSystematics][nPtRelPtBins][nPtRelEtaBins];

  void WriteKinematicWeights(TH1D *&HistoDivide, TString Variable, int PtBin, int EtaBin, TString DataType);

  FitResult FitResultFromTable(TString TableName);
  TString FitResultTableName(TString Tagger, TString EtaBin, TString PtBin, TString Systematic, TString PUWeightingFlag, TString KinWeightingFlag, TString SelectionFlag, TString ProductionFlag);
  TString FitResultTableName(TString Tagger, TString EtaBin, TString PtBin, TString Systematic, TString PUWeightingFlag, TString Configuration, TString ProductionFlag);

  TCanvas *BTagPerformanceCanvas(TString Type);
  TH1D *MCEfficiency, *DataEfficiency, *DataMCSF, *DataMCSFSystematics;
  void FillBTagPerformanceHistograms(TString Tagger, TString EtaBin, TString Configuration, TString Systematic, int ColorIdx, TString DrawedSystematics);

  void ComputeScaleFactorSystematics(TString Tagger, int PtBin, int EtaBin, TString CentralFlag);
  void ComputeScaleFactorSystematics(TString Tagger, int PtBin, int EtaBin, TString Configuration, TString CentralFlag);
  void ComputeScaleFactorSystematics(TString Tagger, TString EtaBin, TString Configuration, TString CentralFlag);

};

