TString CampaignLuminosity = "59.97 fb^{-1} (13 TeV)";
TString CampaignName       = "UL18";

TString JECVersionData = "Summer19UL18_RunXXX_V5";
TString JECVersionMC   = "Summer19UL18_V5";

// Eta bins
const int nPtRelEtaBins = 3;
TString PtRelEtaBin[nPtRelEtaBins] = {"anyEta", "Eta12", "Eta25"};
double PtRelEtaEdge[nPtRelEtaBins] = {0., 1.2, 2.5};

// Pt bins
const int nPtRelPtBins = 20;
TString PtRelPtBin[nPtRelPtBins] = {"Pt2030", "Pt3040", "Pt4050", "Pt5060", "Pt6070", "Pt7080", "Pt80100", "Pt100120", "Pt120140", "Pt140160", "Pt160200", "Pt200260", "Pt260300", "Pt300320", "Pt320400", "Pt400500", "Pt500600", "Pt600800", "Pt8001000", "Pt1000"};
double PtRelPtEdge[nPtRelPtBins+1] = {20., 30., 40., 50., 60., 70., 80., 100., 120., 140., 160., 200., 260., 300., 320., 400., 500., 600., 800., 1000., 1400.};
const int MaxPtRelPtEdge =  1400;

double MuonPtCut[nPtRelPtBins] = {5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.};

const int nMaxPU = 99;
const int nMaxPrescales = 7800;

