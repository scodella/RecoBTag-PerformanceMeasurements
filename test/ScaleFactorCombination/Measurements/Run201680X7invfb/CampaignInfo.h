// General data info
string CampaignNameString = "ICHEP16";
TString CampaignName = CampaignNameString;
TString CampaignLuminosity = "7.7 fb^{-1},"; 
TString CenterOfMassEnergy = " #sqrt{s}=13 TeV, 2016";

// Taggers
int const nTaggingAlgorithms = 2;
string TaggingAlgorithmName[nTaggingAlgorithms] = {"CSVv2", "cMVAv2"};

// Pt bins
int const nBinsCampaign = 7;
float MinPtCampaign, MaxPtCampaign;
float xPt[nBinsCampaign]  = {40., 60., 85., 120., 170., 250., 485.};
float exPt[nBinsCampaign] = {10., 10., 15.,  20.,  30.,  50., 185.};
/* int const nBinsCampaign = 3; */
/* float MinPtCampaign, MaxPtCampaign; */
/* float xPt[nBinsCampaign]  = {45., 90., 220.}; */
/* float exPt[nBinsCampaign] = {15., 30., 100.}; */

// Some general uncertainty options
bool StatisticalCorrelation = true;
bool InflateStatistic       = false;
bool AddSampleDependence    = false;
float SampleDependence      = 0.00; 
float cJetsInflationFactor[3] = {2., 2., 2.};

// Measurements
int const nMeasurements = 9;
string MeasurementName[nMeasurements] = {"PtRel",  "System8", "LT",     "TagCount", "Kin",   "TnP",   "IF",    "TnPH",  "TnPL"};
string MeasurementFlag[nMeasurements] = {"ptrel",  "system8", "lt",     "TagCount", "kin",   "TnP",   "if",    "TnPH",  "TnPL"};
string MeasurementPlot[nMeasurements] = {"_PT",    "_S8",     "_LT",    "_TC",      "_Kin",  "_TnP",  "_IF",   "_TnPH", "_TnPL"};
string TypeMeasurement[nMeasurements] = {"mujets", "mujets",  "mujets", "ttbar",    "ttbar", "ttbar", "ttbar", "ttbar", "ttbar"};
bool UseThisMeasurement[nTaggingAlgorithms][nMeasurements] = { { true,  true,  true,  true,  true,  true,   true,  true,  true},
							       {false, false, false,  true, false,  true,   true,  true,  true} };
bool VetoedMeasurement [nTaggingAlgorithms][nMeasurements] = { {false, false, false,  true, false,  true,   true,  true,  true}, 
							       {false, false, false, false, false,  true,   true, false, false} };
bool PlotTheMeasurement[nTaggingAlgorithms][nMeasurements] = { { true,  true,  true, false,  true, false,  false, false, false},
							       {false, false, false,  true, false,  true,  false,  true,  true} };

// Plot parameters
int GraphXBins[nMeasurements] = {7, 4, 7, 6, 5, 6, 3, 6, 6};
int GraphStyle[nMeasurements] = {3005, -1, -1, -1, -1, -1, -1, -1, -1};
int GraphColor[nMeasurements] = {kBlue, kGreen, kRed, kOrange, kBlack, 28, 46, 29, 30};
int GraphMarker[nMeasurements] = {22, 29, 20, 21, 33, 23, 24, 23, 23};
float GraphSize[nMeasurements] = {1.3, 1.3, 1.3, 1., 1., 1., 1., 1., 1.};
int GraphWidth[nMeasurements] = {6, 6, 4, 4, 4, 4, 4, 4, 4};
int MeasurementOrder[nMeasurements] = {1, 0, 3, 5, 2, 4, 6, 7, 8};


// Uncertainties
int const nUncertainties = 40;
string UncertaintyName[nUncertainties] = {"", "_statistic", "_pileup", "_mupt", "_gluonsplitting", "_bfragmentation", "_cfragmentation", "_jetaway", "_mudr", "_cb", "_jes", "_ttbarmodelling", "_l2c", "_ptrel", "_ipbias", "_ltothers", "_jer", "_tt", "_kt", "_tct", "_ntrk", "_njet", "_jeta", "_dmux", "_ksl", "_ntrkgen", "_btempcorr", "_ltempcorr", "_flavFrac", "_backg", "_scale1", "_ps", "_met", "_hpp", "_topmass", "_hdamp", "_qcdscale", "_sel", "_trig", "_tWth"};

bool IsForBreakdown[nUncertainties] = {false, true,         false,     true,    true,              true,              false,             true,       false,   false, true,   false,             true,   true,     false,     false,       false,  false, false, false,  false,   false,   false,   false,   false,  false,      false,        false,         true,       true,     false,     false,  true,  true,   true,       false,    false,       false,  false,   false};
bool IsPtCorrelated[nUncertainties] = {false, false,        true,      true,    false,             true,              false,             false,      true,     true, false,  false,             false,  false,    false,     true,        false,  false, false, false,  false,   false,   false,   false,   false,  false,      false,        false,         true,       false,    false,     false,  false, false,  true,       false,    false,       false,  false,   false};

int MeasurementSystematic[nMeasurements][nUncertainties] = { {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							     {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
							     {1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},    
							     {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0},  
							     {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0} };

// !!! Fraction of events with statistical correlation: without weighting
float FracPTJP[] = {0.096491, 0.128465, 0.162212, 0.182612, 0.203323, 0.235041, 0.239796};
float FracS8JP[] = {0.206224, 0.230399, 0.310995, 0.331386, 0.356070, 0.403861, 0.432426};
float FracPTS8[] = {0.467892, 0.488014, 0.413078, 0.434579, 0.44793,  0.401653, 0.413572};

// TTbar global measurement
TString Chi2Strategy = "Overall";

// This is for avereging the results on ttbar spectrum:
float TTbarPt[nBinsCampaign] = {85., 50., 36., 26., 14., 2.1, 0.};
float TTbarSpectrumScale;
