// General data info
TString CampaignName = "CTAGMoriond2017NoBins";
string CampaignNameString = "CTAGMoriond2017NoBins";
TString CampaignLuminosity = "36.7 fb^{-1}"; 
TString CenterOfMassEnergy = " (13 TeV, 25 ns)";

// Taggers
int const nTaggingAlgorithms = 3;
string TaggingAlgorithmName[nTaggingAlgorithms] = {"ctag", "csv", "cmva"};

// Pt bins
int const nBinsCampaign = 1;
float MinPtCampaign, MaxPtCampaign;

float xPt[nBinsCampaign]  = {100}; // {30.0, 42.5, 60.0, 135}; //bin center
float exPt[nBinsCampaign] = {100}; // {5, 7.5, 10, 65};        //bin width

// Measurements
int const nMeasurements = 2;
string MeasurementName[nMeasurements] = {"ttbar", "wpluscharm"};
string MeasurementFlag[nMeasurements] = {"ttbar", "wpluscharm"};
string MeasurementPlot[nMeasurements] = {"_ttj",   "_wcc"};
bool UseThisMeasurement[nTaggingAlgorithms][nMeasurements] = { 
	{true, true}, 
	{true, true}, 
	{true, true}, 
};

// Plot parameters
int GraphXBins [nMeasurements] = {7    ,      4};  //, 3, 7, 9, 6, 5, 6};
int GraphStyle [nMeasurements] = {3005 ,     -1};  //, -1, -1, -1, -1, -1, -1};
int GraphColor [nMeasurements] = {kBlue,   kRed};  //, kBlue, kRed, kOrange, kBlack, 28, kBlack};
int GraphMarker[nMeasurements] = {22   ,     29};  //, 23, 20, 33, 21, 25, 21};
float GraphSize[nMeasurements] = {1.3  ,    1.3};  //, 1.3, 1.3, 1., 1., 1., 1.};
int GraphWidth [nMeasurements] = {6    ,      6};  //, 4, 4, 4, 4, 4, 4};
int MeasurementOrder[nMeasurements] = {0, 1};  //, 2, 3, 5, 0, 6};

// Uncertainties
int const nUncertainties = 10;
string UncertaintyName[nUncertainties] = {""   , "_statistic", "_pu", "_brfrag", "_jes", "_ntrack", "_osss", "_zid", "_other", "_sl"};
bool IsForBreakdown[nUncertainties]    = {false, true        , true , true     , true  , true     , true   , true  , true    , true};
bool IsPtCorrelated[nUncertainties]    = {false, false       , true , true     , true  , true     , true   , true  , true    , true};

int MeasurementSystematic[nMeasurements][nUncertainties] = { 
  {1, 1, 1, 0, 1, 0, 0, 0, 1, 0},
  {1, 1, 1, 1, 1, 1, 1, 1, 0, 1},
};

// !!! Fraction of events with statistical correlation: without weighting
float FracPTJP[] = {0.096491, 0.128465, 0.162212, 0.182612, 0.203323, 0.235041, 0.239796};
float FracS8JP[] = {0.206224, 0.230399, 0.310995, 0.331386, 0.356070, 0.403861, 0.432426};
float FracPTS8[] = {0.467892, 0.488014, 0.413078, 0.434579, 0.44793,  0.401653, 0.413572};

// TTbar global measurement
float TTbarScale = 1.;
TString Chi2Strategy = "Overall";

// This is for avereging the results on ttbar spectrum:
string WeightedMeasurement = "ttbar";
float TTbarPt[nBinsCampaign] = {1.};///*{0.01, 0.97, 0.01, 0.01};///*/{3962, 5436, 5025, 6405};
float TTbarSpectrumScale = 1;//*/3962+ 5436+ 5025+ 6405;

//fit function for plots
TString FittingFunction = "[0]";
