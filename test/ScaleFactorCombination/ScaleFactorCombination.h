#include "Measurements/Run201576X/CampaignInfo.h"

bool PrintPNG = true;
bool PrintPDF = true;
bool PrintC   = false;

bool etamin = false;
bool etamax = false;

bool PrintComparison        = false;
bool StatisticalCorrelation = true;
bool InflateStatistic       = false;
bool SystematicBreakdown    = false;

int TaggingAlgorithm = -1;

enum { NBINS = 20 };

float SystematicCorrelation[nUncertainties][NBINS][nMeasurements*(nMeasurements-1)/2];

float MeasuredScaleFactor[nMeasurements][NBINS], MeasuredScaleFactorUncertainty[nMeasurements][NBINS][nUncertainties];

float xpt[nMeasurements][NBINS], expt[nMeasurements][NBINS]; 
float MeasuredScaleFactorValue[nMeasurements][NBINS], MeasuredScaleFactorError[nMeasurements][NBINS], MeasuredScaleFactorStatistic[nMeasurements][NBINS];

int nTotalMeasurements;
int MeasurementBinIndex[1000];
int CampaignBinIndex[1000];

TString TSfun1;

float sf[NBINS], sf_stat[NBINS], sf_error[NBINS], sf_uncerbreak[NBINS][50];
float fun_val[NBINS], fun_err[NBINS], fun_unc[NBINS][50];

double NormalizedChi2;

