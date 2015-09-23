#include "Measurements/Run2015B/CampaignInfo.h"

bool System8  = true;
bool Ptrel    = true;
bool Ip3d     = false;
bool LifeTime = true;
bool JPsi     = false;
bool TTbar    = false;
bool KinTT    = false;
bool TagCntTT = true;

bool etamin = false;
bool etamax = false;

bool PrintComparison = false;

bool StatisticalCorrelation = true;

bool InflateStatistic = false;

bool OfficialCMS = false;

bool SystematicBreakdown = false;

enum { NBINS = 20 };

float SystematicCorrelation[nUncertainties][NBINS][nMeasurements*(nMeasurements-1)/2];

float MeasuredScaleFactor[nMeasurements][NBINS], MeasuredScaleFactorUncertainty[nMeasurements][NBINS][nUncertainties];

float xpt[nMeasurements][NBINS], expt[nMeasurements][NBINS]; 
float MeasuredScaleFactorValue[nMeasurements][NBINS], MeasuredScaleFactorError[nMeasurements][NBINS], MeasuredScaleFactorStatistic[nMeasurements][NBINS];

//int MeasurementBinIndex[nBinsCampaign][nMeasurements];
//int MeasurementsForBin[nBinsCampaign];
int nTotalMeasurements;
int MeasurementBinIndex[1000];
int CampaignBinIndex[1000];

TString TSfun1;

float sf[NBINS], sf_stat[NBINS], sf_syst[NBINS], sf_eror[NBINS], sf_uncerbreak[NBINS][50];
float fun_val[NBINS], fun_err[NBINS], fun_sys[NBINS];

double Chi2Normal[NBINS];
int nBinMeasurements[NBINS];
