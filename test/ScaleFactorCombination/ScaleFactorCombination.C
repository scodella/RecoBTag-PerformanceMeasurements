#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TText.h>
#include <TString.h>
#include <TLatex.h>

#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "ScaleFactorCombination.h"
#include "BTagCalibrationStandalone.cc"

float TakeMaximum(float up, float down) {

  float ThisMaximum = fabs(up);
  if (fabs(down)>ThisMaximum) ThisMaximum = fabs(down);
  if (up<0.) ThisMaximum *= -1.;

  return ThisMaximum;

}

float TakeAverage(float up, float down) {

  float ThisAverage = (fabs(up) + fabs(down))/2.;
  if (up<0.) ThisAverage *= -1.;

  return ThisAverage;

}

#include "Campaigns/CombinedScaleFactors_7TeVLegacyNoTTbar.C"

TMatrixD BuildCovarianceMatrix(TString ErrorCategory = "") {

  TMatrixD MatCov(nTotalMeasurements, nTotalMeasurements);

  MatCov *= 0.;

  int RowIndex = 0;
  for (int im1 = 0; im1<nTotalMeasurements; im1++) {

      int ColumnIndex = 0;
      for (int im2 = 0; im2<nTotalMeasurements; im2++) {

	int jm1 = MeasurementBinIndex[im1];
	int jm2 = MeasurementBinIndex[im2];
	
	int bpt1 = CampaignBinIndex[im1];
	int bpt2 = CampaignBinIndex[im2];
	
	// This term contains systematics specific of a method, which are not correlated with other methods, for which we can ignore the bin-to-bin correlation, and for which we do not need to compute the contribution to the total uncertainty
	float SpecificUncertainty = 0.;
	if (jm1==jm2) SpecificUncertainty = pow(MeasuredScaleFactorUncertainty[jm1][bpt1][0], 2);
	
	for (int is = 1; is<nUncertainties; is++) {
	  
	  // "" if we are making the fit
	  // Specifying a systematic to compute its contribution to the total uncertatiny
	  if (ErrorCategory=="" || ErrorCategory==UncertaintyName[is]) {
	    
	    float ThisMatrixElement = 0.;
	    
	    float ThisMatrixTerm = MeasuredScaleFactorUncertainty[jm1][bpt1][is]*MeasuredScaleFactorUncertainty[jm2][bpt2][is];
	    
	    if (bpt1==bpt2 // Always correlate a systematic in the same pt bin
		|| IsPtCorrelated[is]) { // Systematic to be correlated across the pt bins
	      
	      // !!! Statistical uncertainty needs special treatment for correlation across
	      // the methods
	      if (UncertaintyName[is]=="_statistic") {
		
		if (jm1==jm2) 
		  ThisMatrixElement = fabs(ThisMatrixTerm);
		else {
		  
		  if ((MeasurementName[jm1]=="PtRel" && MeasurementName[jm2]=="System8") ||
		      (MeasurementName[jm2]=="PtRel" && MeasurementName[jm1]=="System8")) {
		    
		    ThisMatrixElement = sqrt(FracPTS8[bpt1])*fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="PtRel" && MeasurementName[jm2]=="IP3D") ||
			     (MeasurementName[jm2]=="PtRel" && MeasurementName[jm1]=="IP3D")) {
		    
		    ThisMatrixElement = fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="PtRel" && MeasurementName[jm2]=="LT") ||
			     (MeasurementName[jm2]=="PtRel" && MeasurementName[jm1]=="LT")) {
		    
		    ThisMatrixElement = sqrt(FracPTS8[bpt1])*FracS8JP[bpt1]*fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="System8" && MeasurementName[jm2]=="IP3D") ||
			     (MeasurementName[jm2]=="System8" && MeasurementName[jm1]=="IP3D")) {
		    
		    ThisMatrixElement = sqrt(FracPTS8[bpt1])*fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="LT" && MeasurementName[jm2]=="IP3D") ||
			     (MeasurementName[jm2]=="LT" && MeasurementName[jm1]=="IP3D")) {
		    
		    ThisMatrixElement = sqrt(FracPTS8[bpt1]*FracS8JP[bpt1])*fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="System8" && MeasurementName[jm2]=="LT") ||
			     (MeasurementName[jm2]=="System8" && MeasurementName[jm1]=="LT")) {
		    
		    ThisMatrixElement = sqrt(FracS8JP[bpt1])*fabs(ThisMatrixTerm);
		    
		  } 
		  
		}
		
		
	      }
	      // Now look at the other systematics. The following list must contains:
	      // - systematics correlated across the pt bin (or between measurements)
	      // - systematics for which we want to provide the specific contribution to the total uncertainty
	      // Any other systematic can be left to the SpecificUncertainty term
	      else if (ThisMatrixTerm!=0. && (jm1!=jm2 || IsPtCorrelated[is] || IsForBreakdown[is])) {

		float Coefficient = 1.;
		if (UncertaintyName[is]=="_ltothers" && bpt1!=bpt2) Coefficient = 0.5; // !!! Special for LT systematics
		
		if (UncertaintyName[is]!="_jetaway" || jm1==jm2) // !!! Special case for jetaway systematic (as in Run1)
		  ThisMatrixElement = Coefficient*ThisMatrixTerm;
		
	      } 
	      
	      MatCov(RowIndex, ColumnIndex) += ThisMatrixElement;
	      
	      if (ThisMatrixElement!=0. && jm1==jm2 && bpt1==bpt2) 
		SpecificUncertainty -= ThisMatrixTerm;
	      
	    }
	    
	  } 
	  
	}
	
	if (ErrorCategory=="")
	  if (SpecificUncertainty>0.00001 && jm1==jm2 && bpt1==bpt2)
	    MatCov(RowIndex, ColumnIndex) += SpecificUncertainty;
	
	ColumnIndex++;
	
      }
      
      RowIndex++;
      
  }
  
  return MatCov;

}

void ScaleFactorsFit() {
  
  cout << "  Build MatU" << endl;

  TMatrixD MatU(nTotalMeasurements, nBinsCampaign);

  MatU *= 0.;

  int RowIndex = 0;
  for (int im = 0; im<nTotalMeasurements; im++) {
    
    int bpt = CampaignBinIndex[im];

    int jm = MeasurementBinIndex[im];
    if (MeasurementName[jm]!="TagCountTT") MatU(RowIndex, bpt) = 1.;
    else {
      
      for (int bpt2 = 0; bpt2<nBinsCampaign; bpt2++) 
	MatU(RowIndex, bpt2) = TTbarPt[bpt2]/TTbarSpectrumScale;

    }
    
    RowIndex++;
    
  }
  
  cout << "  Build MatCov" << endl;

  TMatrixD MatCov = BuildCovarianceMatrix();
  TMatrixD MatCovStatistics = BuildCovarianceMatrix("_statistic");
  
  cout << "  Build VecSF" << endl;
  
  TMatrixD VecSF(nTotalMeasurements, 1);

  VecSF *= 0.;

  RowIndex = 0;
  for (int im = 0; im<nTotalMeasurements; im++) {

    int jm = MeasurementBinIndex[im];
    int bpt = CampaignBinIndex[im];
    
    VecSF(RowIndex, 0) = MeasuredScaleFactor[jm][bpt];
    
    RowIndex++;
    
  }
  
  cout << "  Make fit" << endl;
  
  TMatrixD MatUTr(nBinsCampaign, nTotalMeasurements); 
  MatUTr.Transpose(MatU);
  
  double Determinant;
  TMatrixD MatCovInv(MatCov);
  MatCovInv.Invert(&Determinant);
  if (Determinant==0.)
    cout << "      MatCov not invertible" << endl;
  
  TMatrixD Aux1(nBinsCampaign, nTotalMeasurements);
  Aux1.Mult(MatUTr, MatCovInv);
  
  TMatrixD Aux2(nBinsCampaign, nBinsCampaign);
  Aux2.Mult(Aux1, MatU); 
  Aux2.Invert(&Determinant);

  if (Determinant==0.)
    cout << "    Aux2 not invertible" << endl;
  
  TMatrixD MatCoefficients(nBinsCampaign, nTotalMeasurements);
  MatCoefficients.Mult(Aux2, Aux1);

  cout << "  Get SFs" << endl;

  TMatrixD MatSF(nBinsCampaign, 1);
  MatSF.Mult(MatCoefficients, VecSF); 

  cout << "  Get SF errors" << endl;

  TMatrixD MatCoefficientsTr(nTotalMeasurements, nBinsCampaign);
  MatCoefficientsTr.Transpose(MatCoefficients);

  TMatrixD Aux3(nBinsCampaign, nTotalMeasurements);
  Aux3.Mult(MatCoefficients, MatCov);

  TMatrixD MatSFErr(nBinsCampaign, nBinsCampaign);
  MatSFErr.Mult(Aux3, MatCoefficientsTr);

  TMatrixD Aux3Statistics(nBinsCampaign, nTotalMeasurements);
  Aux3Statistics.Mult(MatCoefficients, MatCovStatistics);

  TMatrixD MatSFErrStatistics(nBinsCampaign, nBinsCampaign);
  MatSFErrStatistics.Mult(Aux3Statistics, MatCoefficientsTr);

  cout << "  Fill fit results" << endl;

  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

    sf[bpt] = MatSF(bpt, 0);

    sf_error[bpt] = sqrt(MatSFErr(bpt, bpt));
    sf_stat[bpt] = sqrt(MatSFErrStatistics(bpt, bpt));

    sf_uncerbreak[bpt][0] = sf_error[bpt];
    sf_uncerbreak[bpt][1] = sf_stat[bpt];

  }

  if (SystematicBreakdown) {

    for (int is = 2; is<nUncertainties; is++) 
      if (IsForBreakdown[is]) {
	
	TMatrixD MatCovSystematic = BuildCovarianceMatrix(UncertaintyName[is]);
	
	TMatrixD Aux3Systematic(nBinsCampaign, nTotalMeasurements);
	Aux3Systematic.Mult(MatCoefficients, MatCovSystematic);
	
	TMatrixD MatSFErrSystematic(nBinsCampaign, nBinsCampaign);
	MatSFErrSystematic.Mult(Aux3Systematic, MatCoefficientsTr);
	
	for (int bpt = 0; bpt<nBinsCampaign; bpt++)
	  sf_uncerbreak[bpt][is] = sqrt(MatSFErrSystematic(bpt, bpt));
	
      }
    
  }

  cout << "  Compute Chi2" << endl;

  NormalizedChi2 = 0.;
  for (int im1 = 0; im1<nTotalMeasurements; im1++)
    for (int im2 = 0; im2<nTotalMeasurements; im2++) 
      for (int bpt1 = 0; bpt1<nBinsCampaign; bpt1++) 
	for (int bpt2 = 0; bpt2<nBinsCampaign; bpt2++) {
	  
	  int jm1 = MeasurementBinIndex[im1];
	  int jm2 = MeasurementBinIndex[im2];
	  
	  if (MeasurementName[jm1]!="TagCountTT" && MeasurementName[jm2]!="TagCountTT") {
	    
	    if (bpt1==CampaignBinIndex[im1] && bpt2==CampaignBinIndex[im2]) 
	      NormalizedChi2 += MatU(im1, bpt1)*(MeasuredScaleFactor[jm1][bpt1] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, bpt2)*(MeasuredScaleFactor[jm2][bpt2] - sf[bpt2]);
	    
	  } else if (MeasurementName[jm1]=="TagCountTT" && MeasurementName[jm2]=="TagCountTT") {
	    NormalizedChi2 += MatU(im1, bpt1)*(MeasuredScaleFactor[jm1][0] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, bpt2)*(MeasuredScaleFactor[jm2][0] - sf[bpt2]);
	    
	  }
	  
	}
  
  NormalizedChi2 /= (nTotalMeasurements - nBinsCampaign);
  cout << "Normalized Chi2 = " << NormalizedChi2 << endl;
  
}

TGraphErrors *MakeTGraphErrors(int nBins, float VectX[NBINS], float VectY[NBINS], float VectXError[NBINS], float VectYError[NBINS], int ThisStyle, int ThisColor, int ThisMarker, float ThisSize, int ThisWidth) {
  
  TGraphErrors *ThisTGraphErrors = new TGraphErrors(nBins, VectX, VectY, VectXError, VectYError);

  if (ThisStyle>=0) ThisTGraphErrors->SetFillStyle(ThisStyle);
  ThisTGraphErrors->SetFillColor(ThisColor);
  ThisTGraphErrors->SetMarkerColor(ThisColor);
  ThisTGraphErrors->SetLineColor(ThisColor);
  ThisTGraphErrors->SetLineStyle(1);
  ThisTGraphErrors->SetMarkerStyle(ThisMarker);
  ThisTGraphErrors->SetMarkerSize(ThisSize);
  ThisTGraphErrors->SetLineWidth(ThisWidth);

  return ThisTGraphErrors;

}

void ReadMeasurements(string  BTagger, TString OP) {
 
  MinPtCampaign = xPt[0] - exPt[0];
  MaxPtCampaign = xPt[nBinsCampaign-1] + exPt[nBinsCampaign-1];

  BTagEntry::OperatingPoint WP;
  if (OP=="Loose") WP = BTagEntry::OP_LOOSE;
  if (OP=="Medium") WP = BTagEntry::OP_MEDIUM;
  if (OP=="Tight") WP = BTagEntry::OP_TIGHT;

  BTagEntry::JetFlavor JF = BTagEntry::FLAV_B;

  nTotalMeasurements = 0;

  for (int nm = 0; nm<nMeasurements; nm++) 
    if (UseThisMeasurement[TaggingAlgorithm][nm]) {

      string CSVFileName = "./Measurements/" + CampaignNameString + "/" + BTagger + "_" + MeasurementName[nm] + ".csv";
      BTagCalibration calib(BTagger, CSVFileName);

      BTagCalibrationReader *CentralReader = new BTagCalibrationReader(&calib, WP, MeasurementFlag[nm], "central");

      int nMeasBins = -1; float LastMeasuredSF = 0., LowBinEdge = 0.;

      for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

        MeasuredScaleFactor[nm][bpt] = CentralReader->eval(JF, 0., xPt[bpt]);
	
	if (MeasuredScaleFactor[nm][bpt]!=LastMeasuredSF) {
	  
	  if (nMeasBins>=0) {
	    
	    xpt[nm][nMeasBins] = (LowBinEdge + xPt[bpt] - exPt[bpt])/2.;
	    expt[nm][nMeasBins] = xpt[nm][nMeasBins] - LowBinEdge;
	    MeasuredScaleFactorValue[nm][nMeasBins] = LastMeasuredSF;
	    
	  }
	  
	  nMeasBins++;

	  if (MeasuredScaleFactor[nm][bpt]!=0.) {
	
	    MeasurementBinIndex[nTotalMeasurements] = nm;
	    CampaignBinIndex[nTotalMeasurements] = bpt;
	    nTotalMeasurements++;
	    
	  }
	  
	  LowBinEdge = xPt[bpt] - exPt[bpt];

	  LastMeasuredSF = MeasuredScaleFactor[nm][bpt];
	  
	}      

      }

      if (LastMeasuredSF>0.) {

	xpt[nm][nMeasBins] = (LowBinEdge + xPt[nBinsCampaign-1] + exPt[nBinsCampaign-1])/2.;
	expt[nm][nMeasBins] = xpt[nm][nMeasBins] - LowBinEdge;
	MeasuredScaleFactorValue[nm][nMeasBins] = LastMeasuredSF;

	nMeasBins++;

      }

      float TotalErrorForFit[50], TotalCSVErrorForFit[50];
      for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
	
	TotalErrorForFit[bpt] = 0.;
	TotalCSVErrorForFit[bpt] = 0.;
	
      }
      
      float TotalError[50], TotalCSVError[50];
      for (int bpt = 0; bpt<nMeasBins; bpt++) {
	
	TotalError[bpt] = 0.;
	TotalCSVError[bpt] = 0.;
	
      }
      
      for (int is = 0; is<nUncertainties; is++) { 

	for (int bpt = 0; bpt<nBinsCampaign; bpt++)
	  MeasuredScaleFactorUncertainty[nm][bpt][is] = 0.;
	
	if (MeasurementSystematic[nm][is]!=0) {
	  
	  string UpFlag = "up" + UncertaintyName[is], DownFlag = "down" + UncertaintyName[is];
	  BTagCalibrationReader *SystUpReader = new BTagCalibrationReader(&calib, WP, MeasurementFlag[nm], UpFlag);
	  BTagCalibrationReader *SystDownReader = new BTagCalibrationReader(&calib, WP, MeasurementFlag[nm], DownFlag);
	  
	  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
	    
	    float UpSF = SystUpReader->eval(JF, 0., xPt[bpt]); UpSF -= MeasuredScaleFactor[nm][bpt];
	    float DownSF = SystDownReader->eval(JF, 0., xPt[bpt]); DownSF -= MeasuredScaleFactor[nm][bpt];
	    
	    float ThisError =  TakeAverage(UpSF, DownSF);
	    MeasuredScaleFactorUncertainty[nm][bpt][is] = ThisError;
	    if (UncertaintyName[is]=="") TotalCSVErrorForFit[bpt] = ThisError;
	    else TotalErrorForFit[bpt] += ThisError*ThisError;
	    
	  }
	  
	  for (int bpt = 0; bpt<nMeasBins; bpt++) {
	    
	    float UpSF = SystUpReader->eval(JF, 0., xpt[nm][bpt]); UpSF -= MeasuredScaleFactorValue[nm][bpt];
	    float DownSF = SystDownReader->eval(JF, 0., xpt[nm][bpt]); DownSF -= MeasuredScaleFactorValue[nm][bpt];
	   
	    float ThisError = TakeAverage(UpSF, DownSF);
	    if (UncertaintyName[is]=="_statistic") MeasuredScaleFactorStatistic[nm][bpt] = ThisError;
	    if (UncertaintyName[is]=="") TotalCSVError[bpt] = ThisError;
	    else TotalError[bpt] += ThisError*ThisError;
	    
	  }

	}
	
      }

      for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
	
	if (TotalCSVErrorForFit[bpt]!=0.) 
	  MeasuredScaleFactorUncertainty[nm][bpt][0] = TotalCSVErrorForFit[bpt];
	else 
	  MeasuredScaleFactorUncertainty[nm][bpt][0] = sqrt(TotalErrorForFit[bpt]);
	
      }

      for (int bpt = 0; bpt<nMeasBins; bpt++) {
	
	if (TotalCSVError[bpt]!=0.) 
	  MeasuredScaleFactorError[nm][bpt] = TotalCSVError[bpt];
	else 
	  MeasuredScaleFactorError[nm][bpt] = sqrt(TotalError[bpt]);
	
      }
      
      if (MeasurementName[nm]=="TagCountTT" || MeasurementName[nm]=="KIN") 
	for (int cpt = 0; cpt<nBinsCampaign; cpt++) 
	  MeasuredScaleFactorUncertainty[nm][cpt][1] *= TTbarScale;
      
      if (InflateStatistic) { // Run1 stuff ...

	for (int bpt = 0; bpt<nMeasBins; bpt++) {

	  int BinMultiplicity = 0;
	  for (int cpt = 0; cpt<nBinsCampaign; cpt++) 
	    if (xPt[cpt]>xpt[nm][bpt]-expt[nm][bpt] && xPt[cpt]<xpt[nm][bpt]+expt[nm][bpt])
	      BinMultiplicity++;
  
	  for (int cpt = 0; cpt<nBinsCampaign; cpt++) 
	    if (xPt[cpt]>xpt[nm][bpt]-expt[nm][bpt] && xPt[cpt]<xpt[nm][bpt]+expt[nm][bpt]) 
	      MeasuredScaleFactorUncertainty[nm][cpt][1] *= sqrt(BinMultiplicity);

	}

      }

      for (int bpt = 0; bpt<nMeasBins; bpt++) {

	float xPtShift = 1.;
	if (expt[nm][bpt]>=15.) xPtShift = 2.;
	if (nm==0) xPtShift *= -1.;
	if (nm==3 || xpt[nm][bpt]>160.) xPtShift = 0.;
	xpt[nm][bpt] += xPtShift;
 
      }
      
    }
  
  // Some old stuff ...
  /*
  for (int is = 2; is<nUncertainties; is++) {
    
    int ThisCounter = 0;
    
    for (int im = 0; im<nMeasurements; im++) 
      for (int jm = im+1; jm<nMeasurements; jm++) {

	for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
	  
	  SystematicCorrelation[is][bpt][ThisCounter] = 0.;
	  
	  if (UseThisMeasurement[TaggingAlgorithm][im] && UseThisMeasurement[TaggingAlgorithm][jm] && MeasurementSystematic[im][is]!=0 && MeasurementSystematic[jm][is]!=0) {
	    SystematicCorrelation[is][bpt][ThisCounter] = MeasuredScaleFactorUncertainty[im][bpt][is]*MeasuredScaleFactorUncertainty[jm][bpt][is]/
	      fabs(MeasuredScaleFactorUncertainty[im][bpt][is]*MeasuredScaleFactorUncertainty[jm][bpt][is]);
	    
	  }
	  
	}

	ThisCounter++;

      }
  
  }
  */

}
  
void SetAlgorithmToUse(string BTagger) {

  TaggingAlgorithm = -1;

  for (int tg = 0; tg<nTaggingAlgorithms; tg++) 
    if (BTagger==TaggingAlgorithmName[tg]) TaggingAlgorithm = tg;
  
}

TStyle* PlotStyle() {

  TStyle *plotStyle = new TStyle("PLOT","Plot style CMS");
  
  plotStyle->SetOptTitle(0);  
  plotStyle->SetOptStat(0); 
  plotStyle->SetOptFit(0);
  
  plotStyle->SetPaperSize(20.,26.);
  
  plotStyle->SetEndErrorSize(2);
  plotStyle->SetErrorX(0.);
  
  plotStyle->SetFrameBorderMode(0);
  plotStyle->SetFrameFillColor(0);
  plotStyle->SetFrameFillStyle(0);
  plotStyle->SetCanvasBorderMode(0);
  plotStyle->SetFillColor(0);
  plotStyle->SetCanvasColor(0);
  plotStyle->SetCanvasBorderSize(2);
  
  plotStyle->SetPadBorderMode(0);
  plotStyle->SetPadColor(0);
  plotStyle->SetPadGridX(false);
  plotStyle->SetPadGridY(false);
  plotStyle->SetGridColor(0);
  plotStyle->SetGridStyle(3);
  plotStyle->SetGridWidth(1);
  
  plotStyle->SetCanvasDefX(0);
  plotStyle->SetCanvasDefY(0);
  
  plotStyle->SetHistLineColor(1);
  plotStyle->SetHistLineStyle(0);
  plotStyle->SetHistLineWidth(1);
  
  plotStyle->SetPadTickX(1);
  plotStyle->SetPadTickY(1);
  
  plotStyle->SetPadLeftMargin(0.16);
  plotStyle->SetPadRightMargin(0.02);
  plotStyle->SetPadTopMargin(0.06);
  plotStyle->SetPadBottomMargin(0.13);
  
  plotStyle->SetLabelFont(42,"x");
  plotStyle->SetTitleFont(42,"x");
  plotStyle->SetLabelFont(42,"y");
  plotStyle->SetTitleFont(42,"y");
  plotStyle->SetLabelFont(42,"z");
  plotStyle->SetTitleFont(42,"z");
  
  plotStyle->SetLabelSize(0.05,"x");
  plotStyle->SetTitleSize(0.06,"x");
  plotStyle->SetLabelSize(0.05,"y");
  plotStyle->SetTitleSize(0.06,"y");
  plotStyle->SetLabelSize(0.035,"z");
  plotStyle->SetTitleSize(0.035,"z");
  
  plotStyle->SetLabelOffset(0.007,"XYZ");
  
  plotStyle->SetNdivisions(510,"XYZ");
  plotStyle->SetStripDecimals(kTRUE);
  plotStyle->SetTickLength(0.03, "XYZ");
  
  plotStyle->SetTitleXOffset(0.9);
  plotStyle->SetTitleYOffset(1.25);
  
  return plotStyle;

}

void ScaleFactorCombination(string BTagger, TString OP) {
  
  // @@@ Build the canvas and set the plot style
  
  TCanvas *c1 = new TCanvas("c1", "plots",200,0,700,700);
  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2);
  
  TPad* pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.98,0.98,21);
  TPad* pad2 = new TPad("pad2","This is pad2",0.02,0.03,0.98,0.49,21);

  // Run2015B Setting
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1->Range(0,0,1,1);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.16);
  c1->SetRightMargin(0.02);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.13);
  c1->SetFrameFillColor(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
  //pad1->SetLogy();
  pad1->SetLogx();
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetLeftMargin(0.16);
  pad1->SetRightMargin(0.02);
  pad1->SetTopMargin(0.065);
  pad1->SetBottomMargin(0.13);
  pad1->SetFrameFillStyle(0);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameFillStyle(0);
  pad1->SetFrameBorderMode(0);
  pad1->Draw();

  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(2);
  //pad2->SetGridy();
  pad2->SetLogx();
  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetLeftMargin(0.16);
  pad2->SetRightMargin(0.02);
  //pad2->SetTopMargin(0.05);
  //pad2->SetBottomMargin(0.31);
  pad2->SetTopMargin(0.065);
  pad2->SetBottomMargin(0.13);
  pad2->SetFrameFillStyle(0);
  pad2->SetFrameBorderMode(0);
  pad2->SetFrameFillStyle(0);
  pad2->SetFrameBorderMode(0);
  pad2->Draw();
  // End Run2015B Setting 

  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");

  /* Kirill's Style
  static TStyle* plotStyle = PlotStyle(); 
  gROOT->SetStyle("PLOT");   
  gROOT->ForceStyle(); 

  pad1->Draw();
  pad2->Draw();
  */// End Kirill's Style 

  // @@@ Define some general parameter and read the measurements
  
  SetAlgorithmToUse(BTagger);

  TString title = BTagger;
  if (OP=="Loose")  title += "L"; 
  if (OP=="Medium") title += "M"; 
  if (OP=="Tight")  title += "T"; 

  float exstat[NBINS];
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) exstat[bpt] = 0.0001;
    
  if (!StatisticalCorrelation) 
    for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
      FracPTJP[bpt] = 1.; 
      FracS8JP[bpt] = 1.; 
      FracPTS8[bpt] = 1.;
    }
  
  ReadMeasurements(BTagger, OP);

  // @@@ Plot the measurements in the top pad

  pad1->cd();
  
  float PlotMaxPtCampaign = MaxPtCampaign;
  //if (OP=="Tight") PlotMaxPtCampaign = 200.;
  //if (OP=="Medium") PlotMaxPtCampaign = 300.;
  float YAxisWidth = 0.3;
  if (OP=="Tight") YAxisWidth = 0.4;
  TH2F* histo = new TH2F("histo","",58,MinPtCampaign,PlotMaxPtCampaign,100,1.-YAxisWidth,1.+YAxisWidth);
  
  histo->Draw(""); 
  histo->SetLabelSize(0.05, "XYZ");
  histo->SetTitleSize(0.06, "XYZ"); 
  histo->SetLabelFont(42, "XYZ"); 
  histo->SetTitleFont(42, "XYZ");
  histo->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  histo->GetYaxis()->SetTitle("Data/Simulation SF_{b}");
  //histo->SetTitleOffset(1.1,"X"); // Ideal for .png
  histo->SetTitleOffset(0.95,"X");
  histo->SetTitleOffset(0.8,"Y");
  histo->SetTickLength(0.06,"X");
  histo->SetNdivisions(509, "XYZ");
  histo->GetXaxis()->SetMoreLogLabels();
  histo->GetXaxis()->SetNoExponent();
  
  TGraphErrors *ScaleFactorStatistic[nMeasurements];
  TGraphErrors *ScaleFactorError[nMeasurements];

  for (int nm = 0; nm<nMeasurements; nm++) {

    ScaleFactorStatistic[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], exstat, MeasuredScaleFactorStatistic[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], GraphWidth[nm]);

    ScaleFactorError[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], expt[nm], MeasuredScaleFactorError[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], 2);

    if (UseThisMeasurement[TaggingAlgorithm][nm]) { 
      
      ScaleFactorStatistic[nm]->Draw("P"); 
      ScaleFactorError[nm]->Draw("P"); 

    }
    
  }
  
  // Legend for the top pad
  TString LegTitle = title;

  float LegXOffset = 0.;
  if (OP=="Medium") LegXOffset = 0.06;
  TLegend *leg = new TLegend(0.48+LegXOffset,0.66,0.70+LegXOffset,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(62);
  leg->SetTextSize(0.05);   
  leg->SetHeader(LegTitle);
  for (int nm = 0; nm<nMeasurements; nm++) 
    if (UseThisMeasurement[TaggingAlgorithm][nm]) {
      TString LegText = MeasurementName[nm];
      leg->AddEntry(ScaleFactorError[nm], LegText, "PL");
    }
  leg->Draw();

  // Run2015B Style
  TLatex *tex = new TLatex(0.2,0.88,"CMS"); 
  tex->SetNDC(); 
  tex->SetTextAlign(13);
  tex->SetTextFont(61);
  tex->SetTextSize(0.07475);
  tex->SetLineWidth(2); 
  tex->Draw();                                                                                       
  TLatex *tex2 = new TLatex(0.2,0.79,"Preliminary"); 
  tex2->SetNDC();
  tex2->SetTextAlign(13);
  tex2->SetTextFont(52);
  tex2->SetTextSize(0.05681);
  tex2->SetLineWidth(2);   
  tex2->Draw();   
 
  TLatex *text1 = new TLatex(0.98,0.95125, CampaignLuminosity + CenterOfMassEnergy); 
  text1->SetNDC();                                              
  text1->SetTextAlign(31);                          
  text1->SetTextFont(42);    
  text1->SetTextSize(0.04875);   
  text1->SetLineWidth(2);    
  text1->Draw(); 
  // End Run2015B Style

  // @@@ Performe the combination fit

  std::pair<double,double> combResult;
  
  cout << "NOW STARTING THE FITS " << endl;

  ScaleFactorsFit();

  cout << "FITS DONE, NOW MAKING PLOTS" << endl;

  // Print the combination result in the top pad
  TGraphErrors* sf0 = new TGraphErrors(nBinsCampaign,xPt,sf,exPt,sf_error);
  sf0->SetFillStyle(3005);
  sf0->SetFillColor(kGray+3);
  sf0->Draw("e2"); // to plot fit band
  
  // Superimpose the single measurements on the combination
  for (int nm = 0; nm<nMeasurements; nm++)
    if (UseThisMeasurement[TaggingAlgorithm][MeasurementOrder[nm]]) {
      ScaleFactorStatistic[MeasurementOrder[nm]]->Draw("P");
      ScaleFactorError[MeasurementOrder[nm]]->Draw("P"); 
    }
  
  // @@@ Now plotting the combination in the bottom pad
  
  pad2->cd();
  
  histo = new TH2F("histo","",58,MinPtCampaign,PlotMaxPtCampaign,100,1.-YAxisWidth,1.+YAxisWidth);

  histo->Draw(""); 
  histo->SetLabelSize(0.05, "XYZ");
  histo->SetTitleSize(0.06, "XYZ"); 
  histo->SetLabelFont(42, "XYZ"); 
  histo->SetTitleFont(42, "XYZ");
  histo->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  histo->GetYaxis()->SetTitle("Data/Simulation SF_{b}");
  //histo->SetTitleOffset(1.1,"X"); // Ideal for .png
  histo->SetTitleOffset(0.95,"X"); // 
  histo->SetTitleOffset(0.8,"Y");
  histo->SetTickLength(0.06,"X");
  histo->SetNdivisions(509, "XYZ");
  histo->GetXaxis()->SetMoreLogLabels();
  histo->GetXaxis()->SetNoExponent();
 
  // Fill vector with results from previous campaign in case a comparison is wished
  float SFb_Comp[nBinsCampaign], SFb_Comp_error[nBinsCampaign];
 
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

    double xx = xPt[bpt];

    if ( title == "CSVL" ) {
      SFb_Comp[bpt]      = funSFb_Comp_CSVL(xx);
      SFb_Comp_error[bpt] = SFb_Comp_error_CSVL[bpt];
    }
    if ( title == "CSVM" ) { 
      SFb_Comp[bpt]      = funSFb_Comp_CSVM(xx);
      SFb_Comp_error[bpt] = SFb_Comp_error_CSVM[bpt];
    }
    if ( title == "CSVT" ) { 
      SFb_Comp[bpt]      = funSFb_Comp_CSVT(xx);
      SFb_Comp_error[bpt] = SFb_Comp_error_CSVT[bpt];
    }
    if ( title == "TCHPT" ) { 
      SFb_Comp[bpt]      = funSFb_Comp_TCHPT(xx);
      SFb_Comp_error[bpt] = SFb_Comp_error_TCHPT[bpt];
    }
    
  }
  
  TGraphErrors* sf_Comp = new TGraphErrors(nBinsCampaign,xPt,SFb_Comp,exPt,SFb_Comp_error);

  sf_Comp->SetFillColor(kYellow-9);
  if (PrintComparison) sf_Comp->Draw("e2");
 
  // Superimpose the combination result
  //sf0->SetFillStyle(3005);
  //sf0->SetFillColor(kBlack);
  sf0->Draw("e2"); 


  
  // Define fitting function
  TString FittingFunction;
  if (CampaignName=="Winter13") {
    if ( title == "CSVL") FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    else if ( title == "CSVM" || title == "CSVT") 
      FittingFunction = "[0]+[1]*x+[2]*x*x";
  } else if (CampaignName=="7TeVLegacy") {
    if ( title == "CSVL") FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    else if ( title == "CSVM" ) 
      FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    else if ( title == "CSVT" ) 
      FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    else if ( title == "TCHPT" ) 
      FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
  } else if (CampaignName=="Run2015B") {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
  } else if (CampaignName=="Run201525ns") {
    if (title.Contains("CSVv2") || title!="JPL")
      FittingFunction = "[0]+[1]*log(x+[3])*log(x+[3])*(3-[2]*log(x+[3]))";
    else 
      FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)"; 
  } else if (CampaignName=="Run201576X") {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)"; 
  } 

  TF1* fun1 = new TF1("fun1", FittingFunction, MinPtCampaign, MaxPtCampaign);

  // Initialise some parameter in the fitting function
  if (CampaignName=="Winter13") {
    if ( title == "CSVM" || title == "CSVT")
      fun1->SetParameters(0.938887,0.00017124,-2.76366e-07);
  } else if (CampaignName=="7TeVLegacy") {
    if ( title == "CSVM") 
      fun1->SetParameters(0.899, 0.0565, 0.0606);
    else if ( title == "CSVT")
      fun1->SetParameters(0.9, 0.07/log(670)/log(670), 2./log(670));    
  } else if (CampaignName=="Run2015B") {
    fun1->SetParameters(0.853, 0.0527, 0.453);
    if (OP!="Loose") {
      fun1->FixParameter(1, 0.);
      fun1->FixParameter(2, 0.);
    }
  } else if (CampaignName=="Run201576X") {
    fun1->SetParameter(0, 0.916293);
    fun1->SetParameter(1, 0.0118628);
    fun1->SetParameter(2, 0.0108104);
  } 

  // Fit the pt dependence of the combination result
  fun1->SetLineColor(kBlack);
  fun1->SetLineWidth(2);
  fun1->SetLineStyle(1);
  float MaxPtFit = MaxPtCampaign;
  //if (OP=="Loose") MaxPtFit = 400.;
  //if (OP=="Medium") MaxPtFit = 400.;
  //if (OP=="Tight") MaxPtFit = 200.;
  sf0->Fit("fun1","0","",MinPtCampaign,MaxPtFit);//MaxPtCampaign);
  sf0->Fit("fun1","rve0","",MinPtCampaign,MaxPtFit);//MaxPtCampaign);
//$$  sf0->Fit("fun1","rvee0"); 
  fun1->Draw("same"); 

  // Compute the function value for the pt bins
  std::cout << std::endl;
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

    fun_val[bpt] = fun1->Integral(xPt[bpt]-1.,xPt[bpt]+1.) / 2.; // Run1 keep
    //fun_val[bpt] = sf[bpt];

    //fun_err[bpt] = fun1->IntegralError(xPt[bpt]-exPt[bpt],xPt[bpt]+exPt[bpt]) / exPt[bpt];
    fun_err[bpt] = sf_error[bpt]*fun_val[bpt]/sf[bpt]; // Run1 keep

    if (SystematicBreakdown)
      for (int uu = 0; uu<nUncertainties; uu++)	 
	fun_unc[bpt][uu] = sf_uncerbreak[bpt][uu]*fun_val[bpt]/sf[bpt];

  }

  // Add the funations values to the bottom pad
  TGraphErrors* fun0 = new TGraphErrors(16,xPt,fun_val,exPt,fun_err);
  fun0->SetMarkerColor(kRed);
  fun0->SetLineColor(kRed);
  fun0->SetLineStyle(1);
  fun0->SetMarkerStyle(24);
  fun0->SetMarkerSize(0.001);
  fun0->SetLineWidth(2);
  fun0->Draw("P"); 
  fun1->SetMarkerColor(kBlack);
  fun1->Draw("same"); 

  // Add legend
  leg = new TLegend(0.48+LegXOffset,0.66,0.70+LegXOffset,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(62);
  leg->SetTextSize(0.05);   
  leg->SetHeader(LegTitle);
  leg->AddEntry(sf0,"weighted average","PF");
  leg->AddEntry(fun1,"fit","L");
  leg->AddEntry(fun0,"fit #pm (stat #oplus syst)","LE");
  if (PrintComparison) leg->AddEntry(sf_Comp,LegComp,"F");
  leg->Draw();

  /* // Run1 Style
   tex = new TLatex(0.17,1,"CMS Preliminary, " + CampaignLuminosity + "  at #sqrt{s} = " + CenterOfMassEnergy);
   tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetLineWidth(2);
   tex->Draw();
  */ // End Run1 Style

  // Run2015B Style
  tex = new TLatex(0.2,0.88,"CMS"); 
  tex->SetNDC(); 
  tex->SetTextAlign(13);
  tex->SetTextFont(61);
  tex->SetTextSize(0.07475);
  tex->SetLineWidth(2); 
  tex->Draw();                                                                                       
  tex2 = new TLatex(0.2,0.79,"Preliminary"); 
  tex2->SetNDC();
  tex2->SetTextAlign(13);
  tex2->SetTextFont(52);
  tex2->SetTextSize(0.05681);
  tex2->SetLineWidth(2);   
  tex2->Draw();      

  text1 = new TLatex(0.98,0.95125, CampaignLuminosity + CenterOfMassEnergy); 
  text1->SetNDC();                                              
  text1->SetTextAlign(31);                          
  text1->SetTextFont(42);    
  text1->SetTextSize(0.04875);   
  text1->SetLineWidth(2);    
  text1->Draw();  
  // End Run2015B Style

  // @@@ Print the combination results

  cout << "Print fit results" << endl << endl;
  TSfun1 = fun1->GetExpFormula("p");
  std::cout << " Tagger: " << title << " within " << MinPtCampaign 
	    << " < pt < " << MaxPtCampaign << " GeV, abs(eta) < 2.4, x = pt" << std::endl;
  std::cout << "  SFb = " << TSfun1 << ";" << std::endl << std::endl;
  std::cout << " Fit details" << std::endl;
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
    std::cout << "   " << xPt[bpt]-exPt[bpt] << "-" << xPt[bpt]+exPt[bpt] << ": SF(bin) = " << sf[bpt] << " +/- " << sf_error[bpt] << "; SF(fun) = " << fun_val[bpt] << " +/- " << fun_err[bpt] << " (" << 100.*fun_err[bpt]/fun_val[bpt] << "%); Stat: " << 100.*sf_stat[bpt]/sf[bpt] << "% (" << 100.*sf_stat[bpt]/sf_error[bpt] << "%)" << std::endl;
  }
  std::cout << std::endl;
  
  cout << " Comparison with average ttbar-based measurements" << endl;
   float ScaleFactorTTbarFun = 0., ScaleFactorTTbarFunError = 0.;
   float ScaleFactorTTbarBin = 0., ScaleFactorTTbarBinError = 0.;
   for (int bpt = 0; bpt < nBinsCampaign; bpt++) {
     ScaleFactorTTbarFun += fun_val[bpt]*TTbarPt[bpt]/TTbarSpectrumScale; 
     ScaleFactorTTbarFunError += fun_err[bpt]*TTbarPt[bpt]/TTbarSpectrumScale;
     ScaleFactorTTbarBin += sf[bpt]*TTbarPt[bpt]/TTbarSpectrumScale; 
     ScaleFactorTTbarBinError += sf_error[bpt]*TTbarPt[bpt]/TTbarSpectrumScale;
   }
   std::cout << "  ttbar-like SF        = " << ScaleFactorTTbarFun << " +/- " << ScaleFactorTTbarFunError << std::endl;
   std::cout << "  ttbar-like SF (bins) = " << ScaleFactorTTbarBin << " +/- " << ScaleFactorTTbarBinError << std::endl << std::endl;

   // @@@ Save the plot
   
   c1->Update();

   TString PlotName = "./Plots/SFb_" + CampaignName + "_" + LegTitle;
   for (int nm = 0; nm<nMeasurements; nm++)
     if (UseThisMeasurement[TaggingAlgorithm][nm]) PlotName += MeasurementPlot[nm];
   if (TTbarScale==100.) {
     PlotName.ReplaceAll("_TT",  "_TTnofit");
     PlotName.ReplaceAll("_KIN", "_KINnofit");
     PlotName.ReplaceAll("_TCT", "_TCTnofit");
   }
   if (PrintPNG) c1->Print(PlotName + ".png");
   if (PrintPDF) c1->Print(PlotName + ".pdf");
   if (PrintC)   c1->Print(PlotName + ".C");
   
}

void StoreScaleFactorCombination(string BTagger, TString OP = "All", string MeasType = "comb", bool StoreSystematicBreakdown = false) {

  BTagCalibration calib(BTagger);

  const int nOperatingPoints = 3;
  string OperatingPoint[nOperatingPoints] = {"Loose", "Medium", "Tight"};
  
  const int nMeasurementTypes = 2;
  string MeasTypeName[nMeasurementTypes] = {"comb", "mujets"};

  SystematicBreakdown = StoreSystematicBreakdown;

  float PtBinEdge[nBinsCampaign+1]; PtBinEdge[0] = xPt[0] - exPt[0];
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) PtBinEdge[bpt+1] = xPt[bpt] + exPt[bpt];
  TH1F *UncertaintyHisto = new TH1F ("UncertaintyHisto", "", nBinsCampaign, PtBinEdge);  
  float DiscMin = 0, DiscMax = 1;
  if (BTagger=="JP") { DiscMin = 0; DiscMax = 5; }

  for (int mt = 0; mt<nMeasurementTypes; mt++) 
    if(MeasType=="All" || MeasType==MeasTypeName[mt]) {

      SetAlgorithmToUse(BTagger);

      // !!! Depende on measurement name in CampaignInfo
      TTbarScale = 1.;
      if (MeasTypeName[mt]=="mujets") {

	TTbarScale = 100.;

	for (int nm = 0; nm<nMeasurements; nm++) 
	if (MeasurementName[nm]!="PtRel" && MeasurementName[nm]!="System8" &&
	    MeasurementName[nm]!="LT") 
	  UseThisMeasurement[TaggingAlgorithm][nm] = false;
	
      }

      for (int op = 0; op<nOperatingPoints; op++) 
	if (OP=="All" || OP.Contains(OperatingPoint[op])) {
		  
	  BTagEntry::OperatingPoint wp = BTagEntry::OP_LOOSE;
	  if (OperatingPoint[op]=="Loose")  wp = BTagEntry::OP_LOOSE;
	  if (OperatingPoint[op]=="Medium") wp = BTagEntry::OP_MEDIUM;
	  if (OperatingPoint[op]=="Tight")  wp = BTagEntry::OP_TIGHT;
	  
	  ScaleFactorCombination(BTagger, OperatingPoint[op]);
	  
	  int maxUncertainty = nUncertainties;
	  if (!StoreSystematicBreakdown) maxUncertainty = 1;
	  
	  for (int ThisFlavour = 4; ThisFlavour<=5; ThisFlavour++) {
	    
	    BTagEntry::JetFlavor jfl;
	    if (ThisFlavour==4) jfl = BTagEntry::FLAV_C;
	    if (ThisFlavour==5) jfl = BTagEntry::FLAV_B;
	    
	    int SysFactor = 6 - ThisFlavour;
	      
	    for (int uu = 0; uu<=maxUncertainty; uu++)
	      if (uu<=1 || IsForBreakdown[uu-1]) 
		for (int vs = 0; vs<2; vs++) {
		  
		  if (uu==0 && vs==1) continue;
		  
		  string SysFlag = "central";
		  
		  TString ThisFunction = TSfun1;
		  
		  if (uu>0) {
		    
		    if (vs==0) { SysFlag = "up"; ThisFunction += " + "; }
		    else if (vs==1) { SysFlag = "down";  ThisFunction += " - "; }
		    
		    UncertaintyHisto->Reset();
		    
		    if (uu==1) {
		      
		      for (int bpt = 0; bpt<nBinsCampaign; bpt++) 
			UncertaintyHisto->SetBinContent(bpt+1, SysFactor*fun_err[bpt]);
		      
		    } else { 
		      
		      SysFlag += UncertaintyName[uu-1]; 
		      for (int bpt = 0; bpt<nBinsCampaign; bpt++) 
			UncertaintyHisto->SetBinContent(bpt+1, SysFactor*fun_unc[bpt][uu-1]);
		      
		    }
		    
		  }
		  
		  /* This is not working yet ...		  
		     BTagEntry::Parameters params(wp, MeasTypeName[mt], SysFlag, jfl, -2.4, 2.4, MinPtCampaign, MaxPtCampaign, 0, 1);
		     
		     if (uu>0) {
		     
		     BTagEntry e2(UncertaintyHisto, params); 
		     ThisFunction += e2.formula;
		     cout << e2.formula << endl;
		     
		     }
		     cout << ThisFunction << endl;
		     const TF1 SFFun("SFFun", ThisFunction, MinPtCampaign, MaxPtCampaign);
		     BTagEntry e1(&SFFun, params);  
		     
		     calib.addEntry(e1);
		  */
		  
		  // ... so doing like this for the time being
		  if (uu==0) {
		    
		    BTagEntry::Parameters params(wp, MeasTypeName[mt], SysFlag, jfl, -2.4, 2.4, MinPtCampaign, MaxPtCampaign, DiscMin, DiscMax);
		    
		    const TF1 SFFun("SFFun", ThisFunction, MinPtCampaign, MaxPtCampaign);
		    BTagEntry e1(&SFFun, params);  
		    
		    calib.addEntry(e1);
		    
		  } else {
		    
		    for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
		      
		      BTagEntry::Parameters params(wp, MeasTypeName[mt], SysFlag, jfl, -2.4, 2.4, PtBinEdge[bpt], PtBinEdge[bpt+1], DiscMin, DiscMax);
		      
		      TString SystFunction = ThisFunction;
		      SystFunction += UncertaintyHisto->GetBinContent(bpt+1);
		      const TF1 SFFun("SFFun", SystFunction, PtBinEdge[bpt], PtBinEdge[bpt+1]);
		      BTagEntry e1(&SFFun, params);  
		      
		      calib.addEntry(e1);
		      
		    }
		    
		  }
		  
		}
	    
	  }	
	  
	}
	  
    }
  
  std::ofstream outFile(BTagger + ".csv");
  calib.makeCSV(outFile);
  outFile.close();
  
}

