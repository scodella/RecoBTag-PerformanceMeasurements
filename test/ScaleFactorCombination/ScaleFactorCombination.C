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

#include "ScaleFactorCombination_Weight.h"
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
/*
#include "Measurements/7TeVLegacy/PtRel.C"
#include "Measurements/7TeVLegacy/System8.C"
#include "Measurements/7TeVLegacy/IP3D.C"
#include "Measurements/7TeVLegacy/LifetimeTagger.C"
#include "Measurements/7TeVLegacy/LifetimeTaggerJPsi.C"
#include "Measurements/7TeVLegacy/LifetimeTaggerTTbarDilepton.C"
#include "Measurements/7TeVLegacy/KinematicTTbarDilepton.C"
#include "Measurements/7TeVLegacy/TagCountingTTbarDilepton.C"
//#include "Measurements/7TeVLegacy/KinematicTTbarLeptonJets.C"
//#include "Measurements/7TeVLegacy/TagCountingTTbarLeptonJets.C"
*/
 
// debug: print sym. matrix
void printMatrix (TMatrixDSym& mat, const char* tit = "matrix")
{
  std::cout << tit << std::endl;
  for ( int i=0; i<mat.GetNrows(); ++i ) {
    for ( int j=0; j<mat.GetNrows(); ++j ) 
      printf("%12.9f ",mat(i,j));
    printf("\n");
  }
}

// debug: print matrix
void printMatrix (TMatrixD& mat, const char* tit = "matrix")
{
  std::cout << tit << std::endl;
  for ( int i=0; i<mat.GetNrows(); ++i ) {
    for ( int j=0; j<mat.GetNrows(); ++j ) 
      printf("%12.9f ",mat(i,j));
    printf("\n");
  }
}

// clear matrix contents
void clearMatrix (TMatrixDSym& mat)
{
  for ( int i=0; i<mat.GetNrows(); ++i ) {
    for ( int j=0; j<mat.GetNrows(); ++j )  mat(i,j) = 0;
  }  
}

// square
inline double sqr (double arg) {return arg*arg;}

// combined variance for a list of uncertainties
double combVar (double e0, double e1=0, double e2=0, double e3=0, double e4=0, 
                double e5=0, double e6=0, double e7=0, double e8=0, double e9=0,
                double e10=0, double e11=0, double e12=0, double e13=0, double e14=0)
{
  double sum(0.);
  sum += sqr(e0);
  sum += sqr(e1);
  sum += sqr(e2);
  sum += sqr(e3);
  sum += sqr(e4);
  sum += sqr(e5);
  sum += sqr(e6);
  sum += sqr(e7);
  sum += sqr(e8);
  sum += sqr(e9);
  sum += sqr(e10);
  sum += sqr(e11);
  sum += sqr(e12);
  sum += sqr(e13);
  sum += sqr(e14);
  return sum;
}

// combined s.d. for a list of uncertainties
double combError (double e0, double e1=0, double e2=0, double e3=0, double e4=0,
                  double e5=0, double e6=0, double e7=0, double e8=0, double e9=0,
                  double e10=0, double e11=0, double e12=0, double e13=0, double e14=0)
{
  return sqrt(combVar(e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14));
}

// fill off-diagonal terms of a covariance matrix assuming 100% (anti-)correlation
void correlate (TMatrixDSym& matrix, float* correlations)
{
  size_t ind = 0;
  for ( int i=0; i<matrix.GetNrows(); ++i ) {
    for ( int j=i+1; j<matrix.GetNrows(); ++j ) {
      matrix(i,j) = matrix(j,i) = correlations[ind++]*sqrt(matrix(i,i)*matrix(j,j));
    }
  }
}

// remove unused rows & columns from matrix
TMatrixDSym compactify (TMatrixDSym& cov, TVectorD& vec)
{
  size_t n(0);
  std::vector<int> flags(nMeasurements,-1);
  for (int nm = 0; nm<nMeasurements; nm++)
    if (vec[nm]>0.) flags[nm] = n++;
 
  TMatrixDSym newCov(n);
  for ( size_t i=0; i<abs(nMeasurements); ++i ) {
    if ( flags[i]<0 )  continue;
    for ( size_t j=0; j<abs(nMeasurements); ++j ) {
      if ( flags[j]<0 )  continue;
      newCov(flags[i],flags[j]) = cov(i,j);
    }
  }
  return newCov;
}

// remove unused rows & columns from vector
TVectorD compactify (TVectorD& vec)
{
  size_t n(0);
  std::vector<int> flags(nMeasurements,-1);
  for (int nm = 0; nm<nMeasurements; nm++)
    if (vec[nm]>0.) flags[nm] = n++;

  TVectorD newVec(n);
  for ( size_t i=0; i<abs(nMeasurements); ++i ) {
    if ( flags[i]<0 )  continue;
    newVec(flags[i]) = vec[i];
  }
  return newVec;
}

//
// LSQ combination of several correlated variables
//
TVectorD BlueCoefficients(0);
std::pair<double, double> combine (const TVectorD& vec, const TMatrixDSym& cov, double& chi2)
{
  std::pair<double,double> result(0.,0.);
  chi2 = -1.;

  // weight matrix
  double det;
  TMatrixDSym wgt(cov);
  wgt.Invert(&det);
  // std::cout << "Det = " << det << std::endl;
  if ( det<0 )  return result;

 // for ( size_t i=0; i<cov.GetNrows(); ++i ) {
 //   printf("%8.3f ; ",vec(i));
 //   for ( size_t j=0; j<cov.GetNrows(); ++j )  printf("%9.6f ",cov(i,j));
 //   printf("\n");
 // }

  // auxiliary matrix (Nx1)
  TVectorD atg(wgt.GetNrows());
  for ( int i=0; i<wgt.GetNrows(); ++i ) {
    for ( int j=0; j<wgt.GetNrows(); ++j )  atg(i) += wgt(i,j);
  }

  double atgm(0);
  double atga(0);
  for ( int i=0; i<wgt.GetNrows(); ++i ) {
    atgm += atg(i)*vec(i);
    atga += atg(i);
  }
  result.first = atgm/atga;
  result.second = sqrt(1./atga);
  
  BlueCoefficients.ResizeTo(wgt.GetNrows());
  for ( int i=0; i<wgt.GetNrows(); ++i )
    BlueCoefficients(i) = atg(i)/atga;
  
  TVectorD dm(vec);
  for ( int i=0; i<wgt.GetNrows(); ++i ) dm(i) -= result.first;
  chi2 = wgt.Similarity(dm);

  double wgtMin(999999);
  for ( int i=0; i<wgt.GetNrows(); ++i ) {
    double wi = atg(i)/atga;
    if ( wi<wgtMin ) wgtMin = wi;
  }
//$$
  if ( wgtMin<0 ) { 
    std::cout << std::endl;
    std::cout << "**************************************" << std::endl;
    std::cout << "**** at least one negative weight ****" << std::endl;
    std::cout << "**************************************" << std::endl;
  }
//$$
  return result;
}

TMatrixDSym BuildErrorMatrix(int k, TString ErrorCategory = "") {

 TMatrixDSym totalCov(nMeasurements);
 clearMatrix(totalCov);
 TMatrixDSym currCov(nMeasurements);

 float SpecificUncertainty[nMeasurements];
 for (int nm = 0; nm<nMeasurements; nm++) 
   SpecificUncertainty[nm] = sqr(MeasuredScaleFactorUncertainty[nm][k][0]);

 for (int is = 1; is<nUncertainties; is++) {

   if (ErrorCategory=="" || ErrorCategory==UncertaintyName[is]) {
     
     bool UsedUncertainty = true;
     
     if (UncertaintyName[is]=="_statistic") {

       // start with stat errors
       // start with stat error (fully correlated)
       // totalCov *= 0.;
       // totalCov(0,0) = sqr(sf_stat_PT[k]);
       // totalCov(1,1) = sqr(sf_stat_S8[k]);
       // totalCov(2,2) = sqr(sf_stat_IP[k]);
       // totalCov(3,3) = sqr(sf_stat_JP[k]);
       // correlate(totalCov,corrStat);
       // assume PT/IP contained in SY, SY contained in JP
       // total stat variance is divided according to the fraction of the subsample
       // for each tagger
       // events used for PT/IP, SY and JP
       
       float PTscale = sqrt(1.);
       float S8scale = sqrt(1.); 
       float IPscale = sqrt(1.);
       float JPscale = sqrt(1.);
       
       clearMatrix(currCov);
       currCov(0,0) = sqr(PTscale*MeasuredScaleFactorUncertainty[0][k][is]);
       currCov(1,1) = frac_PTS8[k]*sqr(S8scale*MeasuredScaleFactorUncertainty[1][k][is]);
       currCov(2,2) = sqr(IPscale*MeasuredScaleFactorUncertainty[2][k][is]);
       currCov(3,3) = frac_PTS8[k]*frac_S8JP[k]*sqr(JPscale*MeasuredScaleFactorUncertainty[3][k][is]);
       correlate(currCov,corrStat);
       totalCov += currCov;
   
       // events used for SY and JP
       clearMatrix(currCov);
       currCov(1,1) = (1-frac_PTS8[k])*sqr(S8scale*MeasuredScaleFactorUncertainty[1][k][is]);
       currCov(3,3) = frac_S8JP[k]*(1-frac_PTS8[k])*sqr(JPscale*MeasuredScaleFactorUncertainty[3][k][is]);
       correlate(currCov,corrStat);
       totalCov += currCov;
   
       // remaining (JP only)
       clearMatrix(currCov);
       currCov(3,3) = (1-frac_S8JP[k])*sqr(JPscale*MeasuredScaleFactorUncertainty[3][k][is]);
       correlate(currCov,corrStat);
       totalCov += currCov;
   
       // PSI (assumed uncorrelated)
       totalCov(4,4) += sqr(MeasuredScaleFactorUncertainty[4][k][is]);

       // TT (assumed uncorrelated, but correlated between them)
       //totalCov(5,5) += sqr(sf_stat_TT[k]);
       clearMatrix(currCov);
       currCov(5,5) = sqr(MeasuredScaleFactorUncertainty[5][k][is]);
       currCov(6,6) = sqr(MeasuredScaleFactorUncertainty[6][k][is]);
       currCov(7,7) = sqr(MeasuredScaleFactorUncertainty[7][k][is]);
       correlate(currCov,corrStat);
       totalCov += currCov;

     } else if (UncertaintyName[is]=="_pileup" || UncertaintyName[is]=="_mupt" || UncertaintyName[is]=="_bfragmentation" || UncertaintyName[is]=="_jetaway" || UncertaintyName[is]=="_mudr" || UncertaintyName[is]=="_ltbias" || UncertaintyName[is]=="_cb" || UncertaintyName[is]=="_jes" || UncertaintyName[is]=="_ttbarmodelling" || UncertaintyName[is]=="_l2c" || UncertaintyName[is]=="_ptrel" || UncertaintyName[is]=="_s8closure" || UncertaintyName[is]=="_ipbias" || UncertaintyName[is]=="_ltaway" || UncertaintyName[is]=="_ltinclusive" || UncertaintyName[is]=="_ltothers" || UncertaintyName[is]=="_jer" || UncertaintyName[is]=="_psiother" || UncertaintyName[is]=="_tt" || UncertaintyName[is]=="_kt" || UncertaintyName[is]=="_tct") {

       // pileup (fully correlated)
       // mupt (fully correlated)
       // gluon (fully correlated)
       // b fragmentation (fully correlated)
       // away (fully correlated between PT and IP; uncorrelated for SY; none for JP)
       // delta R muon (fully correlated between PT/IP/SY; none for JP)
       // bias (fully correlated between JP and PSI, not for TT as only using JP fits)
       // Cb (fully correlated between JP, PSI and TT)
       // jes
       // tt model
       
       currCov *= 0.;
       for (int nm = 0; nm<nMeasurements; nm++)
	 if (MeasurementSystematic[nm][is]!=0)
	   currCov(nm, nm) = sqr(MeasuredScaleFactorUncertainty[nm][k][is]);
       correlate(currCov,SystematicCorrelation[is][k]);
       totalCov += currCov;
       
     } else {
       
       UsedUncertainty = false;
       
     }
     
     if (UsedUncertainty) {

       for (int nm = 0; nm<nMeasurements; nm++)
	 if (MeasurementSystematic[nm][is]!=0)
	   SpecificUncertainty[nm] -= sqr(MeasuredScaleFactorUncertainty[nm][k][is]);;
       
     }
     
   } else { cout << "    BuildErrorMatrix: wrong error category " << ErrorCategory << endl; }

 }

 for (int nm = 0; nm<nMeasurements; nm++)
   if (SpecificUncertainty[nm]>0.00001) {

     currCov *= 0.;
     currCov(nm, nm) = SpecificUncertainty[nm];
     totalCov += currCov;
     
   }
 
 // vector of measurements
 TVectorD sfs(nMeasurements);
 for (int nm = 0; nm<nMeasurements; nm++) {

   sfs[nm] = 0.;

   if (UseThisMeasurement[nm]) 
     sfs[nm] = MeasuredScaleFactor[nm][k];
   
 }
 
 // removal of unused taggers
 TMatrixDSym cov = compactify(totalCov, sfs);

 return cov;

}

// (result,error) for combination of bin "k" & tagger "tagger" for the methods indicated in the list
std::pair<double,double> combine (int k, double& CHI2)
{
 
 TMatrixDSym cov = BuildErrorMatrix(k); 

 // vector of measurements
 TVectorD sfs(nMeasurements);
 for (int nm = 0; nm<nMeasurements; nm++) {

   sfs[nm] = 0.;

   if (UseThisMeasurement[nm]) 
     sfs[nm] = MeasuredScaleFactor[nm][k];
   
 }
 TVectorD vec = compactify(sfs);

 double chi2;
 std::pair<double,double> result = combine(vec,cov,chi2);
 
 //$$
 std::cout << k << " " << result.first << " " << result.second << std::endl;
 double chi2norm(chi2);
 if ( cov.GetNrows()>1 )  chi2norm /= (cov.GetNrows()-1);
 std::cout << "NDF = " << cov.GetNrows()-1 << " chi2 = " << chi2 << " " << chi2norm << std::endl;
 if ( chi2norm > 1. ) result.second *= sqrt(chi2norm);
 if ( chi2norm > 1. ) {
   std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
   std::cout << "  !! chi2/norm too large !!  " << std::endl;
   std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
 }
 if ( chi2norm < 0.01 ) {
   std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
   std::cout << "  !! chi2/norm too small !!  " << std::endl;
   std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
 }
 //$$
 CHI2 = chi2norm;
 //$$
 
 return result;

}

void ScaleFactorsFit(TString ErrorCategory = "") {
  
  //int nTotalMeasurements = 0;
  //for (int bpt = 0; bpt<nBinsCampaign; bpt++) nTotalMeasurements += MeasurementsForBin[bpt];

  cout << "  Build MatU" << endl;

  TMatrixD MatU(nTotalMeasurements, nBinsCampaign);

  MatU *= 0.;

  int RowIndex = 0;
  //for (int bpt = 0; bpt<nBinsCampaign; bpt++) 
  //for (int im = 0; im<MeasurementsForBin[bpt]; im++) {
  for (int im = 0; im<nTotalMeasurements; im++) {

    int bpt = CampaignBinIndex[im];

    int jm = MeasurementBinIndex[im];
    if (MeasurementName[jm]!="TagCountTT") MatU(RowIndex, bpt) = 1.;
    else {
      
      for (int bpt2 = 0; bpt2<nBinsCampaign; bpt2++) 
	MatU(RowIndex, bpt2) = ttbar_pt[bpt2]/ttbar_tot;

    }
    
    RowIndex++;
    
  }
  
  cout << "  Build MatCov" << endl;

  TMatrixD MatCov(nTotalMeasurements, nTotalMeasurements);

  MatCov *= 0.;

  RowIndex = 0;
  //for (int bpt1 = 0; bpt1<nBinsCampaign; bpt1++) 
  //for (int im1 = 0; im1<MeasurementsForBin[bpt1]; im1++) {
  for (int im1 = 0; im1<nTotalMeasurements; im1++) {

      int ColumnIndex = 0;
      //for (int bpt2 = 0; bpt2<nBinsCampaign; bpt2++) 
      //for (int im2 = 0; im2<MeasurementsForBin[bpt2]; im2++) {
      for (int im2 = 0; im2<nTotalMeasurements; im2++) {

	  int jm1 = MeasurementBinIndex[im1];
	  int jm2 = MeasurementBinIndex[im2];

	  int bpt1 = CampaignBinIndex[im1];
	  int bpt2 = CampaignBinIndex[im2];

	  float SpecificUncertainty = 0.;
	  if (jm1==jm2) SpecificUncertainty = sqr(MeasuredScaleFactorUncertainty[jm1][bpt1][0]);
	  
	  for (int is = 1; is<nUncertainties; is++) {
	    
	    if (ErrorCategory=="" || ErrorCategory==UncertaintyName[is]) {
	      
	      bool UsedUncertainty = true;
      
	      if (bpt1==bpt2 
		  //|| UncertaintyName[is]=="_gluonsplitting" || UncertaintyName[is]=="_bfragmentation" 
		  //|| UncertaintyName[is]=="_mupt" 
		  //|| UncertaintyName[is]=="_jetaway" 
		  //|| UncertaintyName[is]=="_mudr" 
		  || UncertaintyName[is]=="_ptrel" 
		  //|| UncertaintyName[is]=="_ltothers"
		  ) { // Ignore bin to bin correlations for the moment (just like Run 1)
		
		if (UncertaintyName[is]=="_statistic") {
		  
		  float ThisCoefficient = fabs(MeasuredScaleFactorUncertainty[jm1][bpt1][is]*MeasuredScaleFactorUncertainty[jm2][bpt2][is]);
		  
		  if (jm1==jm2) 
		    MatCov(RowIndex, ColumnIndex) += ThisCoefficient;
		  else {
		    
		    if ((MeasurementName[jm1]=="PtRel" && MeasurementName[jm2]=="System8") ||
			(MeasurementName[jm2]=="PtRel" && MeasurementName[jm1]=="System8")) {
			  
		      MatCov(RowIndex, ColumnIndex) += sqrt(frac_PTS8[bpt1])*ThisCoefficient;

		    } else if ((MeasurementName[jm1]=="PtRel" && MeasurementName[jm2]=="IP3D") ||
			       (MeasurementName[jm2]=="PtRel" && MeasurementName[jm1]=="IP3D")) {
		      
		      MatCov(RowIndex, ColumnIndex) += ThisCoefficient;
		      
		    } else if ((MeasurementName[jm1]=="PtRel" && MeasurementName[jm2]=="LT") ||
			       (MeasurementName[jm2]=="PtRel" && MeasurementName[jm1]=="LT")) {
		      
		      MatCov(RowIndex, ColumnIndex) += sqrt(frac_PTS8[bpt1])*frac_S8JP[bpt1]*ThisCoefficient;

		    } else if ((MeasurementName[jm1]=="System8" && MeasurementName[jm2]=="IP3D") ||
			       (MeasurementName[jm2]=="System8" && MeasurementName[jm1]=="IP3D")) {
		      
		      MatCov(RowIndex, ColumnIndex) += sqrt(frac_PTS8[bpt1])*ThisCoefficient;

		    } else if ((MeasurementName[jm1]=="LT" && MeasurementName[jm2]=="IP3D") ||
			       (MeasurementName[jm2]=="LT" && MeasurementName[jm1]=="IP3D")) {
		      
		      MatCov(RowIndex, ColumnIndex) += sqrt(frac_PTS8[bpt1]*frac_S8JP[bpt1])*ThisCoefficient;

		    } else if ((MeasurementName[jm1]=="System8" && MeasurementName[jm2]=="LT") ||
			       (MeasurementName[jm2]=="System8" && MeasurementName[jm1]=="LT")) {
		      
		      MatCov(RowIndex, ColumnIndex) += sqrt(frac_S8JP[bpt1])*ThisCoefficient;

		    } 
		    
		  }
		  
		} else if (UncertaintyName[is]=="_pileup" || UncertaintyName[is]=="_mupt" || UncertaintyName[is]=="_gluonsplitting" || UncertaintyName[is]=="_bfragmentation" || UncertaintyName[is]=="_jetaway" || UncertaintyName[is]=="_mudr" || UncertaintyName[is]=="_ltbias" || UncertaintyName[is]=="_cb" || UncertaintyName[is]=="_jes" || UncertaintyName[is]=="_ttbarmodelling" || UncertaintyName[is]=="_l2c" || UncertaintyName[is]=="_ptrel" || UncertaintyName[is]=="_s8closure" || UncertaintyName[is]=="_ipbias" || UncertaintyName[is]=="_ltaway" || UncertaintyName[is]=="_ltinclusive" || UncertaintyName[is]=="_ltothers" || UncertaintyName[is]=="_jer" || UncertaintyName[is]=="_psiother" || UncertaintyName[is]=="_tt" || UncertaintyName[is]=="_kt" || UncertaintyName[is]=="_tct") {

		  // pileup (fully correlated)
		  // mupt (fully correlated)
		  // gluon (fully correlated)
		  // b fragmentation (fully correlated)
		  // away (fully correlated between PT and IP; uncorrelated for SY; none for JP)
		  // delta R muon (fully correlated between PT/IP/SY; none for JP)
		  // bias (fully correlated between JP and PSI, not for TT as only using JP fits)
		  // Cb (fully correlated between JP, PSI and TT)
		  // jes
		  // tt model

		  MatCov(RowIndex, ColumnIndex) += MeasuredScaleFactorUncertainty[jm1][bpt1][is]*MeasuredScaleFactorUncertainty[jm2][bpt2][is];

		} else {
		  
		  UsedUncertainty = false;
		  
		}
		
		if (UsedUncertainty && jm1==jm2 && bpt1==bpt2) {
		  
		  SpecificUncertainty -= sqr(MeasuredScaleFactorUncertainty[jm1][bpt1][is]);
		  
		}

	      }

	    } else { cout << "    BuildErrorMatrix: wrong error category " << ErrorCategory << endl; }
	      
	  }
	    
	  if (SpecificUncertainty>0.00001 && jm1==jm2 && bpt1==bpt2)
	    MatCov(RowIndex, ColumnIndex) += SpecificUncertainty;

	  ColumnIndex++;
	  
	}
      
      RowIndex++;
      
    }

  cout << "  Build VecSF" << endl;

  TMatrixD VecSF(nTotalMeasurements, 1);

  VecSF *= 0.;

  RowIndex = 0;
  for (int im = 0; im<nTotalMeasurements; im++) {
    //for (int bpt = 0; bpt<nBinsCampaign; bpt++) 
    //for (int im = 0; im<MeasurementsForBin[bpt]; im++) {

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
  cout << "MatCov " << Determinant << endl;
  if (Determinant==0.)
    cout << "MatCov not invertible" << endl;

  TMatrixD Aux1(nBinsCampaign, nTotalMeasurements);
  Aux1.Mult(MatUTr, MatCovInv);

  TMatrixD Aux2(nBinsCampaign, nBinsCampaign);
  Aux2.Mult(Aux1, MatU); 
  Aux2.Invert(&Determinant);

  if (Determinant==0.)
    cout << "Aux2 not invertible" << endl;
  
  TMatrixD MatCoefficients(nBinsCampaign, nTotalMeasurements);
  MatCoefficients.Mult(Aux2, Aux1);

  cout << "  Get SFs" << endl;

  TMatrixD MatSF(nBinsCampaign, 1);
  MatSF.Mult(MatCoefficients, VecSF);

  TMatrixD MatCoefficientsTr(nTotalMeasurements, nBinsCampaign);
  MatCoefficientsTr.Transpose(MatCoefficients);

  TMatrixD Aux3(nBinsCampaign, nTotalMeasurements);
  Aux3.Mult(MatCoefficients, MatCov); 

  cout << "  Get SF errors" << endl;

  TMatrixD MatSFErr(nBinsCampaign, nBinsCampaign);
  MatSFErr.Mult(Aux3, MatCoefficientsTr);

  cout << "  Fill fit results" << endl;

  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

    sf[bpt] = MatSF(bpt, 0);
    sf_eror[bpt] = sqrt(MatSFErr(bpt, bpt));
    sf_stat[bpt] = 0.;

  }

  cout << "  Compute Chi2" << endl;

  // Bin by bin
  
  if (Chi2str=="Bin") {
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

    TVectorD VecMeasurements(nTotalMeasurements);
    
    VecMeasurements *= 0.;

    nBinMeasurements[bpt] = 0;

    RowIndex = 0;
    for (int im = 0; im<nTotalMeasurements; im++) {
      
      int jm = MeasurementBinIndex[im];
      int bpt2 = CampaignBinIndex[im];

      if (MeasurementName[jm]!="TagCountTT") 
	if (bpt2==bpt) {
	  
	  VecMeasurements(RowIndex) = MeasuredScaleFactor[jm][bpt2] - sf[bpt];
	  nBinMeasurements[bpt]++;

	}
      
      RowIndex++;
      
    }

    Chi2Normal[bpt] = MatCovInv.Similarity(VecMeasurements);

    if (nBinMeasurements[bpt]>1)
      cout << "Bin Chi2/ndof = " << Chi2Normal[bpt]/(nBinMeasurements[bpt]-1) << " " << endl;

  }


  float sf_ttbar = 0., sf_ttbar_eror = 0.;
  float ttbar_spectrum_scale = ttbar_tot;
  for (int k = 0; k < nBinsCampaign; k++) {
    sf_ttbar += sf[k]*ttbar_pt[k]/ttbar_spectrum_scale; 
    sf_ttbar_eror += sf_eror[k]*ttbar_pt[k]/ttbar_spectrum_scale;
    //std::cout << "k: " << k << " " << xPt[k] << " " << sf_eror[k] << " " << ttbar_pt[k] << " " << fun_val[k] << " " << sf[k] << " " << Chi2Normal[k] << std::endl;
  }
  
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
    
    Chi2Normal[bpt] += sqr((MeasuredScaleFactor[7][0] - sf_ttbar))/(sqr(MeasuredScaleFactorUncertainty[7][0][1]/TTbarScale)+sqr(MeasuredScaleFactorUncertainty[7][0][19]));
    nBinMeasurements[bpt]++;
    
    if (nBinMeasurements[bpt]>1.) Chi2Normal[bpt] /= (nBinMeasurements[bpt]-1);
    if (Chi2Normal[bpt]>1.) sf_eror[bpt] *= sqrt(Chi2Normal[bpt]);
    
  }
  }
  // Overall
  if (Chi2str=="Overall") {

  float TotalChi2 = 0.;
  for (int im1 = 0; im1<nTotalMeasurements; im1++)
    for (int im2 = 0; im2<nTotalMeasurements; im2++) 
      for (int bpt1 = 0; bpt1<nBinsCampaign; bpt1++) 
	for (int bpt2 = 0; bpt2<nBinsCampaign; bpt2++) {
	  
	  int jm1 = MeasurementBinIndex[im1];
	  int jm2 = MeasurementBinIndex[im2];

	  if (MeasurementName[jm1]!="TagCountTT" && MeasurementName[jm2]!="TagCountTT") {
	    
	    if (bpt1==CampaignBinIndex[im1] && bpt2==CampaignBinIndex[im2]) 
	      TotalChi2 += MatU(im1, bpt1)*(MeasuredScaleFactor[jm1][bpt1] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, bpt2)*(MeasuredScaleFactor[jm2][bpt2] - sf[bpt2]);
	    
	  } else if (MeasurementName[jm1]=="TagCountTT" && MeasurementName[jm2]=="TagCountTT") {
	    TotalChi2 += MatU(im1, bpt1)*(MeasuredScaleFactor[jm1][0] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, bpt2)*(MeasuredScaleFactor[jm2][0] - sf[bpt2]);
	   
	  }

	}
  
  TotalChi2 /= (nTotalMeasurements - nBinsCampaign);
  cout << "TotalChi2 = " << TotalChi2 << endl;
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
    
    Chi2Normal[bpt] = TotalChi2;
    if (Chi2Normal[bpt]>1.) sf_eror[bpt] *= sqrt(Chi2Normal[bpt]);
    
  }

  float ttbar_spectrum_scale = 0.; float ErrCut = 0.15;
  for (int k = 0; k < nBinsCampaign; k++) 
    if (sf_eror[k]<ErrCut) ttbar_spectrum_scale += ttbar_pt[k];

  float sf_ttbar = 0., sf_ttbar_eror = 0., sf_ttbar_eror_up = 0.;
  for (int k = 0; k < nBinsCampaign; k++) {
    if (sf_eror[k]>=ErrCut) continue; 
    sf_ttbar += sf[k]*ttbar_pt[k]/ttbar_spectrum_scale; 
    sf_ttbar_eror += sf_eror[k]*ttbar_pt[k]/ttbar_spectrum_scale;
    sf_ttbar_eror_up += (sf[k]+sf_eror[k])*ttbar_pt[k]/ttbar_spectrum_scale;
  }
  sf_ttbar_eror_up -= sf_ttbar;
  //cout << "E " << sf_ttbar_eror << " " << sf_ttbar_eror_up << endl;

  float Scaling = fabs(sf_ttbar - MeasuredScaleFactor[7][0])/sqrt(sqr(sf_ttbar_eror)+sqr(MeasuredScaleFactorUncertainty[7][0][1]/TTbarScale)+sqr(MeasuredScaleFactorUncertainty[7][0][19]));
  if (Scaling<1.) Scaling = 1.;

  cout << "  Scaling = " << Scaling << endl;

  float sf_ttbar_eror2 = 0., sf_ttbar_eror3 = 0.;
  for (int k = 0; k < nBinsCampaign; k++) {
    if (sf_eror[k]>=ErrCut) continue; 
    sf_ttbar_eror2 += sf_eror[k]*Scaling*ttbar_pt[k]/ttbar_spectrum_scale;
    sf_ttbar_eror3 += sf_eror[k]*(1.+ttbar_pt[k]/ttbar_spectrum_scale)*ttbar_pt[k]/ttbar_spectrum_scale;
  }
  
  for (int k = 0; k < nBinsCampaign; k++) {
    Chi2Normal[k] = sqr(1.+ttbar_pt[k]/ttbar_spectrum_scale*((sf_ttbar_eror2-sf_ttbar_eror)/(sf_ttbar_eror3-sf_ttbar_eror)));
    sf_eror[k] *= sqrt(Chi2Normal[k]);
  }


  } // Overall

  // Bin overall
  /*
  for (int bpt1 = 0; bpt1<nBinsCampaign; bpt1++) {

    float TotalChi2 = 0.;

    for (int im1 = 0; im1<nTotalMeasurements; im1++)
      for (int im2 = 0; im2<nTotalMeasurements; im2++) 
	for (int bpt2 = 0; bpt2<nBinsCampaign; bpt2++) {
	  
	  int jm1 = MeasurementBinIndex[im1];
	  int jm2 = MeasurementBinIndex[im2];
	  
	  if (MeasurementName[jm1]!="TagCountTT" && MeasurementName[jm2]!="TagCountTT") {
	    
	    if (bpt1==CampaignBinIndex[im1] && bpt2==CampaignBinIndex[im2]) 
	      TotalChi2 += MatU(im1, bpt1)*(MeasuredScaleFactor[jm1][bpt1] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, bpt2)*(MeasuredScaleFactor[jm2][bpt2] - sf[bpt2]);
	    
	  } else if (MeasurementName[jm1]=="TagCountTT" && MeasurementName[jm2]=="TagCountTT") {
	    TotalChi2 += MatU(im1, bpt1)*(MeasuredScaleFactor[jm1][0] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, bpt2)*(MeasuredScaleFactor[jm2][0] - sf[bpt2]);
	   
	  }
	  
	}
    
    if (bpt1<4)
      TotalChi2 /= (4 - 1);
        
    Chi2Normal[bpt1] = TotalChi2;
    if (Chi2Normal[bpt1]>1.) sf_eror[bpt1] *= sqrt(Chi2Normal[bpt1]);
    
  }
  */

  // Average
  if (Chi2str=="Average") {
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

    TVectorD VecMeasurements(nTotalMeasurements);
    
    VecMeasurements *= 0.;

    nBinMeasurements[bpt] = 0;

    RowIndex = 0;
    for (int im = 0; im<nTotalMeasurements; im++) {
      
      int jm = MeasurementBinIndex[im];
      int bpt2 = CampaignBinIndex[im];

      if (MeasurementName[jm]!="TagCountTT") 
	if (bpt2==bpt) {
	  
	  VecMeasurements(RowIndex) = MeasuredScaleFactor[jm][bpt2] - sf[bpt];
	  nBinMeasurements[bpt]++;

	}
      
      RowIndex++;
      
    }

    Chi2Normal[bpt] = MatCovInv.Similarity(VecMeasurements);

  }

  float sf_ttbar = 0., sf_ttbar_eror = 0.;
  float ttbar_spectrum_scale = ttbar_tot;
  for (int k = 0; k < nBinsCampaign; k++) {
    sf_ttbar += sf[k]*ttbar_pt[k]/ttbar_spectrum_scale; 
    sf_ttbar_eror += sf_eror[k]*ttbar_pt[k]/ttbar_spectrum_scale;
  }
  
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
    
    if (TTbarScale!=11.)
      Chi2Normal[bpt] += sqr((MeasuredScaleFactor[7][0] - sf_ttbar))/(sqr(MeasuredScaleFactorUncertainty[7][0][1]/TTbarScale)+sqr(MeasuredScaleFactorUncertainty[7][0][19]));
    else 
      Chi2Normal[bpt] = sqr((MeasuredScaleFactor[7][0] - sf_ttbar))/(sqr(MeasuredScaleFactorUncertainty[7][0][1]/TTbarScale)+sqr(MeasuredScaleFactorUncertainty[7][0][19]));
    nBinMeasurements[bpt]++;
    
    if (TTbarScale!=11.)
      if (nBinMeasurements[bpt]>1.) Chi2Normal[bpt] /= (nBinMeasurements[bpt]-1);
    //if (Chi2Normal[bpt]>1.) sf_eror[bpt] *= sqrt(Chi2Normal[bpt]);
    
  }

  float sf_ttbar_eror2 = 0., sf_ttbar_eror3 = 0.;
  for (int k = 0; k < nBinsCampaign; k++) {
    float ThisScaling = 1.;
    if (Chi2Normal[k]>1.) ThisScaling = sqrt(Chi2Normal[k]);
    sf_ttbar_eror2 += sf_eror[k]*ThisScaling*ttbar_pt[k]/ttbar_spectrum_scale;
    sf_ttbar_eror3 += sf_eror[k]*(1.+ttbar_pt[k]/ttbar_spectrum_scale)*ttbar_pt[k]/ttbar_spectrum_scale;
  }

  for (int k = 0; k < nBinsCampaign; k++) {
    Chi2Normal[k] = sqr(1.+ttbar_pt[k]/ttbar_spectrum_scale*((sf_ttbar_eror2-sf_ttbar_eror)/(sf_ttbar_eror3-sf_ttbar_eror)));
    sf_eror[k] *= sqrt(Chi2Normal[k]);
  }
  }

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
    if (UseThisMeasurement[nm]) {

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
	    
	    //MeasurementBinIndex[bpt][MeasurementsForBin[bpt]] = nm;
	    //MeasurementsForBin[bpt]++;
	    MeasurementBinIndex[nTotalMeasurements] = nm;
	    CampaignBinIndex[nTotalMeasurements] = bpt;
	    nTotalMeasurements++;
	    
	  }
	  
	  LowBinEdge = xPt[bpt] - exPt[bpt];

	  LastMeasuredSF = MeasuredScaleFactor[nm][bpt];
	  
	}      }

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
	cout << MeasurementName[nm] << " " << xpt[nm][bpt] << " " << MeasuredScaleFactorError[nm][bpt] << endl;
      }

      if (MeasurementName[nm]=="TagCountTT") 
	for (int cpt = 0; cpt<nBinsCampaign; cpt++) 
	  MeasuredScaleFactorUncertainty[nm][cpt][1] *= TTbarScale;

      if (InflateStatistic) {

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

  for (int is = 2; is<nUncertainties; is++) {
 
    int ThisCounter = 0;

    for (int im = 0; im<nMeasurements; im++) 
      for (int jm = im+1; jm<nMeasurements; jm++) {

	for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

	  SystematicCorrelation[is][bpt][ThisCounter] = 0.;

	  /*
	    if (UncertaintyName[is]=="_statistic") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrStat[ThisCounter];
	    else if (UncertaintyName[is]=="_pileup") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrSystPU[ThisCounter];
	    else if (UncertaintyName[is]=="_mupt") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrSystMuPt[ThisCounter];
	    else if (UncertaintyName[is]=="_gluonsplitting") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrSystGluon[ThisCounter];
	    else if (UncertaintyName[is]=="_bfragmentation") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrSystBfrag[ThisCounter];
	    else if (UncertaintyName[is]=="_jetaway") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrSystAway[ThisCounter];
	    else if (UncertaintyName[is]=="_mudr") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrSystMudr[ThisCounter];
	    else if (UncertaintyName[is]=="_cb") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrSystLife[ThisCounter];
	    else if (UncertaintyName[is]=="_jes") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrSystJES[ThisCounter];
	    else if (UncertaintyName[is]=="_ttbarmodelling") 
	    SystematicCorrelation[is][bpt][ThisCounter] = corrSystTTModel[ThisCounter];
	    else 
	    SystematicCorrelation[is][bpt][ThisCounter] = 0.;
	  */

	  /*
	  float TestCorr = 0.;
	  if (UncertaintyName[is]=="_statistic") 
	    TestCorr = corrStat[ThisCounter];
	  else if (UncertaintyName[is]=="_pileup") 
	    TestCorr = corrSystPU[ThisCounter];
	  else if (UncertaintyName[is]=="_mupt") 
	    TestCorr = corrSystMuPt[ThisCounter];
	  else if (UncertaintyName[is]=="_gluonsplitting") 
	    TestCorr = corrSystGluon[ThisCounter];
	  else if (UncertaintyName[is]=="_bfragmentation") 
	    TestCorr = corrSystBfrag[ThisCounter];
	  else if (UncertaintyName[is]=="_jetaway") 
	    TestCorr = corrSystAway[ThisCounter];
	  else if (UncertaintyName[is]=="_mudr") 
	    TestCorr = corrSystMudr[ThisCounter];
	  else if (UncertaintyName[is]=="_cb") 
	    TestCorr = corrSystLife[ThisCounter];
	  else if (UncertaintyName[is]=="_jes") 
	    TestCorr = corrSystJES[ThisCounter];
	  else if (UncertaintyName[is]=="_ttbarmodelling") 
	    TestCorr = corrSystTTModel[ThisCounter];
	  else 
	    TestCorr = 0.;
	  */
	  
	  if (UseThisMeasurement[im] && UseThisMeasurement[jm] && MeasurementSystematic[im][is]!=0 && MeasurementSystematic[jm][is]!=0) {
	    SystematicCorrelation[is][bpt][ThisCounter] = MeasuredScaleFactorUncertainty[im][bpt][is]*MeasuredScaleFactorUncertainty[jm][bpt][is]/
	      fabs(MeasuredScaleFactorUncertainty[im][bpt][is]*MeasuredScaleFactorUncertainty[jm][bpt][is]);
	    //if (bpt==0) cout << UncertaintyName[is] << " " << " " << im << " " << jm << " " << MeasuredScaleFactorUncertainty[im][bpt][is] << " " << MeasuredScaleFactorUncertainty[jm][bpt][is] << " " << SystematicCorrelation[is][bpt][ThisCounter] <<  " " << TestCorr << endl;
	  }
	  
	}

	ThisCounter++;

      }
  
  }

}
  
void SetMeasurementsToUse() {

  UseThisMeasurement[0] = Ptrel;
  UseThisMeasurement[1] = System8;
  UseThisMeasurement[2] = Ip3d;
  UseThisMeasurement[3] = LifeTime;
  UseThisMeasurement[4] = JPsi;
  UseThisMeasurement[5] = TTbar;
  UseThisMeasurement[6] = KinTT;
  UseThisMeasurement[7] = TagCntTT;

}

TStyle* PlotStyle()
{
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

  SetMeasurementsToUse();

  TString title = BTagger;
  if (OP=="Loose")  title += "L"; 
  if (OP=="Medium") title += "M"; 
  if (OP=="Tight")  title += "T"; 
  
  int stati=0;
  bool fit= 0;
  bool logx=1;

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  TCanvas *c1 = new TCanvas("c1", "plots",200,0,700,700);
  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2);
  
  TPad* pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.98,0.98,21);
  TPad* pad2 = new TPad("pad2","This is pad2",0.02,0.03,0.98,0.49,21);

  /* Run1 Setting  
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetFrameFillColor(10);
  pad1->Draw();
  pad1->SetLogx(logx);
  // pad1->SetGrid();
  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0.15);
  pad1->SetRightMargin(0.04);
  pad1->SetLeftMargin(0.16);
  pad1->SetTickx();
  pad1->SetTicky();
  
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetFrameFillColor(10);
  pad2->Draw();
  pad2->SetLogx(logx);
  // pad2->SetGrid();
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.15);
  pad2->SetRightMargin(0.04);
  pad2->SetLeftMargin(0.16);
  pad2->SetTickx();
  pad2->SetTicky();

  gStyle->SetOptDate(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleFont(62);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleW(0.4);
  gStyle->SetTitleH(0.09);
  // gStyle->SetTitleX(0); // Set the position of the title box
  // gStyle->SetTitleY(0.985); // Set the position of the title box
  // gStyle->SetTitleStyle(Style_t style = 1001);
  // gStyle->SetTitleBorderSize(2);
  gStyle->SetOptStat(stati);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  // gStyle->SetPadGridX(true); gStyle->SetPadGridY(true);
  gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);
  */ // End Run1 Setting

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

  /* Kirill's Style
  static TStyle* plotStyle = PlotStyle(); 
  gROOT->SetStyle("PLOT");   
  gROOT->ForceStyle(); 

  pad1->Draw();
  pad2->Draw();
  */// End Kirill's Style 

  if (fit) {
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.1);
    gStyle->SetOptFit(111);
  } else {
    gStyle->SetStatW(0.3);
    gStyle->SetStatH(0.2);
    gStyle->SetOptFit(0);
  }  

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  float exstat[20];
  for (int k = 0; k <= nBinsCampaign; k++) exstat[k] = 0.0001;
  
  // This is for avereging the results on ttbar spectrum (PAS BTV 11-00-04):
  //for (int k=0; k<=nBinsCampaign; k++) {
  //ttbar_pt[k] = ttbar_pt[k] / ttbar_tot;
  //}
  
  if (!StatisticalCorrelation) 
    for (int k = 0; k <= 14; k++) {
      frac_PTJP[k] = 1.; 
      frac_S8JP[k] = 1.; 
      frac_PTS8[k] = 1.;
    }

  //for (int bpt = 0; bpt<nBinsCampaign; bpt++) 
  //MeasurementsForBin[bpt] = 0;

  ReadMeasurements(BTagger, OP);

  // Define TGraph 

  pad1->cd();
  
  float PlotMaxPtCampaign = MaxPtCampaign;
  if (OP=="Tight") PlotMaxPtCampaign = 200.;
  if (OP=="Medium") PlotMaxPtCampaign = 300.;
  float YAxisWidth = 0.3;
  if (OP=="Tight") YAxisWidth = 0.4;
  TH2F* histo = new TH2F("histo","",58,MinPtCampaign,PlotMaxPtCampaign,100,1.-YAxisWidth,1.+YAxisWidth);
  
  histo->Draw(""); 
  histo->SetLabelSize(0.05, "XYZ");
  histo->SetTitleSize(0.06, "XYZ"); 
  histo->SetLabelFont(42, "XYZ"); 
  histo->SetTitleFont(42, "XYZ");
  histo->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  histo->GetYaxis()->SetTitle("Data/Sim. b-tag SF");
  histo->SetTitleOffset(1.1,"X");
  histo->SetTitleOffset(0.8,"Y");
  histo->SetTickLength(0.06,"X");
  histo->SetNdivisions(509, "XYZ");
  histo->GetXaxis()->SetMoreLogLabels();
  histo->GetXaxis()->SetNoExponent();
  
  MeasuredScaleFactor[0][8] = 0.;
  MeasuredScaleFactor[1][8] = 0.;
  if (BTagger=="TCHP") {
    MeasuredScaleFactor[1][6] = 0.;
    MeasuredScaleFactor[1][7] = 0.;
    MeasuredScaleFactorValue[1][3] = 0.;
  }

  TGraphErrors *ScaleFactorStatistic[nMeasurements];
  TGraphErrors *ScaleFactorError[nMeasurements];

  //for (int bpt = 0; bpt<20; bpt++) 
  //cout << bpt << " " << xpt[3][bpt] << " " << MeasuredScaleFactorValue[3][bpt] << endl;

  for (int nm = 0; nm<nMeasurements; nm++) {

    ScaleFactorStatistic[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], exstat, MeasuredScaleFactorStatistic[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], GraphWidth[nm]);

    ScaleFactorError[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], expt[nm], MeasuredScaleFactorError[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], 2);

    if (UseThisMeasurement[nm] && (MeasurementName[nm]!="TagCountTT" || TTbarScale==1.)) { 

      ScaleFactorStatistic[nm]->Draw("P"); 
      ScaleFactorError[nm]->Draw("P"); 

    }
    
  }

  TString LegTitle = title;

  //TLegend*      leg = new TLegend(0.28,0.65,0.50,0.93);
  //if ( !TTbar ) leg = new TLegend(0.28,0.69,0.50,0.93);
  float LegXOffset = 0.;
  if (OP=="Medium") LegXOffset = 0.06;
  TLegend *leg = new TLegend(0.48+LegXOffset,0.66,0.70+LegXOffset,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(62);
  leg->SetTextSize(0.05);   
  leg->SetHeader(LegTitle);
  for (int nm = 0; nm<nMeasurements; nm++) 
    if (UseThisMeasurement[nm] && (MeasurementName[nm]!="TagCountTT" || TTbarScale==1.)) {
      TString LegText = MeasurementName[nm];
      leg->AddEntry(ScaleFactorError[nm], LegText, "PL");
    }
  leg->Draw();

  /* // Run1 Style
  TLatex* tex = new TLatex(0.17,1,"CMS Preliminary, " + CampaignLuminosity + "  at #sqrt{s} = " + CenterOfMassEnergy);
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetLineWidth(2);
  tex->Draw();
  */ // End Run1 Style

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

  TLatex *text1 = new TLatex(0.98,0.95125, CampaignLuminosity + " (13 TeV, 50 ns)"); 
  text1->SetNDC();                                              
  text1->SetTextAlign(31);                          
  text1->SetTextFont(42);    
  text1->SetTextSize(0.04875);   
  text1->SetLineWidth(2);    
  text1->Draw(); 
  // End Run2015B Style

  //include the official CMS label
  if (OfficialCMS) {

    //TPaveText* pt = new TPaveText(0.55,0.21,0.90,0.31,"brNDC");   
    TPaveText* pt = new TPaveText(0.55,0.85,0.90,0.95,"brNDC");   
    pt->SetBorderSize(0);   
    pt->SetFillStyle(0);  
    pt->SetTextAlign(13);   
    pt->SetTextFont(22);   
    pt->SetTextSize(0.06);   
    pt->AddText("CMS prelim. at " + CenterOfMassEnergy + ", " + CampaignLuminosity);   
    pt->Draw(""); 

  }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  std::pair<double,double> combResult;

  //AlignPtBins();
  
  double CHI2;

  cout << "NOW STARTING THE FITS " << endl;

  ScaleFactorsFit();

  /*  
  for (int k = 0; k<nBinsCampaign; k++) {

    combResult = combine(k, CHI2);
 
    sf[k] = combResult.first;
    sf_eror[k] = combResult.second;
    sf_stat[k] = 0.;
    Chi2Normal[k] = CHI2;

    std::cout << "   FIT RESULTS FOR BIN " << k << " " << xPt[k] << " " << sf[k] <<  " +_ " << sf_eror[k] <<  " CHI2 " << CHI2 << std::endl;
    
    // Breakdown of systematics
    if (SystematicBreakdown) {
      
      cout << "     Breakdown of the systematics" << endl;

      float TotalUncertainty = 0.;

      for (int uu = 0; uu<nUncertainties; uu++) {
	
	TMatrixDSym ErrorMatrix = BuildErrorMatrix(k, UncertaintyName[uu]); 
	
	double Variance = 0.;
	for (int ii = 0; ii<BlueCoefficients.GetNoElements(); ii++)
	  for (int jj = 0; jj<BlueCoefficients.GetNoElements(); jj++)
	    Variance += BlueCoefficients(ii)*ErrorMatrix(ii,jj)*BlueCoefficients(jj);

	if (CHI2>1.) Variance *= CHI2;

	sf_uncerbreak[k][uu] = sqrt(Variance);

	cout << "       " << UncertaintyName[uu] << " " << sqrt(Variance) << endl;
	
	if (uu>0) TotalUncertainty += Variance;

      }
     
      cout << "       Total uncertainty " << sqrt(TotalUncertainty) << endl;
      
    }
    
  }
  */
  cout << "FITS DONE, NOW MAKING PLOTS" << endl;
  
  sf_eror[4] = MeasuredScaleFactorUncertainty[3][4][0];
  TGraphErrors* sf0 = new TGraphErrors(nBinsCampaign,xPt,sf,exPt,sf_eror);
  sf0->SetFillStyle(3005);
  sf0->SetFillColor(kGray+3);
  sf0->Draw("e2"); // to plot fit band
  
  //for (int nm = 0; nm<nMeasurements; nm++)
  //if (UseThisMeasurement[nm]) 
  //ScaleFactorError[nm]->Draw("P"); 
  // Keep the same sequence as in the old macro
  if (UseThisMeasurement[4]) ScaleFactorError[4]->Draw("P"); 
  if (UseThisMeasurement[1]) ScaleFactorError[1]->Draw("P");
  if (UseThisMeasurement[2]) ScaleFactorError[2]->Draw("P");
  if (UseThisMeasurement[3]) ScaleFactorError[3]->Draw("P");
  if (UseThisMeasurement[5]) ScaleFactorError[5]->Draw("P");
  if (UseThisMeasurement[0]) ScaleFactorError[0]->Draw("P");

  // #############################################################################
 
  float SFb_Comp[nBinsCampaign], SFb_Comp_error[nBinsCampaign];
 
  for (int k = 0; k<nBinsCampaign; k++) {

    double xx = xPt[k];

    if ( title == "CSVL" ) {
      SFb_Comp[k]      = funSFb_Comp_CSVL(xx);
      SFb_Comp_error[k] = SFb_Comp_error_CSVL[k];
    }
    if ( title == "CSVM" ) { 
      SFb_Comp[k]      = funSFb_Comp_CSVM(xx);
      SFb_Comp_error[k] = SFb_Comp_error_CSVM[k];
    }
    if ( title == "CSVT" ) { 
      SFb_Comp[k]      = funSFb_Comp_CSVT(xx);
      SFb_Comp_error[k] = SFb_Comp_error_CSVT[k];
    }
    if ( title == "TCHPT" ) { 
      SFb_Comp[k]      = funSFb_Comp_TCHPT(xx);
      SFb_Comp_error[k] = SFb_Comp_error_TCHPT[k];
    }
    
  }
  
  TGraphErrors* sf_Comp = new TGraphErrors(nBinsCampaign,xPt,SFb_Comp,exPt,SFb_Comp_error);

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  pad2->cd();
  cout << "Limit " << MinPtCampaign << " " << MaxPtCampaign << endl;
  TF1* fun1;
  if (CampaignName=="7TeVLegacy") {
    if ( title == "CSVL") fun1 = new TF1("fun1","[0]*(1.+[1]*x)/(1.+[2]*x)",MinPtCampaign,MaxPtCampaign);
    else if ( title == "CSVM" ) {
      fun1 = new TF1("fun1","[0]*(1.+[1]*x)/(1.+[2]*x)",MinPtCampaign,MaxPtCampaign);
      fun1->SetParameters(0.899, 0.0565, 0.0606);
    } else if ( title == "CSVT" ) {
      fun1 = new TF1("fun1","[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))",MinPtCampaign,MaxPtCampaign);
      fun1->SetParameters(0.9, 0.07/log(670)/log(670), 2./log(670));
    } else if ( title == "TCHPT" ) {
      fun1 = new TF1("fun1","[0]*(1.+[1]*x)/(1.+[2]*x)",MinPtCampaign,MaxPtCampaign);
    } 
  } else if (CampaignName=="Run2015B") {
    /*if ( title == "CSVv2L") fun1 = new TF1("fun1","[0]*(1.+[1]*x)/(1.+[2]*x)",MinPtCampaign,MaxPtCampaign);
    else if ( title == "CSVv2M" ) {
      fun1 = new TF1("fun1","[0]*(1.+[1]*x)/(1.+[2]*x)",MinPtCampaign,MaxPtCampaign);
      fun1->SetParameters(0.899, 0.0565, 0.0606);
    } else if ( title == "CSVv2T" ) {
      fun1 = new TF1("fun1","[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))",MinPtCampaign,MaxPtCampaign);
      fun1->SetParameters(0.9, 0.07/log(670)/log(670), 2./log(670));
    }*/
    //fun1 = new TF1("fun1","[0]");
    fun1 = new TF1("fun1","[0]*(1.+[1]*x)/(1.+[2]*x)",MinPtCampaign,MaxPtCampaign);
    fun1->SetParameters(0.853, 0.0527, 0.453);
    if (OP!="Loose") {
      fun1->FixParameter(1, 0.);
      fun1->FixParameter(2, 0.);
    }
  } else if (CampaignName=="Winter13") {
    if ( title == "CSVL") fun1 = new TF1("fun1","[0]*(1.+[1]*x)/(1.+[2]*x)",MinPtCampaign,MaxPtCampaign);
    else if ( title == "CSVM" || title == "CSVT") {
      fun1 = new TF1("fun1","[0]+[1]*x+[2]*x*x",MinPtCampaign,MaxPtCampaign);
      fun1->SetParameters(0.938887,0.00017124,-2.76366e-07);
    }
  }
  
  histo = new TH2F("histo","",58,MinPtCampaign,PlotMaxPtCampaign,100,1.-YAxisWidth,1.+YAxisWidth);

  histo->Draw(""); 
  histo->SetLabelSize(0.05, "XYZ");
  histo->SetTitleSize(0.06, "XYZ"); 
  histo->SetLabelFont(42, "XYZ"); 
  histo->SetTitleFont(42, "XYZ");
  histo->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  histo->GetYaxis()->SetTitle("Data/Sim. b-tag SF");
  histo->SetTitleOffset(1.1,"X");
  histo->SetTitleOffset(0.8,"Y");
  histo->SetTickLength(0.06,"X");
  histo->SetNdivisions(509, "XYZ");
  histo->GetXaxis()->SetMoreLogLabels();
  histo->GetXaxis()->SetNoExponent();
 
  sf_Comp->SetFillColor(kYellow-9);
  if (PrintComparison) sf_Comp->Draw("e2");
 
//   TGraphErrors* sf0stat = new TGraphErrors(15,x,sf,ex,sf_stat);
//   sf0stat->Draw("P"); 
//   sf0stat->SetMarkerColor(kBlack);
//   sf0stat->SetLineColor(kBlack);
//   sf0stat->SetLineStyle(1);
//   sf0stat->SetMarkerStyle(20);
//   sf0stat->SetMarkerSize(1.3);
//   sf0stat->SetLineWidth(2);
//   TGraphErrors* sf0 = new TGraphErrors(15,x,sf,ex,sf_eror);
  sf0->Draw("e2"); 
//   sf0->SetMarkerColor(kBlack);
//   sf0->SetLineColor(kBlack);
//   sf0->SetLineStyle(1);
//   sf0->SetMarkerStyle(20);
//   sf0->SetMarkerSize(1.3);
//   sf0->SetLineWidth(2);
  sf0->SetFillStyle(3005);
  sf0->SetFillColor(kBlack);

  fun1->SetLineColor(kBlack);
  fun1->SetLineWidth(2);
  fun1->SetLineStyle(1);
  float MaxPtFit = MaxPtCampaign;
  if (OP=="Loose") MaxPtFit = 400.;
  //if (OP=="Medium") MaxPtFit = 400.;
  //if (OP=="Tight") MaxPtFit = 200.;
  sf0->Fit("fun1","0","",MinPtCampaign,MaxPtFit);//MaxPtCampaign);
  sf0->Fit("fun1","rve0","",MinPtCampaign,MaxPtFit);//MaxPtCampaign);
//$$  sf0->Fit("fun1","rvee0"); 
  fun1->Draw("same"); 

  std::cout << std::endl;
  for (int k=0; k<nBinsCampaign; k++) {
    fun_val[k] = fun1->Integral(xPt[k]-1.,xPt[k]+1.) / 2.; //keep
    fun_err[k] = fun1->IntegralError(xPt[k]-exPt[k],xPt[k]+exPt[k]) / exPt[k];
    //std::cout << " Function " << k << " " << fun_val[k] << " " << fun_err[k] << std::endl;
    fun_sys[k] = sf_eror[k] * fun_val[k] / sf[k];
  }
  std::cout << std::endl;
  TGraphErrors* fun0 = new TGraphErrors(16,xPt,fun_val,exPt,fun_sys);
//   fun0->Draw("Pe2"); 
  fun0->Draw("P"); 
  fun0->SetMarkerColor(kRed);
  fun0->SetLineColor(kRed);
  fun0->SetLineStyle(1);
  fun0->SetMarkerStyle(24);
  fun0->SetMarkerSize(0.001);
  fun0->SetLineWidth(2);
//   fun0->SetFillStyle(3004);
//   fun0->SetFillColor(kRed);
  fun1->SetMarkerColor(kBlack);
  fun1->Draw("same"); 
  
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

  text1 = new TLatex(0.98,0.95125, CampaignLuminosity + " (13 TeV, 50 ns)"); 
  text1->SetNDC();                                              
  text1->SetTextAlign(31);                          
  text1->SetTextFont(42);    
  text1->SetTextSize(0.04875);   
  text1->SetLineWidth(2);    
  text1->Draw();  
  // End Run2015B Style

   //include the official CMS label
   if (OfficialCMS) {

     //TPaveText* pt = new TPaveText(0.55,0.21,0.90,0.31,"brNDC");   
     TPaveText* pt = new TPaveText(0.55,0.85,0.90,0.95,"brNDC");   
     pt->SetBorderSize(0);   
     pt->SetFillStyle(0);  
     pt->SetTextAlign(13);   
     pt->SetTextFont(22);   
     pt->SetTextSize(0.06);   
     pt->AddText("CMS prelim. at " + CampaignLuminosity + ", " + CenterOfMassEnergy);   
     pt->Draw(""); 

   }

   // Print results for payloads
   cout << "Results for payloads" << endl;
   TSfun1 = fun1->GetExpFormula("p");
   std::cout << " Tagger: " << title << " within " << MinPtCampaign 
	     << " < pt < " << MaxPtCampaign << " GeV, abs(eta) < 2.4, x = pt" << std::endl;
   std::cout << "  SFb = " << TSfun1 << ";"<< std::endl;
   std::cout << "  SFb_error[" << nBinsCampaign << "] = {" << std::endl;
   for (int k = 0; k < nBinsCampaign; k++) {
     //fun_sys[k] *= sqrt(Chi2Normal[k]);
     //sf_eror[k] *= sqrt(Chi2Normal[k]);
     if (k < 14 ) std::cout << "   " << fun_sys[k] << " " << std::endl;
     else         std::cout << "   " << fun_sys[k] << " };" << std::endl;
   }
   std::cout << std::endl;

   // Comparison with average ttbar-based measurements
   float sf_ttbar = 0., sf_ttbar_eror = 0.;
   float ttbar_spectrum_scale = ttbar_tot; 
   for (int k = 0; k < nBinsCampaign; k++) {
     sf_ttbar += fun_val[k]*ttbar_pt[k]/ttbar_spectrum_scale; 
     sf_ttbar_eror += sf_eror[k]*ttbar_pt[k]/ttbar_spectrum_scale * fun_val[k]/sf[k];
     std::cout << "k: " << k << " " << xPt[k] << " " << sf_eror[k] << " " << ttbar_pt[k] << " " << fun_val[k] << " " << sf[k] << " " << Chi2Normal[k] << std::endl;
   }
   std::cout << title << std::endl;
   std::cout << " ttbar-like SF = " << sf_ttbar << " +_ " << sf_ttbar_eror << std::endl;

   // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   
   c1->Update();

   TString PlotName = "./Plots/SFb_" + CampaignName + "_" + LegTitle + ".png";
   if (!Ptrel) PlotName.ReplaceAll(".png", "_noPT.png");
   if (!System8) PlotName.ReplaceAll(".png", "_noS8.png");
   if (!Ip3d) PlotName.ReplaceAll(".png", "_noIP.png");
   if (!LifeTime) PlotName.ReplaceAll(".png", "_noLT.png");
   if (!JPsi && CampaignName=="Winter13") PlotName.ReplaceAll(".png", "_noPSI.png");
   if (!TTbar) PlotName.ReplaceAll(".png", "_noTT.png");
   if (!KinTT && CampaignName!="Winter13") PlotName.ReplaceAll(".png", "_noKT.png");
   if (!TagCntTT && CampaignName!="Winter13") PlotName.ReplaceAll(".png", "_noTCT.png");
   if (TagCntTT && TTbarScale==100.) PlotName.ReplaceAll(".png", "_noTCT.png");
   PlotName.ReplaceAll(".png", "_Weight3.png");
   c1->Print(PlotName);
   PlotName.ReplaceAll(".png", ".pdf");
   c1->Print(PlotName);
   PlotName.ReplaceAll(".pdf", ".C");
   c1->Print(PlotName);

   cout << "Breakdown of the systematics" << endl;
   if (SystematicBreakdown) {
      
     for (int uu = 0; uu<nUncertainties; uu++) {

       TString VectorName = "  sf_uncer_" + UncertaintyName[uu] + "[";
       VectorName += nBinsCampaign; VectorName += "] = { "; 

       cout << endl << VectorName;  

       TString SpaceName = "";
       for (int ss = 0; ss<VectorName.Sizeof()-1; ss++) SpaceName += " ";

       for (int k = 0; k<nBinsCampaign; k++)
	 if (k<nBinsCampaign-1) cout << sf_uncerbreak[k][uu]*fun_val[k]/sf[k] << "," << endl << SpaceName;
	 else cout << sf_uncerbreak[k][uu]*fun_val[k]/sf[k] << " };" << endl;
       
     }

   }
  
}

void StoreScaleFactorCombination(string BTagger, TString OP = "All", string MeasType = "comb", bool StoreSystematicBreakdown = false) {

  BTagCalibration calib(BTagger);

  const int nOperatingPoints = 3;
  string OperatingPoint[nOperatingPoints] = {"Loose", "Medium", "Tight"};
  
  const int nMeasurementTypes = 2;
  string MeasName[nMeasurementTypes] = {"mujets", "comb"};

  SystematicBreakdown = StoreSystematicBreakdown;

  float PtBinEdge[nBinsCampaign+1]; PtBinEdge[0] = xPt[0] - exPt[0];
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) PtBinEdge[bpt+1] = xPt[bpt] + exPt[bpt];
  TH1F *UncertaintyHisto = new TH1F ("UncertaintyHisto", "", nBinsCampaign, PtBinEdge);  

  for (int mt = 0; mt<nMeasurementTypes; mt++) 
    if(MeasType=="All" || MeasType==MeasName[mt]) {

      if (MeasName[mt]=="mujets") {
	System8 = Ptrel = LifeTime = true; 
	JPsi = TTbar = KinTT = false;
	TTbarScale = 100.;
	if (CampaignName=="Winter13") JPsi = true;
      } else if (MeasName[mt]=="comb") {
	System8 = Ptrel = LifeTime = TagCntTT = true; 
	JPsi = KinTT = TTbar = false; 
	if (CampaignName=="Winter13") JPsi = true;
	TTbarScale = 1.;
      } 

      for (int op = 0; op<nOperatingPoints; op++) 
	if (OP=="All" || OP.Contains(OperatingPoint[op])) {
		  
	  BTagEntry::OperatingPoint wp;
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
		      UncertaintyHisto->SetBinContent(bpt+1, SysFactor*fun_sys[bpt]);

		  } else { 

		    SysFlag += "_"; SysFlag += UncertaintyName[uu-1]; 
		    for (int bpt = 0; bpt<nBinsCampaign; bpt++) 
		      UncertaintyHisto->SetBinContent(bpt+1, SysFactor*sf_uncerbreak[bpt][uu-1]*fun_val[bpt]/sf[bpt]);
		    
		  }

		}

		/* This is not working yet ...		  
		BTagEntry::Parameters params(wp, MeasName[mt], SysFlag, jfl, -2.4, 2.4, MinPtCampaign, MaxPtCampaign, 0, 1);
		
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

		  BTagEntry::Parameters params(wp, MeasName[mt], SysFlag, jfl, -2.4, 2.4, MinPtCampaign, MaxPtCampaign, 0, 1);

		  const TF1 SFFun("SFFun", ThisFunction, MinPtCampaign, MaxPtCampaign);
		  BTagEntry e1(&SFFun, params);  
		  
		  calib.addEntry(e1);
		  
		} else {

		  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

		    BTagEntry::Parameters params(wp, MeasName[mt], SysFlag, jfl, -2.4, 2.4, PtBinEdge[bpt], PtBinEdge[bpt+1], 0, 1);
		    
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

