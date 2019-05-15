const int nPtRelPtBins = 10;
TString PtRelPtBin[nPtRelPtBins] = {"Pt2030", "Pt3050", "Pt5070", "Pt70100", "Pt100140", "Pt140200", "Pt200300", "Pt300600", "Pt6001000", "Pt1000"};
double PtRelPtEdge[nPtRelPtBins+1] = {20., 30., 50., 70., 100., 140., 200., 300., 600., 1000., 1400.};

const int nPtRelEtaBins = 3;
TString PtRelEtaBin[nPtRelEtaBins] = {"anyEta", "Eta12", "Eta24"};
double PtRelEtaEdge[nPtRelEtaBins] = {0., 1.2, 2.4};

const int nTaggers = 6;
TString TaggerName[nTaggers] = {"DeepCSVL", "DeepCSVM", "DeepCSVT", "CSVv2L", "CSVv2M", "CSVv2T", };

const int nSystematics = 15; 
TString SystematicName[nSystematics] = {"_Central", "_GluonSplittingBDown", "_GluonSplittingBUp", "_MuPt6", "_MuPt8", "_MuDRMinus", "_MuDRPlus", "_JBPL", "_JBPT", "_BTemplatesMinus", "_BTemplatesPlus", "_PileUpDown", "_PileUpUp", "_JEUDown", "_JEUUp"};

//const int nSystematics = 1; 
//TString SystematicName[nSystematics] = {"_TCHPM"};

void System8Input(TString TemplateName = "_PSRun2017EFMoriond18_KinPtBinsCentral_LowPtAway_Run2016Production") {
  
  TString TemplateFileName = "Templates/System8_Templates" + TemplateName + ".root";
  TFile *InputFile = TFile::Open(TemplateFileName);

  TString DataName[2] = {"Data", "QCD"};

  for (int dt = 0; dt<2; dt++) 
    for (int is = 0; is<nSystematics; is++) 
      for (int tg = 0; tg<nTaggers; tg++) 
	for (int nb = 0; nb<1/*nPtRelEtaBins*/; nb++) {
	  
	  TString OutputFileName = "/afs/cern.ch/user/s/scodella/work/public/BTagging/System8/Inputs/Run2017Moriond18/Final/Run2017EF/S8_" + DataName[dt] + SystematicName[is] + "_" + TaggerName[tg] + "_" + PtRelEtaBin[nb] + ".root";
	  TFile *OutputFile = new TFile(OutputFileName, "recreate");
	  
	  TH2D *n_pT       = new TH2D("n_pT",       "n_pT",       nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  TH2D *ntag_pT    = new TH2D("ntag_pT",    "ntag_pT",    nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  TH2D *p_pT       = new TH2D("p_pT",       "p_pT",       nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  TH2D *ptag_pT    = new TH2D("ptag_pT",    "ptag_pT",    nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  
	  TH2D *n_pT_b     = new TH2D("n_pT_b",     "n_pT_b",     nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  TH2D *ntag_pT_b  = new TH2D("ntag_pT_b",  "ntag_pT_b",  nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  TH2D *p_pT_b     = new TH2D("p_pT_b",     "p_pT_b",     nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  TH2D *ptag_pT_b  = new TH2D("ptag_pT_b",  "ptag_pT_b",  nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  
	  TH2D *n_pT_cl    = new TH2D("n_pT_cl",    "n_pT_cl",    nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  TH2D *ntag_pT_cl = new TH2D("ntag_pT_cl", "ntag_pT_cl", nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  TH2D *p_pT_cl    = new TH2D("p_pT_cl",    "p_pT_cl",    nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  TH2D *ptag_pT_cl = new TH2D("ptag_pT_cl", "ptag_pT_cl", nPtRelPtBins, PtRelPtEdge, 70, 0., 7.);
	  
	  if (dt==0) {

	    for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {
	      
	      TString ThisHistoLabel = PtRelPtBin[ptb] + "_" + PtRelEtaBin[nb] + SystematicName[is];

	      TH1D *PtRel_n    = (TH1D*) InputFile->Get("n_pT_"    + ThisHistoLabel);
	      TH1D *PtRel_ntag = (TH1D*) InputFile->Get("ntag_pT_" + ThisHistoLabel + "_" + TaggerName[tg]);
	      TH1D *PtRel_p    = (TH1D*) InputFile->Get("p_pT_"    + ThisHistoLabel);
	      TH1D *PtRel_ptag = (TH1D*) InputFile->Get("ptag_pT_" + ThisHistoLabel + "_" + TaggerName[tg]);

	      for (int ptr = 1; ptr<=70; ptr++) {

		n_pT   ->SetBinContent(ptb+1, ptr, PtRel_n   ->GetBinContent(ptr));
		n_pT   ->SetBinError  (ptb+1, ptr, PtRel_n   ->GetBinError  (ptr));

		ntag_pT->SetBinContent(ptb+1, ptr, PtRel_ntag->GetBinContent(ptr));
		ntag_pT->SetBinError  (ptb+1, ptr, PtRel_ntag->GetBinError  (ptr));

		p_pT   ->SetBinContent(ptb+1, ptr, PtRel_p   ->GetBinContent(ptr));
		p_pT   ->SetBinError  (ptb+1, ptr, PtRel_p   ->GetBinError  (ptr));

		ptag_pT->SetBinContent(ptb+1, ptr, PtRel_ptag->GetBinContent(ptr));
		ptag_pT->SetBinError  (ptb+1, ptr, PtRel_ptag->GetBinError  (ptr));
 
	      }

	    }

	  } else {

	    for (int ptb = 0; ptb<nPtRelPtBins; ptb++) {

	      TString ThisHistoLabel = PtRelPtBin[ptb] + "_" + PtRelEtaBin[nb] + SystematicName[is];

	      TH1D *PtRel_n_b     = (TH1D*) InputFile->Get("n_pT_"    + ThisHistoLabel + "_b");
	      TH1D *PtRel_ntag_b  = (TH1D*) InputFile->Get("ntag_pT_" + ThisHistoLabel + "_" + TaggerName[tg] + "_b");
	      TH1D *PtRel_p_b     = (TH1D*) InputFile->Get("p_pT_"    + ThisHistoLabel + "_b");
	      TH1D *PtRel_ptag_b  = (TH1D*) InputFile->Get("ptag_pT_" + ThisHistoLabel + "_" + TaggerName[tg] + "_b");

	      TH1D *PtRel_n_cl    = (TH1D*) InputFile->Get("n_pT_"    + ThisHistoLabel + "_cl");
	      TH1D *PtRel_ntag_cl = (TH1D*) InputFile->Get("ntag_pT_" + ThisHistoLabel + "_" + TaggerName[tg] + "_cl");
	      TH1D *PtRel_p_cl    = (TH1D*) InputFile->Get("p_pT_"    + ThisHistoLabel + "_cl");
	      TH1D *PtRel_ptag_cl = (TH1D*) InputFile->Get("ptag_pT_" + ThisHistoLabel + "_" + TaggerName[tg] + "_cl");

	      for (int ptr = 1; ptr<=70; ptr++) {

		n_pT_b    ->SetBinContent(ptb+1, ptr, PtRel_n_b    ->GetBinContent(ptr));
		n_pT_b    ->SetBinError  (ptb+1, ptr, PtRel_n_b    ->GetBinError  (ptr));

		ntag_pT_b ->SetBinContent(ptb+1, ptr, PtRel_ntag_b ->GetBinContent(ptr));
		ntag_pT_b ->SetBinError  (ptb+1, ptr, PtRel_ntag_b ->GetBinError  (ptr));

		p_pT_b    ->SetBinContent(ptb+1, ptr, PtRel_p_b    ->GetBinContent(ptr));
		p_pT_b    ->SetBinError  (ptb+1, ptr, PtRel_p_b    ->GetBinError  (ptr));

		ptag_pT_b ->SetBinContent(ptb+1, ptr, PtRel_ptag_b ->GetBinContent(ptr));
		ptag_pT_b ->SetBinError  (ptb+1, ptr, PtRel_ptag_b ->GetBinError  (ptr));

		n_pT_cl   ->SetBinContent(ptb+1, ptr, PtRel_n_cl   ->GetBinContent(ptr));
		n_pT_cl   ->SetBinError  (ptb+1, ptr, PtRel_n_cl   ->GetBinError  (ptr));

		ntag_pT_cl->SetBinContent(ptb+1, ptr, PtRel_ntag_cl->GetBinContent(ptr));
		ntag_pT_cl->SetBinError  (ptb+1, ptr, PtRel_ntag_cl->GetBinError  (ptr));

		p_pT_cl   ->SetBinContent(ptb+1, ptr, PtRel_p_cl   ->GetBinContent(ptr));
		p_pT_cl   ->SetBinError  (ptb+1, ptr, PtRel_p_cl   ->GetBinError  (ptr));

		ptag_pT_cl->SetBinContent(ptb+1, ptr, PtRel_ptag_cl->GetBinContent(ptr));
		ptag_pT_cl->SetBinError  (ptb+1, ptr, PtRel_ptag_cl->GetBinError  (ptr));
		
	      }

	    }

	    n_pT   ->Add( n_pT_b    ); n_pT   ->Add( n_pT_cl    );
	    ntag_pT->Add( ntag_pT_b ); ntag_pT->Add( ntag_pT_cl );
	    p_pT   ->Add( p_pT_b    ); p_pT   ->Add( p_pT_cl    );
	    ptag_pT->Add( ptag_pT_b ); ptag_pT->Add( ptag_pT_cl );

	  }

	  TDirectory *lepton_in_jet = OutputFile->mkdir("lepton_in_jet");
	  lepton_in_jet->cd();
	  
	  n_pT   ->Write();
	  ntag_pT->Write();
	  p_pT   ->Write();
	  ptag_pT->Write();
 
	  if (dt==1) {

	    TDirectory *MCTruth = OutputFile->mkdir("MCTruth");
	    MCTruth->cd();
	    
	    n_pT_b    ->Write();
	    ntag_pT_b ->Write();
	    p_pT_b    ->Write();
	    ptag_pT_b ->Write();
	    
	    n_pT_cl   ->Write();
	    ntag_pT_cl->Write();
	    p_pT_cl   ->Write();
	    ptag_pT_cl->Write();
	  
	  }
 
	  OutputFile->Close();

	}
  
  InputFile->Close();
	  
}
