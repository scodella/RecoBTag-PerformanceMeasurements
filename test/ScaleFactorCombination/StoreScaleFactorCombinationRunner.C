void StoreScaleFactorCombinationRunner(bool Compile, string BTagger, string OP = "All", string MeasType = "comb", string StoreSystematicBreakdown = "false") {
  
  if (Compile)
    gROOT->ProcessLine(".L ScaleFactorCombination.C++");
  else 
    gSystem->Load("ScaleFactorCombination_C.so");

  TString StoreCommand =  "StoreScaleFactorCombination(\"" + BTagger + "\", \"" + OP + "\", \"" + MeasType + "\", " + StoreSystematicBreakdown + ")";
  gROOT->ProcessLine(StoreCommand);
  
}
  
