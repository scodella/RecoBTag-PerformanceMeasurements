void StoreScaleFactorCombinationRunner(bool Compile, string BTagger, string OP = "All", string MeasType = "mujets", string JetFlavour = "6", string StoreSystematicBreakdown = "NONE") {
  
  if (Compile)
    gROOT->ProcessLine(".L ScaleFactorCombination.C++");
  else 
    gSystem->Load("ScaleFactorCombination_C.so");

  TString StoreCommand =  "StoreScaleFactorCombination(\"" + BTagger + "\", \"" + OP + "\", \"" + MeasType + "\", " + JetFlavour + ", \"" + StoreSystematicBreakdown + "\")";
  
  gROOT->ProcessLine(StoreCommand);
  
}

