void ScaleFactorCombinationRunner(bool Compile, TString Algo, TString WP) {
  
  if (Compile)
    gROOT->ProcessLine(".L ScaleFactorCombination.C++");
  else 
    gSystem->Load("ScaleFactorCombination_C.so");

  if (WP=="Loose"  || WP=="All") gROOT->ProcessLine(".x ScaleFactorCombinationExec.C(\"" + Algo + "\", \"Loose\")");
  if (WP=="Medium" || WP=="All") gROOT->ProcessLine(".x ScaleFactorCombinationExec.C(\"" + Algo + "\", \"Medium\")");
  if (WP=="Tight"  || WP=="All") gROOT->ProcessLine(".x ScaleFactorCombinationExec.C(\"" + Algo + "\", \"Tight\")");

}
  
