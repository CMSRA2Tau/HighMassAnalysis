#include "TFitter.h"

int main(){
  TFitter theFitter;
  theFitter.SetNExp(10);  
  theFitter.SetLumi(50.);
  theFitter.SetLumiErr(12.5);//enter "0." if you do not wish to smear luminosity
  theFitter.NMCInt(1000); //0 if wish to run without systematics
  theFitter.SetSystSettings(0);
  //enter 0 if you want to use smearing procedure only
  //enter 1 if you want to use morphing only
  //enter 2 if you want to smear and morph simulatenously (default)
  //enter 3 if you want to smear luminosity only (and not smear normalization due to cummulative errors)
  theFitter.SetSignalMassXsection(500., 1.914); // xsection in pb
  theFitter.Loop();
  return 0;
}
