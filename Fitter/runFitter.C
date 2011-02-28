#include "TFitter.h"

int main(){
  TFitter theFitter;
  theFitter.SetNExp(0);  //if zero, no pseudodata will be generated.  Data will be used for the fit
  theFitter.SetRebinFactor(1);//(int)rebin*10 = binsize;
  theFitter.SetLumi(36.15);
  theFitter.SetLumiErr(3.95);//enter "0." if you do not wish to smear luminosity
  theFitter.NMCInt(1000); //0 if wish to run without systematics
  theFitter.SetSystSettings(0);
  //enter 0 if you want to use smearing procedure only & smear luminosity
  //enter 1 if you want to use morphing only & smear luminosity
  //enter 2 if you want to smear and morph simulatenously (default) & smear luminosity
  //enter 3 if you want to smear luminosity only (and not smear normalization due to cummulative errors)
  theFitter.SetZprimeTempName("Zprime400TemplateDirectory");
  theFitter.SetZprimeTreeName("Zprime400Tree");
  theFitter.SetSignalMassXsection(400., 4.98); // xsection in pb
  theFitter.Loop();
  return 0;
}
