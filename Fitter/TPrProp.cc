#include "TPrProp.h"

Float_t TPrProp::Scaling(TPrProp fChProcess){

Float_t totscl1=0;

totscl1  = (Float_t)gRandom->Gaus(fChProcess.GetBR(), fChProcess.GetbranchingFractionSystError())
         * (Float_t)gRandom->Gaus(fChProcess.GetFilterEfficiency(), 0)
         * (Float_t)gRandom->Gaus(fChProcess.GetacceptanceEfficiency(), fChProcess.GetacceptanceEfficiencySystError())
         * (Float_t)gRandom->Gaus(fChProcess.GetelectronIDEfficiency(), fChProcess.GetelectronIDEfficiencySystError())
         * (Float_t)gRandom->Gaus(fChProcess.GetmuonIDEfficiency(), fChProcess.GetmuonIDEfficiencySystError())
         * (Float_t)gRandom->Gaus(fChProcess.GettopologyEfficiency(), fChProcess.GettopologyEfficiencySystError())
         * (Float_t)gRandom->Gaus(fChProcess.GettauIDEfficiency(), fChProcess.GettauIDEfficiencySystError())
         * (Float_t)gRandom->Gaus(50.0, 5.0);

//cout<<totscl1<<endl;
return(totscl1);

}
