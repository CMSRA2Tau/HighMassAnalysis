/*********************************************************
********************************************************** 
***  October 10,2010					**
***  Indara Suarez					**
***  Texas A&M University				**
***	   TSystematicStateGenerator constructor   	**
***  version 1.0    	    				**
**********************************************************
*********************************************************/
#include "TH1F.h"
#include "TPrProp.h"
#include "TChProp.h"

#ifndef __TSystematicStateGenerator_h__
#define __TSystematicStateGenerator_h__


class TSystematicStateGenerator{
private:  


public:
 TSystematicStateGenerator(){};
 virtual ~TSystematicStateGenerator(){};
 
 vector<float> GenerateSystematicState(int n);
 
 TH1F* GetNormalizedDistribution(int setting, TPrProp ProcessProps);
 Float_t SmearLumi(float lumi_mean, float lumi_err);
  
  
};

#endif
