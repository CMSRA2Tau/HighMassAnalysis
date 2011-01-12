ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#######################################
# -- DEFINE ARCH to something sensible!
#######################################
 
#
CXX           = g++
CXXFLAGS      = -g -Wall -fPIC
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared


CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

runFitter: runFitter.o TFitter.o TChProp.o TPrProp.o TSystematicStateGenerator.o
	$(CXX) $(CXXFLAGS) -c $<
	$(LD) $(LDFLAGS) -o runFitter runFitter.o TFitter.o TChProp.o TPrProp.o TSystematicStateGenerator.o $(GLIBS)

runFitter.o: runFitter.C TFitter.h
	$(CXX) $(CXXFLAGS) -c runFitter.C

TFitter.o: TFitter.cc TFitter.h 
	$(CXX) $(CXXFLAGS) -c TFitter.cc

TChProp.o: TChProp.cc TChProp.h TPrProp.h
	$(CXX) $(CXXFLAGS) -c TChProp.cc
	
TPrProp.o: TPrProp.cc TPrProp.h
	$(CXX) $(CXXFLAGS) -c TPrProp.cc 

TSystematicStateGenerator.o: TSystematicStateGenerator.cc TSystematicStateGenerator.h
	$(CXX) $(CXXFLAGS) -c TSystematicStateGenerator.cc 

