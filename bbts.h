#ifndef BBTS_H
#define BBTS_H

//ROOT Includes
#include "TH1F.h"
#include "TFile.h"
#include "Riostream.h"
#include "TFractionFitter.h"
#include "TTreeReader.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFitter.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TCut.h"
#include "TChain.h"
#include "TPaveText.h"
#include "TRatioPlot.h"
//VEGAS Includes
#include "VASkyMap.h"
#include "VSOptions.hpp"
#include "VSDataConverter.hpp"
#include "VATime.h"
#include "VASimpleCuts.h"
//CPP Includes
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
//C Includes
#include <ctime>
#include <cstring>
#include <cstdlib>
//Project Includes
#include "TRatioPlot_BetterError.h"

#define timestamp                                           \
{                                                           \
  Long_t __timestamp = std::time(0);                        \
  std::cout << std::asctime(std::localtime(&__timestamp));  \
}
#define PRINTSPACE std::left << std::setw(PW)
#define NBIN 20
#define MSWLOW 0.8
#define MSWHIGH 1.3
#define PW 15

extern double ZABINS[];
extern double EBINS[];
extern int TBINS[];
extern double AZBINS[];
extern double OBINS[];

//Histograms
/*
extern TH1F* DAT_HIST;
extern TH1F* BKG_HIST;
extern TH1F* SRC_HIST;
*/
enum class Format_t{Toy, Csv, Vegas, Sample};
struct args_t {
  args_t() : format(Format_t::Toy),
  hist(0),
  output(0),
  verbosity(-1),
  bin_vars(7),
  graphics(0),
  bidir(0) {}
  Format_t format;
  int hist;
  int output;
  int verbosity;
  int bin_vars;
  int graphics;
  bool bidir;
  std::string op_info;
};
struct indices_t{
  int za;
  int e;
  int tel;
  int az;
  int off;
};
struct cuts_t{
  cuts_t() : za(0), e(0), tel(0), az(0), off(0), msw(0),
  src(0), passed(0), read(0) {}
  int za;
  int e;
  int tel;
  int az;
  int off;
  int msw;
  int src;
  int passed;
  int read;
};
struct hists_t{
  TH1F* dat_hist;
  TH1F* bkg_hist;
  TH1F* src_hist;
  TH2F* dat_2hist;
  TH2F* bkg_2hist;
};


void print_cuts(std::string action, cuts_t* cuts);







#endif