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

//Bin boundaries
double ZABINS[] = {25.46,32.59,37.57,42.56,47.55};
double EBINS[] = {316,630,1259,2512,5012};
int TBINS[] = {3,4};
double AZBINS[] = {0,45,90,135,180,225,270,315,360};
double OBINS[] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};
enum class Format_t{Toy, Csv, Vegas, Sample};
struct args_t {
  args_t() : format(Format_t::Toy),
  hist(0),
  output(0),
  verbosity(-1),
  bin_vars(7),
  graphics(0) {}
  Format_t format;
  int hist;
  int output;
  int verbosity;
  int bin_vars;
  int graphics;
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










#endif