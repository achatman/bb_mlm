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
std::string OUTPATH;
std::string data_excl_file = "data_src.txt";
std::string bkg_excl_file = "bkg_src.txt";
//Histograms
TH1F* DAT_HIST;
TH1F* BKG_HIST;
TH1F* SRC_HIST;
//Bin boundaries
double ZABINS[] = {25.46,32.59,37.57,42.56,47.55};
double EBINS[] = {316,630,1259,2512,5012};
int TBINS[] = {3,4};
double AZBINS[] = {0,45,90,135,180,225,270,315,360};
double OBINS[] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};
enum class Format_t{Toy, Csv, Vegas};
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

//extra function declarations
int parse_command_line(int argc, char* argv[], args_t* args);
void printRawData();
void histogram_raw_data(indices_t ins);
void histogram_fit_data(double fracs[6], indices_t ins);
void calculate_errors(double Pb, double Ps, double sigma_Pb, double sigma_Ps, indices_t ins, double alpha);
void print_cuts(std::string action, cuts_t* cuts);
void fit_manual_bin();
void map_likelihood(double Pb, double Ps, std::string title_tag, indices_t ins, args_t args);


//Li Ma Significance Calculation
//Eq 17, Li Ma (1983)
double lima_sig(double Pb, double Ps){
  double alpha = DAT_HIST->Integral() / BKG_HIST->Integral();
  double pb = Pb * DAT_HIST->Integral() / BKG_HIST->Integral();
  double ps = Ps * DAT_HIST->Integral() / SRC_HIST->Integral();
  double N_on = ps * SRC_HIST->Integral();
  double N_off = pb * BKG_HIST->Integral();
  double S = TMath::Sqrt(2)
           * TMath::Sqrt(
              N_on * TMath::Log( ((1 + alpha) / alpha) * (N_on / (N_on + N_off)) )
            + N_off * TMath::Log( (1 + alpha) * (N_off / (N_on + N_off)) )
            );
  return S;
}

void loadsrc_csv(indices_t ins, args_t args){
  std::cout << "Loading Source" << std::endl;
  SRC_HIST =  new TH1F("SrcHist" , "SRC" , NBIN, MSWLOW, MSWHIGH);
  timestamp;
  std::ifstream srcif("src.list");
  std::string srcstr;
  std::stringstream srcPath;
  std::getline(srcif, srcstr);
  Int_t index;
  index = srcstr.find("MCSRC_");
  srcPath << srcstr.substr(0, index+5);
  if(args.bin_vars & 1){
    srcPath << "_ZA_" << ZABINS[ins.za] << "_" << ZABINS[ins.za + 1];
  }
  if(args.bin_vars & 2){
    srcPath << "_E" << ins.e;
  }
  if(args.bin_vars & 4){
    srcPath << "_T" << TBINS[ins.tel];
  }
  if(args.bin_vars & 8){
    srcPath << "_A" << ins.az;
  }
  if(args.bin_vars & 16){
    srcPath << "_O" << ins.off;
  }
  srcPath << ".csv";

  std::string line;
  std::ifstream srcdatin(srcPath.str().c_str());
  while(std::getline(srcdatin, line)){
    index = line.find(",");
    line = line.substr(index+1,line.size() - index);
    SRC_HIST->Fill(atof(line.c_str()));
  }
  srcdatin.close();
  std::cout << SRC_HIST->Integral()
            << " events loaded for source." << std::endl;
}

void loadData(indices_t ins, args_t args, double* alpha){
  std::string line;
  TFile *f;
  std::ifstream data_flist("data.list");
  std::ifstream bkg_flist("bkg.list");

  TH1::SetDefaultSumw2();
  DAT_HIST = new TH1F("DataHist", "Data", NBIN, MSWLOW, MSWHIGH);
  BKG_HIST = new TH1F("BkgHist" , "BKG" , NBIN, MSWLOW, MSWHIGH);

  //Check each point against cuts
  cuts_t data_cuts;
  double livetime_data = 0, livetime_bkg = 0;
  int data_overflow = 0, bkg_overflow = 0;

  //Construct Data Histogram
  std::cout << "Loading Data" << std::endl;
  timestamp;
  while(std::getline(data_flist, line)){
    f = TFile::Open(line.c_str());
    if(!f){
      std::cout << "Cannot open " << line << ". Moving swiftly onward." << std::endl;
      continue;
    }

    TTreeReader reader((TTree*)f->Get("ToyMCStudies/ToyMCSample_1/DataTree_1"));
    TTreeReaderValue<double> za(reader, "Elevation_Deg");
    TTreeReaderValue<double> energy(reader, "EnergyLower_TeV");
    TTreeReaderValue<double> telescopes(reader, "TelescopeMultiplicity");
    TTreeReaderValue<double> msw(reader, "MeanScaledWidth");
    TTreeReaderValue<double> live(reader, "Livetime_Seconds");
    TTreeReaderValue<double> az(reader, "Azimuth_Deg");
    TTreeReaderValue<double> event_dec(reader, "EventDec_Deg");
    TTreeReaderValue<double> event_ra(reader, "EventRA_Deg");
    TTreeReaderValue<double> tracking_dec(reader, "TrackingDec_Deg");
    TTreeReaderValue<double> tracking_ra(reader, "TrackingRA_Deg");

    while(reader.Next()){
      data_cuts.read++;
      //Telescope Cut
      if(args.bin_vars & 4){
        if(*telescopes != TBINS[ins.tel]){
          data_cuts.tel++;
          continue;
        }
      }

      //Energy Cut
      if(args.bin_vars & 2){
        double en = *energy * 1000;
        if(en > EBINS[ins.e+1] || en < EBINS[ins.e]){
          data_cuts.e++;
          continue;
        }
      }

      //ZA Cut
      if(args.bin_vars & 1){
        if(90-(*za) > ZABINS[ins.za + 1] || 90-(*za) < ZABINS[ins.za]){
          data_cuts.za++;
          continue;
        }
      }

      //MSW Cut
      if(*msw < MSWLOW || *msw > MSWHIGH){
        data_cuts.msw++;
        continue;
      }

      //AZ Cut
      if(args.bin_vars & 8){
        if(*az > AZBINS[ins.az + 1] || *az < AZBINS[ins.az]){
          data_cuts.az++;
          continue;
        }
      }
      //Offset Cut
      if(args.bin_vars & 16){
        VACoordinatePair shower_coords = VACoordinatePair(
          *event_dec, *event_ra, VACoordinates::J2000, VACoordinates::Deg
        );
        VACoordinatePair tracking_coords = VACoordinatePair(
          *tracking_dec, *tracking_ra, VACoordinates::J2000, VACoordinates::Deg
        );
        double offset = tracking_coords.angularSeparation_Deg(shower_coords);
        if(offset > 2.0) data_overflow++;
        if(offset > OBINS[ins.off + 1] || offset < OBINS[ins.off]){
          data_cuts.off++;
          continue;
        }
      }

      DAT_HIST->Fill(*msw);
      livetime_data += *live;
    }
    f->Close();
  }
  data_cuts.passed = DAT_HIST->Integral();
  std::cout << data_cuts.passed << " passed cuts." << std::endl;
  std::cout << data_cuts.read << " points read." << std::endl;
  std::cout << data_cuts.tel << " failed tel cut." << std::endl;
  std::cout << data_cuts.e << " failed energy cut." << std::endl;
  std::cout << data_cuts.za << " failed za cut." << std::endl;
  std::cout << data_cuts.msw << " failed msw cut." << std::endl;
  std::cout << data_cuts.az << " failed az cut." << std::endl;
  std::cout << data_cuts.off << " failed off cut." << std::endl;

  //Construct Background Histogram
  cuts_t bkg_cuts;
  std::cout << "Loading Background" << std::endl;
  timestamp;
  while(std::getline(bkg_flist,line)){
    f = TFile::Open(line.c_str());
    if(!f){
      std::cout << "Cannot open " << line << ". Moving swiftly onward." << std::endl;
      continue;
    }

    TTreeReader reader((TTree*)f->Get("ToyMCStudies/ToyMCSample_1/DataTree_1"));
    TTreeReaderValue<double> za(reader, "Elevation_Deg");
    TTreeReaderValue<double> energy(reader, "EnergyLower_TeV");
    TTreeReaderValue<double> telescopes(reader, "TelescopeMultiplicity");
    TTreeReaderValue<double> msw(reader, "MeanScaledWidth");
    TTreeReaderValue<double> live(reader, "Livetime_Seconds");
    TTreeReaderValue<double> az(reader, "Azimuth_Deg");
    TTreeReaderValue<double> event_dec(reader, "EventDec_Deg");
    TTreeReaderValue<double> event_ra(reader, "EventRA_Deg");
    TTreeReaderValue<double> tracking_dec(reader, "TrackingDec_Deg");
    TTreeReaderValue<double> tracking_ra(reader, "TrackingRA_Deg");

    while(reader.Next()){
      bkg_cuts.read++;
      //Telescope Cut
      if(args.bin_vars & 4){
        if(*telescopes != TBINS[ins.tel]){
          bkg_cuts.tel++;
          continue;
        }
      }

      //Energy Cut
      if(args.bin_vars & 2){
        double en = *energy * 1000;
        if(en > EBINS[ins.e+1] || en < EBINS[ins.e]){
          bkg_cuts.e++;
          continue;
        }
      }

      //ZA Cut
      if(args.bin_vars & 1){
        if(90-(*za) > ZABINS[ins.za + 1] || 90-(*za) < ZABINS[ins.za]){
          bkg_cuts.za++;
          continue;
        }
      }

      //MSW Cut
      if(*msw < MSWLOW || *msw > MSWHIGH){
        bkg_cuts.msw++;
        continue;
      }

      //AZ Cut
      if(args.bin_vars & 8){
        if(*az > AZBINS[ins.az + 1] || *az < AZBINS[ins.az]){
          bkg_cuts.az++;
          continue;
        }
      }
      //Offset Cut
      if(args.bin_vars & 16){
        VACoordinatePair shower_coords = VACoordinatePair(
          *event_dec, *event_ra, VACoordinates::J2000, VACoordinates::Deg
        );
        VACoordinatePair tracking_coords = VACoordinatePair(
          *tracking_dec, *tracking_ra, VACoordinates::J2000, VACoordinates::Deg
        );
        double offset = tracking_coords.angularSeparation_Deg(shower_coords);
        if(offset > 2.0) bkg_overflow++;
        if(offset > OBINS[ins.off + 1] || offset < OBINS[ins.off]){
          bkg_cuts.off++;
          continue;
        }
      }

      BKG_HIST->Fill(*msw);
      livetime_bkg += *live;
    }
  }
  bkg_cuts.passed = BKG_HIST->Integral();
  std::cout << bkg_cuts.passed << " passed cuts." << std::endl;
  std::cout << bkg_cuts.read << " cuts read." << std::endl;
  std::cout << bkg_cuts.tel << " failed tel cut." << std::endl;
  std::cout << bkg_cuts.e << " failed energy cut." << std::endl;
  std::cout << bkg_cuts.za << " failed za cut." << std::endl;
  std::cout << bkg_cuts.msw << " failed msw cut." << std::endl;
  std::cout << bkg_cuts.az << " failed az cut." << std::endl;
  std::cout << bkg_cuts.off << " failed off cut." << std::endl;

  //Construct Source histogram
  loadsrc_csv(ins, args);
  if(args.output & 8){
    print_cuts("data", &data_cuts);
    print_cuts("bkg", &bkg_cuts);
    print_cuts("src", 0);
  }
  std::cout << "Data Overflow: " << data_overflow << std::endl;
  std::cout << "Bkg Overflow: " << bkg_overflow << std::endl;
  *alpha = livetime_data / livetime_bkg;
}

void loadData_csv(indices_t ins, args_t args, double* alpha){
    std::ifstream dataif("data.list");
    std::ifstream bkgif("bkg.list");
    std::ifstream srcif("src.list");
    std::string datastr, bkgstr, srcstr;
    std::stringstream dataPath, bkgPath, srcPath;
    std::getline(dataif, datastr);
    std::getline(bkgif,  bkgstr );
    std::getline(srcif,  srcstr );
    Int_t index;
    index = datastr.find("ZA");
    dataPath << datastr.substr(0, index) << "ZA_" << ZABINS[ins.za] << "_" << ZABINS[ins.za + 1] << "_E" << ins.e << "_T" << TBINS[ins.tel] << ".csv";
    index = bkgstr.find("ZA");
    bkgPath << bkgstr.substr(0, index) << "ZA_" << ZABINS[ins.za] << "_" << ZABINS[ins.za + 1] << "_E" << ins.e << "_T" << TBINS[ins.tel] << ".csv";
    index = srcstr.find("ZA");
    srcPath << srcstr.substr(0, index) << "ZA_" << ZABINS[ins.za] << "_" << ZABINS[ins.za + 1] << "_E" << ins.e << "_T" << TBINS[ins.tel] << ".csv";

    //Create TTrees and open files
    TTree *dataTree= new TTree("DATA","datatree");
    TTree *bkgTree = new TTree("BKG","bkgtree");
    TTree *srcTree = new TTree("SRC","srctree");

    //If any file fails to open, execution is aborted
    if(!(dataTree->ReadFile(dataPath.str().c_str(),"ZA:MSW",','))){
        throw 2;
    }
    if(!(bkgTree->ReadFile(bkgPath.str().c_str(),"ZA:MSW",','))){
        throw 1;
    }
    if(!(srcTree->ReadFile(srcPath.str().c_str(),"ZA:MSW",','))){
        throw 3;
    }

    TH1::SetDefaultSumw2();
    DAT_HIST = new TH1F("DataHist", "Data", NBIN, MSWLOW, MSWHIGH);
    BKG_HIST =  new TH1F("BkgHist" , "BKG" , NBIN, MSWLOW, MSWHIGH);
    SRC_HIST =  new TH1F("SrcHist" , "SRC" , NBIN, MSWLOW, MSWHIGH);

    TTreeReader dataReader(dataTree);
    TTreeReaderValue<float> datamsw(dataReader,"MSW");
    TTreeReader bkgReader(bkgTree);
    TTreeReaderValue<float> bkgmsw(bkgReader,"MSW");
    TTreeReader srcReader(srcTree);
    TTreeReaderValue<float> srcmsw(srcReader,"MSW");

    while(dataReader.Next()){
        DAT_HIST->Fill(*datamsw);
    }
    while(bkgReader.Next()){
        BKG_HIST->Fill(*bkgmsw);
    }
    while(srcReader.Next()){
        SRC_HIST->Fill(*srcmsw);
    }

    std::cout << "Data: " << DAT_HIST->Integral() << " counts loaded." << std::endl;
    std::cout << "Bkg: " << BKG_HIST->Integral() << " counts loaded." << std::endl;
    std::cout << "Src: " << SRC_HIST->Integral() << " counts loaded." << std::endl;
}

void loadData_vegas(indices_t ins, args_t args, double* alpha){
  //Check for source exclusion files
  if(access(data_excl_file.c_str(), F_OK)){
    std::cerr << "No source exclusion file found at "
              << data_excl_file
              << ". This may cause problems."
              << std::endl;
  }
  if(access(bkg_excl_file.c_str(), F_OK)){
    std::cerr << "No source exclusion file found at "
              << bkg_excl_file
              << ". This may cause problems."
              << std::endl;
  }
  int data_overflow = 0, bkg_overflow = 0;
  typedef struct SourceCut
  {
    VACoordinatePair coords;
    double cut_rad;
    SourceCut(std::vector<std::string> line_in) :
      coords(std::stod(line_in[1]),std::stod(line_in[2]),VACoordinates::J2000,VACoordinates::Deg),
      cut_rad(std::stod(line_in[3]))
      {}

    bool InsideExclRadius(VACoordinatePair evt_coords) {
      double dist = coords.angularSeparation_Deg(evt_coords);
      return (dist < cut_rad);
    }
  } SourceCut_t;

  std::vector<SourceCut_t> data_sourcecuts;
  std::vector<SourceCut_t> bkg_sourcecuts;
  std::string line;
  std::vector<std::string> line_fields;
  std::ifstream data_srclist(data_excl_file.c_str());
  std::ifstream bkg_srclist(bkg_excl_file.c_str());
  while(std::getline(data_srclist, line)){
    boost::split(line_fields, line, boost::is_any_of(","));
    data_sourcecuts.emplace_back(line_fields);
    line_fields.clear();
  }
  while(std::getline(bkg_srclist, line)){
    boost::split(line_fields, line, boost::is_any_of(","));
    bkg_sourcecuts.emplace_back(line_fields);
    line_fields.clear();
  }

  //Setup data chains
  TChain* data_chain = new TChain("SelectedEvents/CombinedEventsTree");
  TChain* bkg_chain = new TChain("SelectedEvents/CombinedEventsTree");
  std::ifstream data_flist("data.list");
  std::ifstream bkg_flist("bkg.list");

  while(std::getline(data_flist, line)){
    data_chain->Add(line.c_str());
  }
  while(std::getline(bkg_flist, line)){
    bkg_chain->Add(line.c_str());
  }
  VAShowerData* data_shower = new VAShowerData;
  VAShowerData* bkg_shower = new VAShowerData;
  data_chain->SetBranchAddress("S", &data_shower);
  bkg_chain->SetBranchAddress("S", &bkg_shower);

  //Setup histograms
  TH1::SetDefaultSumw2();
  DAT_HIST = new TH1F("DataHist", "Data", NBIN, MSWLOW, MSWHIGH);
  BKG_HIST =  new TH1F("BkgHist" , "BKG", NBIN, MSWLOW, MSWHIGH);

  cuts_t data_cuts, bkg_cuts;
  timestamp;
  std::cout << "Loading Data" << std::endl;
  //Construct Data Histogram
  for(int i = 0; i < data_chain->GetEntries(); i++){
    data_chain->GetEntry(i);
    data_cuts.read++;
    //Source Cut
    Double_t eventRA = data_shower->fDirectionRA_J2000_Rad * TMath::RadToDeg();
    Double_t eventDec = data_shower->fDirectionDec_J2000_Rad * TMath::RadToDeg();
    VACoordinatePair eventCoord = VACoordinatePair(eventRA,eventDec,VACoordinates::J2000,VACoordinates::Deg);
    bool fail_src_cuts;
    for(auto &it: data_sourcecuts){
      fail_src_cuts = it.InsideExclRadius(eventCoord);
      if(fail_src_cuts) break;
    }
    if(fail_src_cuts){
      data_cuts.src++;
      continue;
    }

    //Tel cut
    if(args.bin_vars & 4){
      auto tels_used = data_shower->fTelUsedInReconstruction;
      int eventNTel = count(tels_used.begin(), tels_used.end(), 1);
      if(eventNTel != TBINS[ins.tel]){
        data_cuts.tel++;
        continue;
      }
    }

    //Energy cut
    if(args.bin_vars & 2){
      if(data_shower->fEnergy_GeV > EBINS[ins.e + 1] || data_shower->fEnergy_GeV < EBINS[ins.e]){
        data_cuts.e++;
        continue;
      }
    }

    //ZA cut
    if(args.bin_vars & 1){
      Double_t ZA = 90.0 - (data_shower->fDirectionElevation_Rad * TMath::RadToDeg());
      if(ZA > ZABINS[ins.za + 1] || ZA < ZABINS[ins.za]){
        data_cuts.za++;
        continue;
      }
    }

    //MSW cut
    if(data_shower->fMSW < MSWLOW || data_shower->fMSW > MSWHIGH){
      data_cuts.msw++;
      continue;
    }

    //AZ cut
    if(args.bin_vars & 8){
      double az = data_shower->fDirectionAzimuth_Rad * TMath::RadToDeg();
      if(az > AZBINS[ins.az + 1] || az < AZBINS[ins.az]){
        data_cuts.az++;
        continue;
      }
    }
    //Offset cut
    if(args.bin_vars & 16){
      VACoordinatePair shower_coords = VACoordinatePair(
        data_shower->fDirectionRA_J2000_Rad,
        data_shower->fDirectionDec_J2000_Rad,
        VACoordinates::J2000,
        VACoordinates::Rad
      );
      VACoordinatePair tracking_coords = VACoordinatePair(
        data_shower->fArrayTrackingRA_J2000_Rad,
        data_shower->fArrayTrackingDec_J2000_Rad,
        VACoordinates::J2000,
        VACoordinates::Rad
      );
      double offset = tracking_coords.angularSeparation_Deg(shower_coords);
      if(offset > 2.0) data_overflow++;
      if(offset > OBINS[ins.off + 1] || offset < OBINS[ins.off]){
        data_cuts.off++;
        continue;
      }
    }

    DAT_HIST->Fill(data_shower->fMSW);
  }
  data_cuts.passed = DAT_HIST->Integral();
  std::cout << data_cuts.passed << " passed cuts." << std::endl;
  std::cout << data_cuts.read << " cuts read." << std::endl;
  std::cout << data_cuts.src << " failed src cut." << std::endl;
  std::cout << data_cuts.tel << " failed tel cut." << std::endl;
  std::cout << data_cuts.e << " failed energy cut." << std::endl;
  std::cout << data_cuts.za << " failed za cut." << std::endl;
  std::cout << data_cuts.msw << " failed msw cut." << std::endl;
  std::cout << data_cuts.az << " failed az cut." << std::endl;
  std::cout << data_cuts.off << " failed off cut." << std::endl;

  //Construct Bkg Histogram
  std::cout << "Loading Background" << std::endl;
  for(int i = 0; i < bkg_chain->GetEntries(); i++){
    bkg_chain->GetEntry(i);
    //Source Cuts
    bkg_cuts.read++;
    Double_t eventRA = bkg_shower->fDirectionRA_J2000_Rad * TMath::RadToDeg();
    Double_t eventDec = bkg_shower->fDirectionDec_J2000_Rad * TMath::RadToDeg();
    VACoordinatePair eventCoord = VACoordinatePair(eventRA,eventDec,VACoordinates::J2000,VACoordinates::Deg);
    bool fail_src_cuts;
    for(auto &it: bkg_sourcecuts){
      fail_src_cuts = it.InsideExclRadius(eventCoord);
      if(fail_src_cuts) break;
    }
    if(fail_src_cuts){
      bkg_cuts.src++;
      continue;
    }

    //Tel cut
    if(args.bin_vars & 4){
      auto tels_used = bkg_shower->fTelUsedInReconstruction;
      int eventNTel = count(tels_used.begin(), tels_used.end(), 1);
      if(eventNTel != TBINS[ins.tel]){
        bkg_cuts.tel++;
        continue;
      }
    }

    //Energy cut
    if(args.bin_vars & 2){
      if(bkg_shower->fEnergy_GeV > EBINS[ins.e + 1] || bkg_shower->fEnergy_GeV < EBINS[ins.e]){
        bkg_cuts.e++;
        continue;
      }
    }

    //ZA cut
    if(args.bin_vars & 1){
      Double_t ZA = 90.0 - (bkg_shower->fDirectionElevation_Rad * TMath::RadToDeg());
      if(ZA > ZABINS[ins.za + 1] || ZA < ZABINS[ins.za]){
        bkg_cuts.za++;
        continue;
      }
    }

    //MSW cut
    if(bkg_shower->fMSW < MSWLOW || bkg_shower->fMSW > MSWHIGH){
      bkg_cuts.msw++;
      continue;
    }

    //AZ cut
    if(args.bin_vars & 8){
      double az = bkg_shower->fDirectionAzimuth_Rad * TMath::RadToDeg();
      if(az > AZBINS[ins.az + 1] || az < AZBINS[ins.az]){
        bkg_cuts.az++;
        continue;
      }
    }
    //Offset cut
    if(args.bin_vars & 16){
      VACoordinatePair shower_coords = VACoordinatePair(
        bkg_shower->fDirectionRA_J2000_Rad,
        bkg_shower->fDirectionDec_J2000_Rad,
        VACoordinates::J2000,
        VACoordinates::Rad
      );
      VACoordinatePair tracking_coords = VACoordinatePair(
        bkg_shower->fArrayTrackingRA_J2000_Rad,
        bkg_shower->fArrayTrackingDec_J2000_Rad,
        VACoordinates::J2000,
        VACoordinates::Rad
      );
      double offset = tracking_coords.angularSeparation_Deg(shower_coords);
      if(offset > 2.0) bkg_overflow++;
      if(offset > OBINS[ins.off + 1] || offset < OBINS[ins.off]){
        bkg_cuts.off++;
        continue;
      }
    }

    BKG_HIST->Fill(bkg_shower->fMSW);
  }
  bkg_cuts.passed = BKG_HIST->Integral();
  std::cout << bkg_cuts.passed << " passed cuts." << std::endl;
  std::cout << bkg_cuts.read << " cuts read." << std::endl;
  std::cout << bkg_cuts.src << " failed src cut." << std::endl;
  std::cout << bkg_cuts.tel << " failed tel cut." << std::endl;
  std::cout << bkg_cuts.e << " failed energy cut." << std::endl;
  std::cout << bkg_cuts.za << " failed za cut." << std::endl;
  std::cout << bkg_cuts.msw << " failed msw cut." << std::endl;
  std::cout << bkg_cuts.az << " failed az cut." << std::endl;
  std::cout << bkg_cuts.off << " failed off cut." << std::endl;

  //Construct Source Histogram
  loadsrc_csv(ins, args);

  if(args.output & 8){
    print_cuts("data", &data_cuts);
    print_cuts("bkg", &bkg_cuts);
    print_cuts("src", 0);
  }

  std::cout << "Data Offset Overflow: " << data_overflow << std::endl;
  std::cout << "Bkg Offset Overflow: " << bkg_overflow << std::endl;

  //TODO
  *alpha = 1.0;
}

double nosrc_noBB(double Pb, bool print = false, TH1F* F = 0){
  int i;
  double pb;
  double di, bi, fi;
  double dtot = 0, btot = 0, ftot = 0;
  std::ofstream f;
  double lnL = 0;
  pb = Pb * DAT_HIST->Integral() / BKG_HIST->Integral();
  if(print){
    std::stringstream path;
    path << "BIN_src0_bb0_" << OUTPATH << ".txt";
    f.open(path.str());
    f.precision(5);
    f << "Bin: " << OUTPATH << std::endl;
    f << std::scientific;
    f << "bkgfrac    " << Pb << std::endl;
    f << "srcfrac    -----" << std::endl;
    f << PRINTSPACE << "Counts:";
    f << PRINTSPACE << "Data";
    f << PRINTSPACE << "Bkg";
    f << std::endl << std::defaultfloat;
    f << PRINTSPACE << "";
    f << PRINTSPACE << DAT_HIST->Integral();
    f << PRINTSPACE << BKG_HIST->Integral();
    f << std::endl << std::endl;
    f << PRINTSPACE << "i";
    f << PRINTSPACE << "di";
    f << PRINTSPACE << "fi";
    f << PRINTSPACE << "bi" << std::endl;
  }
  for(i = 1; i <= NBIN; i++){
    di = DAT_HIST->GetBinContent(i);
    bi = BKG_HIST->GetBinContent(i);
    fi = pb * bi;
    if(di == 0) lnL -= fi;
    else        lnL += di * TMath::Log(fi) - fi;
    dtot += di; btot += bi; ftot += fi;
    if(print){
      f << PRINTSPACE << i;
      f << PRINTSPACE << di;
      f << PRINTSPACE << fi;
      f << PRINTSPACE << bi << std::endl;
    }
    if(F){
      F->Sumw2(false);
      F->SetBinContent(i, fi);
      F->Sumw2(true);
    }
  }
  if(print){
    f << PRINTSPACE << "Sum";
    f << PRINTSPACE << dtot;
    f << PRINTSPACE << ftot;
    f << PRINTSPACE << btot << std::endl;
    f.close();
  }
  return -lnL;
}

double src_noBB(double Pb, double Ps, bool print = false, TH1F* F = 0){
  int i;
  double pb, ps;
  double di, bi, si, fi;
  double dtot = 0, btot = 0, ftot = 0, stot = 0;
  std::ofstream f;
  double lnL = 0;
  pb = Pb * DAT_HIST->Integral() / BKG_HIST->Integral();
  ps = Ps * DAT_HIST->Integral() / SRC_HIST->Integral();
  if(print){
    std::stringstream path;
    path << "BIN_src1_bb0_" << OUTPATH << ".txt";
    f.open(path.str());
    f.precision(5);
    f << "Bin: " << OUTPATH << std::endl;
    f << std::scientific;
    f << "bkgfrac    " << Pb << std::endl;
    f << "srcfrac    " << Ps << std::endl;
    f << PRINTSPACE << "Counts:";
    f << PRINTSPACE << "Data";
    f << PRINTSPACE << "Bkg";
    f << PRINTSPACE << "Src";
    f << std::endl << std::defaultfloat;
    f << PRINTSPACE << "";
    f << PRINTSPACE << DAT_HIST->Integral();
    f << PRINTSPACE << BKG_HIST->Integral();
    f << PRINTSPACE << SRC_HIST->Integral();
    f << std::endl << std::endl;
    f << PRINTSPACE << "i";
    f << PRINTSPACE << "di";
    f << PRINTSPACE << "fi";
    f << PRINTSPACE << "bi";
    f << PRINTSPACE << "si" << std::endl;
  }
  for(i = 1; i <= NBIN; i++){
    di = DAT_HIST->GetBinContent(i);
    bi = BKG_HIST->GetBinContent(i);
    si = SRC_HIST->GetBinContent(i);
    fi = pb * bi + ps * si;
    if(di == 0) lnL -= fi;
    else        lnL += di * TMath::Log(fi) - fi;

    dtot += di; btot += bi; ftot += fi; stot += si;
    if(print){
      f << PRINTSPACE << i;
      f << PRINTSPACE << di;
      f << PRINTSPACE << fi;
      f << PRINTSPACE << bi;
      f << PRINTSPACE << si << std::endl;
    }
    if(F){
      F->Sumw2(false);
      F->SetBinContent(i, fi);
      F->Sumw2(true);
    }
  }
  if(print){
    f << PRINTSPACE << "Sum";
    f << PRINTSPACE << dtot;
    f << PRINTSPACE << ftot;
    f << PRINTSPACE << btot;
    f << PRINTSPACE << stot << std::endl;
    f.close();
  }
  return -lnL;
}

double nosrc_BB(double Pb, bool print = false, TH1F* F = 0, TH1F* B = 0){
  int i;
  double pb;
  double di, bi, Bi, fi;
  double Y, Z, ti;
  double dtot = 0, btot = 0, ftot = 0, Btot = 0;
  std::ofstream f;
  double lnL = 0;
  pb = Pb * DAT_HIST->Integral() / BKG_HIST->Integral();
  if(print){
    std::stringstream path;
    path << "BIN_src0_bb1_" << OUTPATH << ".txt";
    f.open(path.str());
    f.precision(5);
    f << "Bin: " << OUTPATH << std::endl;
    f << std::scientific;
    f << "bkgfrac    " << Pb << std::endl;
    f << "srcfrac    -----" << std::endl;
    f << PRINTSPACE << "Counts:";
    f << PRINTSPACE << "Data";
    f << PRINTSPACE << "Bkg";
    f << std::endl << std::defaultfloat;
    f << PRINTSPACE << "";
    f << PRINTSPACE << DAT_HIST->Integral();
    f << PRINTSPACE << BKG_HIST->Integral();
    f << std::endl << std::endl;
    f << PRINTSPACE << "i";
    f << PRINTSPACE << "di";
    f << PRINTSPACE << "fi";
    f << PRINTSPACE << "bi";
    f << PRINTSPACE << "Bi" << std::endl;
  }
  for(i = 1; i <= NBIN; i++){
    di = DAT_HIST->GetBinContent(i);
    bi = BKG_HIST->GetBinContent(i);

    Y = pb * di + pb * bi;
    Z = di - pb * bi;
    ti = - Z / Y;
    if(di == 0) ti = 1; //because floating point
    Bi = bi / (pb * ti + 1);
    //Check for zero in background
    if(bi == 0 && (di / (1 + pb)) > 0){
      ti = - 1 / pb;
      Bi = di / (1 + pb);
    }
    fi = pb * Bi;
    if(di == 0) lnL -= fi;
    else        lnL += di * TMath::Log(fi) - fi;
    if(bi == 0) lnL -= Bi;
    else        lnL += bi * TMath::Log(Bi) - Bi;

    dtot += di; btot += bi; ftot += fi; Btot += Bi;
    if(print){
      f << PRINTSPACE << i;
      f << PRINTSPACE << di;
      f << PRINTSPACE << fi;
      f << PRINTSPACE << bi;
      f << PRINTSPACE << Bi << std::endl;
    }
    if(F){
      F->Sumw2(false);
      F->SetBinContent(i, fi);
      F->Sumw2(true);
    }
    if(B){
      B->Sumw2(false);
      B->SetBinContent(i, Bi);
      B->Sumw2(true);
    }
  }
  if(print){
    f << PRINTSPACE << "Sum";
    f << PRINTSPACE << dtot;
    f << PRINTSPACE << ftot;
    f << PRINTSPACE << btot;
    f << PRINTSPACE << Btot << std::endl;
    f.close();
  }
  return -lnL;
}

double src_BB(double Pb, double Ps, bool print = false, TH1F* F = 0, TH1F* B = 0){
  int i;
  double pb, ps;
  double di, bi, si, Bi, fi;
  double dtot = 0, btot = 0, ftot = 0, stot = 0, Btot = 0;
  std::ofstream f;
  double X, Y, Z, ti;
  double lnL = 0;
  pb = Pb * DAT_HIST->Integral() / BKG_HIST->Integral();
  ps = Ps * DAT_HIST->Integral() / SRC_HIST->Integral();
  if(print){
    std::stringstream path;
    path << "BIN_src1_bb1_" << OUTPATH << ".txt";
    f.open(path.str());
    f.precision(5);
    f << "Bin: " << OUTPATH << std::endl;
    f << std::scientific;
    f << "bkgfrac    " << Pb << std::endl;
    f << "srcfrac    " << Ps << std::endl;
    f << PRINTSPACE << "Counts:";
    f << PRINTSPACE << "Data";
    f << PRINTSPACE << "Bkg";
    f << PRINTSPACE << "Src";
    f << std::endl << std::defaultfloat;
    f << PRINTSPACE << "";
    f << PRINTSPACE << DAT_HIST->Integral();
    f << PRINTSPACE << BKG_HIST->Integral();
    f << PRINTSPACE << SRC_HIST->Integral();
    f << std::endl << std::endl;
    f << PRINTSPACE << "i";
    f << PRINTSPACE << "di";
    f << PRINTSPACE << "fi";
    f << PRINTSPACE << "bi";
    f << PRINTSPACE << "Bi";
    f << PRINTSPACE << "si" << std::endl;
  }
  for(i = 1; i <= NBIN; i++){
    di = DAT_HIST->GetBinContent(i);
    bi = BKG_HIST->GetBinContent(i);
    si = SRC_HIST->GetBinContent(i);

    X = pb * ps * si;
    Y = pb * di + pb * bi + ps * si - pb * ps * si;
    Z = di - pb * bi - ps * si;
    ti = (-Y + TMath::Sqrt( TMath::Power(Y, 2) - 4 * X * Z)) / (2 * X);
    if(di == 0) ti = 1; //because floating point
    Bi = bi / (pb * ti + 1);
    //Check for zero in background
    if(bi == 0 && pb > ps && (di / (1 + pb) - (ps * si) / (pb - ps)) > 0){
      ti = - 1 / pb;
      Bi = di / (1 + pb) - (ps * si) / (pb - ps);
    }
    //Check for zero in source
    if(si == 0 && ps > pb && (di / (1 + ps) - (pb * bi) / (ps - pb)) > 0){
      ti = - 1 / ps;
      Bi = bi / (pb * ti + 1);
    }

    fi = pb * Bi + ps * si;
    if(di == 0) lnL -= fi;
    else        lnL += di * TMath::Log(fi) - fi;
    if(bi == 0) lnL -= Bi;
    else        lnL += bi * TMath::Log(Bi) - Bi;

    dtot += di; btot += bi; ftot += fi; stot += si; Btot += Bi;
    if(print){
      f << PRINTSPACE << i;
      f << PRINTSPACE << di;
      f << PRINTSPACE << fi;
      f << PRINTSPACE << bi;
      f << PRINTSPACE << Bi;
      f << PRINTSPACE << si << std::endl;
    }
    if(F){
      F->Sumw2(false);
      F->SetBinContent(i, fi);
      F->Sumw2(true);
    }
    if(B){
      B->Sumw2(false);
      B->SetBinContent(i, Bi);
      B->Sumw2(true);
    }
  }
  if(print){
    f << PRINTSPACE << "Sum";
    f << PRINTSPACE << dtot;
    f << PRINTSPACE << ftot;
    f << PRINTSPACE << btot;
    f << PRINTSPACE << Btot;
    f << PRINTSPACE << stot << std::endl;
    f.close();
  }
  return -lnL;
}

void wrapper_nosrc_noBB(Int_t &nDim, Double_t *gout, Double_t &result, Double_t par[], Int_t flag){ result = nosrc_noBB(par[0]); }

void wrapper_src_noBB(Int_t &nDim, Double_t *gout, Double_t &result, Double_t par[], Int_t flag){ result = src_noBB(par[0],par[1]); }

void wrapper_nosrc_BB(Int_t &nDim, Double_t *gout, Double_t &result, Double_t par[], Int_t flag){ result = nosrc_BB(par[0]); }

void wrapper_src_BB(Int_t &nDim, Double_t *gout, Double_t &result, Double_t par[], Int_t flag){ result = src_BB(par[0], par[1]); }

void fit(indices_t ins, args_t args, double alpha){
  //Set up fitters
  TFitter* fit_nosrc_nobb = new TFitter(1);
  TFitter* fit_src_nobb   = new TFitter(2);
  TFitter* fit_nosrc_bb   = new TFitter(1);
  TFitter* fit_src_bb     = new TFitter(2);
  //Set Printout
  double p1 = args.verbosity;
  fit_nosrc_nobb->ExecuteCommand("SET PRINTOUT", &p1, 1);
  fit_src_nobb  ->ExecuteCommand("SET PRINTOUT", &p1, 1);
  fit_nosrc_bb  ->ExecuteCommand("SET PRINTOUT", &p1, 1);
  fit_src_bb    ->ExecuteCommand("SET PRINTOUT", &p1, 1);
  //Set fcn function
  fit_nosrc_nobb->SetFCN(wrapper_nosrc_noBB);
  fit_src_nobb  ->SetFCN(wrapper_src_noBB);
  fit_nosrc_bb  ->SetFCN(wrapper_nosrc_BB);
  fit_src_bb    ->SetFCN(wrapper_src_BB);
  //Set parameters
  fit_nosrc_nobb->SetParameter(0, "bkgFrac", .5, .5, 0, 1);
  fit_src_nobb  ->SetParameter(0, "bkgFrac", .5, .5, 0, 1);
  fit_nosrc_bb  ->SetParameter(0, "bkgFrac", .5, .5, 0, 1);
  fit_src_bb    ->SetParameter(0, "bkgFrac", .5, .5, 0, 1);
  fit_src_nobb  ->SetParameter(1, "srcFrac", .5, .5, 0, 1);
  fit_src_bb    ->SetParameter(1, "srcFrac", .5, .5, 0, 1);
  //Run SIMPLEX
  fit_nosrc_nobb->ExecuteCommand("SIMPLEX", 0, 0);
  fit_src_nobb  ->ExecuteCommand("SIMPLEX", 0, 0);
  fit_nosrc_bb  ->ExecuteCommand("SIMPLEX", 0, 0);
  fit_src_bb    ->ExecuteCommand("SIMPLEX", 0, 0);
  //Run MIGRAD
  fit_nosrc_nobb->ExecuteCommand("MIGRAD", 0, 0);
  fit_src_nobb  ->ExecuteCommand("MIGRAD", 0, 0);
  fit_nosrc_bb  ->ExecuteCommand("MIGRAD", 0, 0);
  fit_src_bb    ->ExecuteCommand("MIGRAD", 0, 0);

  if(args.hist & 2){
    double fracs[] = {fit_src_nobb->GetParameter(0),
                      fit_src_nobb->GetParameter(1),
                      fit_src_bb->GetParameter(0),
                      fit_src_bb->GetParameter(1)};
    histogram_fit_data(fracs , ins);
  }

  //Get likelihood and TS
  bool output_bins = (args.output & 1) && (DAT_HIST->Integral() + BKG_HIST->Integral());
  double lnL_nosrc_nobb = -nosrc_noBB(fit_nosrc_nobb->GetParameter(0), output_bins);
  double lnL_src_nobb = -src_noBB(fit_src_nobb->GetParameter(0), fit_src_nobb->GetParameter(1), output_bins);
  double lnL_nosrc_bb = -nosrc_BB(fit_nosrc_bb->GetParameter(0), output_bins);
  double lnL_src_bb = -src_BB(fit_src_bb->GetParameter(0), fit_src_bb->GetParameter(1), output_bins);
  double TS_nobb = -2 * (lnL_src_nobb - lnL_nosrc_nobb);
  double TS_bb = -2 * (lnL_src_bb - lnL_nosrc_bb);
  double lima_nobb = lima_sig(fit_src_nobb->GetParameter(0), fit_src_nobb->GetParameter(1));
  double lima_bb = lima_sig(fit_src_bb->GetParameter(0), fit_src_bb->GetParameter(1));

  std::ofstream f("fitstats_nobb.csv", std::ios::out | std::ios::app);
  if(args.bin_vars & 1) f << ins.za << ",";
  if(args.bin_vars & 2) f << ins.e << ",";
  if(args.bin_vars & 4) f << ins.tel << ",";
  if(args.bin_vars & 8) f << ins.az << ",";
  if(args.bin_vars & 16) f << ins.off << ",";
  f << std::scientific
    << fit_nosrc_nobb->GetParameter(0) << ","
    << fit_src_nobb->GetParameter(0) << ","
    << fit_src_nobb->GetParameter(1) << ","
    << std::defaultfloat
    << DAT_HIST->Integral() << ","
    << BKG_HIST->Integral() << ","
    << SRC_HIST->Integral() << ","
    << lnL_nosrc_nobb << ","
    << lnL_src_nobb << ","
    << TS_nobb << std::endl;
  f.close();

  f.open("fitstats_bb.csv", std::ios::out | std::ios::app);
  if(args.bin_vars & 1) f << ins.za << ",";
  if(args.bin_vars & 2) f << ins.e << ",";
  if(args.bin_vars & 4) f << ins.tel << ",";
  if(args.bin_vars & 8) f << ins.az << ",";
  if(args.bin_vars & 16) f << ins.off << ",";
  f << std::scientific
    << fit_nosrc_bb->GetParameter(0) << ","
    << fit_src_bb->GetParameter(0) << ","
    << fit_src_bb->GetParameter(1) << ","
    << std::defaultfloat
    << DAT_HIST->Integral() << ","
    << BKG_HIST->Integral() << ","
    << SRC_HIST->Integral() << ","
    << lnL_nosrc_bb << ","
    << lnL_src_bb << ","
    << TS_bb << std::endl;
  f.close();

  f.open("summary.csv", std::ios::out | std::ios::app);
  if(args.bin_vars & 1) f << ins.za << ",";
  if(args.bin_vars & 2) f << ins.e << ",";
  if(args.bin_vars & 4) f << ins.tel << ",";
  if(args.bin_vars & 8) f << ins.az << ",";
  if(args.bin_vars & 16) f << ins.off << ",";
  f << BKG_HIST->Integral() / DAT_HIST->Integral() << ","
    << TS_nobb << ","
    << TS_bb << ","
    << lima_nobb << ","
    << lima_bb << std::endl;
  f.close();

  if(args.output & 4){
    calculate_errors(fit_src_bb->GetParameter(0),
                    fit_src_bb->GetParameter(1),
                    fit_src_bb->GetParError(0),
                    fit_src_bb->GetParError(1),
                    ins, alpha);
  }
  if(args.graphics & 1) map_likelihood(fit_src_nobb->GetParameter(0), fit_src_nobb->GetParameter(1), "Std", ins, args);
  if(args.graphics & 2) map_likelihood(fit_src_bb->GetParameter(0), fit_src_bb->GetParameter(1), "BB", ins, args);
}

int main(int argc, char* argv[]){
  static void (*loaddata[])(indices_t ins, args_t args, double* alpha) = {loadData, loadData_csv, loadData_vegas};
  args_t* args = new args_t;
  if(parse_command_line(argc, argv, args)) return 0;

  //Setup output files
  std::ofstream f1("fitstats_bb.csv");
  std::ofstream f2("fitstats_nobb.csv");
  std::ofstream f3("summary.csv");
  if(args->bin_vars & 1)  {f1 << "ZA,"; f2 << "ZA,"; f3 << "ZA,";}
  if(args->bin_vars & 2)  {f1 << "E,"; f2 << "E,"; f3 << "E,";}
  if(args->bin_vars & 4)  {f1 << "T,"; f2 << "T,"; f3 << "T,";}
  if(args->bin_vars & 8)  {f1 << "A,"; f2 << "A,"; f3 << "A,";}
  if(args->bin_vars & 16) {f1 << "O,"; f2 << "O,"; f3 << "O,";}

  f1 << "bkgfrac_nosrc, bkgfrac, srcfrac, dataCt, bkgCt, srcCt, lnL_nosrc, lnL_src, TS" << std::endl;
  f2 << "bkgfrac_nosrc, bkgfrac, srcfrac, dataCt, bkgCt, srcCt, lnL_nosrc, lnL_src, TS" << std::endl;
  f3 << "ct_ratio,TS_noBB,TS_BB,Lima_std,Lima_bb" << std::endl;
  f1.close(); f2.close(); f3.close();

  if(args->output & 8) print_cuts("reset", 0);
  int zi = 1, ei = 0, ti = 0, ai = 0, oi = 0;
  indices_t indices;
  for(indices.za = zi; indices.za < 2; indices.za++){
    for(indices.e = ei; indices.e < 4; indices.e++){
      for(indices.tel = ti; indices.tel < 2; indices.tel++){
        for(indices.az = ai; indices.az < 8; indices.az++){
          for(indices.off = oi; indices.off < 8; indices.off++){
            //This takes care of optional binning in a rather crude manner
            std::stringstream path;
            if(!(args->bin_vars & 1)){
              if(indices.za != zi) continue;
            }
            else path << "ZA" << indices.za;
            if(!(args->bin_vars & 2)){
              if(indices.e != ei) continue;
            }
            else path << "E" << indices.e;
            if(!(args->bin_vars & 4)){
              if(indices.tel != ti) continue;
            }
            else path << "T" << TBINS[indices.tel];
            if(!(args->bin_vars & 8)){
              if(indices.az != ai) continue;
            }
            else path << "A" << indices.az;
            if(!(args->bin_vars & 16)){
              if(indices.off != oi) continue;
            }
            else path << "O" << indices.off;

            //Run Fit
            std::cout << path.str() << std::endl;
            OUTPATH = path.str();
            double alpha = 1;
            loaddata[((int)args->format)](indices, *args, &alpha);
            if(!DAT_HIST || !BKG_HIST || !SRC_HIST) throw 407;
            fit(indices, *args, alpha);
            if(args->output & 2) printRawData();
            if(args->hist & 1) histogram_raw_data(indices);
            delete DAT_HIST;
            delete BKG_HIST;
            delete SRC_HIST;
          }
        }
      }
    }
  }
  delete args;
}

//Extras
int parse_command_line(int argc, char* argv[], args_t* args){
  for(int i = 0; i < argc; i++){
    if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
      std::cout << "OPTIONS:" << std::endl
                << "  -d FORMAT, --data-format FORMAT" << std::endl
                << "    Format of the imput data. Available: toy, csv, vegas. Default: toy." << std::endl
                << "  -h, --help" << std::endl
                << "    Print this message." << std::endl
                << "  -hist DATA, --histogram DATA" << std::endl
                << "    Triggers output of histograms. Available: none, raw, fit, all. Default: none." << std::endl
                << "  -out DATA, --output DATA" << std::endl
                << "    Triggers output of ascii files. Available: none, bins, raw, errors, cuts, all. Default: none." << std::endl
                << "  -v VERBOSITY, --verbosity VERBOSITY" << std::endl
                << "    Controls the verbosity of MINUIT. Default: -1." << std::endl
                << "  -b VAR, --bin-variable VAR" << std::endl
                << "    Sets which variables to bin over. Available: zenith, energy, telescope, azimuth, offset, all. Default: zenith, energy, telescope." << std::endl
                << "  -g GRAPHICS, --graphics GRAPHICS" << std::endl
                << "    Triggers output of graphics files. Available: none, stdlnL, bblnL, all. Default: none." << std::endl;
      return 1;
    }
    if(!strcmp(argv[i], "-d") || !strcmp(argv[i], "--data-format")){
      if(i < argc - 1 && !strcmp(argv[i+1], "toy")){
        args->format = Format_t::Toy;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "csv")){
        args->format = Format_t::Csv;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "vegas")){
        args->format = Format_t::Vegas;
      }
    }
    if(!strcmp(argv[i], "-hist") || !strcmp(argv[i], "--histogram")){
      if(i < argc - 1 && !strcmp(argv[i+1], "none")){
        args->hist = 0;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "raw")){
        args->hist += 1;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "fit")){
        args->hist += 2;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "all")){
        args->hist = 3;
      }
    }
    if(!strcmp(argv[i], "-out") || !strcmp(argv[i], "--output")){
      if(i < argc - 1 && !strcmp(argv[i+1], "none")){
        args->output = 0;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "bins")){
        args->output += 1;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "raw")){
        args->output += 2;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "errors")){
        args->output += 4;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "cuts")){
        args->output += 8;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "all")){
        args->output = 15;
      }
    }
    if(!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbosity")){
      if(i < argc - 1){
        args->verbosity = atoi(argv[++i]);
      }
    }
    if(!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bin-variable")){
      if(i < argc - 1 && !strcmp(argv[i+1], "zenith") && !(args->bin_vars & 1)){
        args->bin_vars += 1;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "energy") && !(args->bin_vars & 2)){
        args->bin_vars += 2;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "telescope") && !(args->bin_vars & 4)){
        args->bin_vars += 4;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "azimuth") && !(args->bin_vars & 8)){
        args->bin_vars += 8;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "offset") && !(args->bin_vars & 16)){
        args->bin_vars += 16;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "all")){
        args->bin_vars = 31;
      }
    }
    if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--graphics")){
      if(i < argc - 1 && !strcmp(argv[i+1], "none")){
        args->graphics = 0;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "stdlnL") && !(args->graphics & 1)){
        args->graphics += 1;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "bblnL") && !(args->graphics & 2)){
        args->graphics += 2;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "all")){
        args->graphics = 3;
      }
    }
  }
  return 0;
}

void printRawData(){
  if(!(DAT_HIST->Integral() + BKG_HIST->Integral())) return;
  std::stringstream path;
  path << "raw_data_" << OUTPATH << ".csv";
  std::ofstream f(path.str().c_str());
  f << "MSW, di, bi, si" << std::endl;
  for(int i = 1; i <= NBIN; i++){
    f << DAT_HIST->GetBinCenter(i) << ","
      << DAT_HIST->GetBinContent(i) << ","
      << BKG_HIST->GetBinContent(i) << ","
      << SRC_HIST->GetBinContent(i) << std::endl;
  }
}

void histogram_raw_data(indices_t ins){
  if(!(DAT_HIST->Integral() + BKG_HIST->Integral())) return;
  std::stringstream filepath;
  filepath << "HIST_RAW_" << OUTPATH << ".png";

  //Copy the global histograms so that we can edit them
  TH1F* dathist = new TH1F(*DAT_HIST);
  TH1F* bkghist = new TH1F(*BKG_HIST);
  TH1F* srchist = new TH1F(*SRC_HIST);

  //Set up
  TCanvas *c1 = new TCanvas("","",1200,1200);
  c1->cd();
  dathist->SetLineColor(4);
  bkghist->SetLineColor(6);
  srchist->SetLineColor(1);

  dathist->SetStats(false);
  bkghist->SetStats(false);
  srchist->SetStats(false);

  dathist->SetTitle(OUTPATH.c_str());

  bkghist->Scale(DAT_HIST->Integral() / BKG_HIST->Integral());
  srchist->Scale(DAT_HIST->Integral() / SRC_HIST->Integral());

  TLegend *legend;
  legend = new TLegend(0.12, 0.8, 0.3, 0.9);
  legend->AddEntry(dathist, "Data");
  legend->AddEntry(bkghist, "Background Template");
  legend->AddEntry(srchist, "Source Template");

  dathist->SetMinimum(0);
  dathist->SetMaximum(std::max(dathist->GetMaximum(),
                      std::max(bkghist->GetMaximum(),
                               srchist->GetMaximum())) * 1.1);

  //Draw
  dathist->Draw();
  dathist->Draw("sameE0");
  bkghist->Draw("same");
  bkghist->Draw("sameE0");
  srchist->Draw("same hist");
  legend->Draw();

  //Save and clean up
  c1->SaveAs(filepath.str().c_str());
  c1->Clear();
  delete c1;
  delete dathist;
  delete bkghist;
  delete srchist;
  delete legend;
}

void histogram_fit_data(double fracs[4], indices_t ins){
  if(!(DAT_HIST->Integral() + BKG_HIST->Integral())) return;
  std::stringstream filepath;
  filepath << "HIST_FIT_" << OUTPATH << ".png";

  //Copy the global histograms so that we can edit them
  TH1F* dathist = new TH1F(*DAT_HIST);
  TH1F* bkghist = new TH1F(*BKG_HIST);
  TH1F* srchist = new TH1F(*SRC_HIST);
  bkghist->Scale(DAT_HIST->Integral() / BKG_HIST->Integral());
  srchist->Scale(DAT_HIST->Integral() / SRC_HIST->Integral());

  //Canvas set up
  TH1F* F0 = new TH1F("F0Hist", "Std Fit", NBIN, MSWLOW, MSWHIGH);
  TH1F* F1 = new TH1F("F1Hist", "BB Fit", NBIN, MSWLOW, MSWHIGH);
  TH1F* B1 = new TH1F("B1Hist", "BB Fit", NBIN, MSWLOW, MSWHIGH);
  TCanvas *c1 = new TCanvas("",OUTPATH.c_str(),1600,1600);
  c1->Divide(2,2);

  dathist->SetLineColor(4);
  bkghist->SetLineColor(6);
  srchist->SetLineColor(1);
  F0->SetLineColor(8);
  F1->SetLineColor(8);
  B1->SetLineColor(6);

  dathist->SetStats(false);
  bkghist->SetStats(false);
  srchist->SetStats(false);
  F0->SetStats(false);
  F1->SetStats(false);
  B1->SetStats(false);

  TLegend *legend0 = new TLegend(0.15, 0.8, 0.25, 0.9);
  TLegend *legend1 = new TLegend(0.15, 0.8, 0.25, 0.9);

  //Draw Std Src Fit
  bkghist->Scale(fracs[0]);
  srchist->Scale(fracs[1]);
  src_noBB(fracs[0], fracs[1], false, F0);
  c1->cd(1);
  F0->Draw("hist");
  dathist->Draw("same");
  bkghist->Draw("same hist");
  srchist->Draw("same hist");

  F0->SetMinimum(0);
  F0->SetMaximum(std::max(
    std::max(dathist->GetMaximum(), F0->GetMaximum()),
    std::max(srchist->GetMaximum(), bkghist->GetMaximum())
  ) * 1.1);

  legend0->AddEntry(F0, "fi");
  legend0->AddEntry(bkghist, "bi");
  legend0->AddEntry(dathist, "di");
  legend0->AddEntry(srchist, "Si");
  legend0->Draw();

  srchist->Scale(1/fracs[1]);

  //Draw BB Src
  srchist->Scale(fracs[3]);
  src_BB(fracs[2], fracs[3], false, F1, B1);
  B1->Scale(DAT_HIST->Integral() / BKG_HIST->Integral());
  B1->Scale(fracs[2]);
  c1->cd(2);
  F1->Draw("hist");
  dathist->Draw("same");
  B1->Draw("same hist");
  srchist->Draw("same hist");

  F1->SetMinimum(0);
  F1->SetMaximum(std::max(
    std::max(dathist->GetMaximum(), F1->GetMaximum()),
    std::max(srchist->GetMaximum(), B1->GetMaximum())
  ) * 1.1);

  legend1->AddEntry(F1, "fi");
  legend1->AddEntry(B1, "Bi");
  legend1->AddEntry(dathist, "di");
  legend1->AddEntry(srchist, "Si");
  legend1->Draw();

  //Draw Residual
  c1->cd(3);
  TH1F* res = new TH1F("Residual", "Residual", NBIN, MSWLOW, MSWHIGH);
  res->Add(F0, F1, 1, -1);
  res->Draw("");

  //Save and Clean up
  c1->SaveAs(filepath.str().c_str());
  c1->Clear();
  delete c1;
  delete dathist;
  delete bkghist;
  delete srchist;
  delete F0;
  delete F1;
  delete B1;
  delete res;
  delete legend0;
  delete legend1;
}

void calculate_errors(double Pb_m, double Ps_m, double sigma_Pb_m, double sigma_Ps_m, indices_t ins, double alpha){
  if(!(DAT_HIST->Integral() + BKG_HIST->Integral())) return;
  double N_D = DAT_HIST->Integral();
  double N_B = BKG_HIST->Integral();
  double B = alpha * N_B;
  double S = N_D - B;
  double P_S = S / N_D;
  double P_B = B / N_D;

  double sigma_N_D = TMath::Sqrt(DAT_HIST->Integral());
  double sigma_B = TMath::Sqrt(alpha) * TMath::Sqrt(B);
  double sigma_S = TMath::Sqrt( TMath::Power(sigma_N_D, 2) + TMath::Power(sigma_B, 2) );
  double sigma_P_S = TMath::Abs(P_S) * TMath::Sqrt( TMath::Power(sigma_S / S, 2) + TMath::Power(sigma_N_D / N_D, 2) );
  double sigma_P_B = TMath::Abs(P_B) * TMath::Sqrt( TMath::Power(sigma_B / B, 2) + TMath::Power(sigma_N_D / N_D, 2) );

  std::stringstream filename;
  filename << "errors" << OUTPATH << ".csv";
  std::ofstream f(filename.str().c_str());
  f << "Minuit Values:" << std::endl
    << "Ps: " << Ps_m << std::endl
    << "sigma_Ps: " << sigma_Ps_m << std::endl
    << "Pb: " << Pb_m << std::endl
    << "sigma_Pb: " << sigma_Pb_m << std::endl;
  f << std::endl;
  f << "Calculated Values:" << std::endl
    << "Ps: " << P_S << std::endl
    << "sigma_Ps: " << sigma_P_S << std::endl
    << "Pb: " << P_B << std::endl
    << "sigma_Pb: " << sigma_P_B << std::endl;
  f << std::endl;
  f << "Relative Differences:" << std::endl
    << "Ps: " << TMath::Abs(P_S - Ps_m) / Ps_m << std::endl
    << "sigma_Ps: " << TMath::Abs(sigma_Ps_m - sigma_P_S) / sigma_Ps_m << std::endl
    << "Pb: " << TMath::Abs(P_B - Pb_m) / Pb_m << std::endl
    << "sigma_Pb: " << TMath::Abs(sigma_Pb_m - sigma_P_B) / sigma_Pb_m << std::endl;
  f << std::endl;
  f << "N_D: " << N_D << std::endl
    << "sigma_N_D: " << sigma_N_D << std::endl
    << "B: " << B << std::endl
    << "sigma_B: " << sigma_B << std::endl
    << "S: " << S << std::endl
    << "sigma_S: " << sigma_S << std::endl
    << "alpha: " << alpha << std::endl;

  f.close();
}

void print_cuts(std::string action, cuts_t* cuts){
  if(!action.compare("reset")){
    std::ofstream f("cuts_data.txt");
    f << PRINTSPACE << "Bin"
      << PRINTSPACE << "passed"
      << PRINTSPACE << "read"
      << PRINTSPACE << "source"
      << PRINTSPACE << "telescope"
      << PRINTSPACE << "energy"
      << PRINTSPACE << "zenith angle"
      << PRINTSPACE << "msw"
      << PRINTSPACE << "azimuth"
      << PRINTSPACE << "offset" << std::endl;
    f.close();
    f.open("cuts_bkg.txt");
    f << PRINTSPACE << "Bin"
      << PRINTSPACE << "passed"
      << PRINTSPACE << "read"
      << PRINTSPACE << "source"
      << PRINTSPACE << "telescope"
      << PRINTSPACE << "energy"
      << PRINTSPACE << "zenith angle"
      << PRINTSPACE << "msw"
      << PRINTSPACE << "azimuth"
      << PRINTSPACE << "offset" << std::endl;
    f.close();
    f.open("cuts_src.txt");
    f << PRINTSPACE << "Bin"
      << PRINTSPACE << "Events" << std::endl;
    f.close();
  }
  else if(!action.compare("data") || !action.compare("bkg")){
    std::stringstream fname;
    fname << "cuts_" << action << ".txt";
    std::ofstream f(fname.str().c_str(), std::ios::out | std::ios::app);
    f << PRINTSPACE << OUTPATH
      << PRINTSPACE << cuts->passed
      << PRINTSPACE << cuts->read
      << PRINTSPACE << cuts->src
      << PRINTSPACE << cuts->tel
      << PRINTSPACE << cuts->e
      << PRINTSPACE << cuts->za
      << PRINTSPACE << cuts->msw
      << PRINTSPACE << cuts->az
      << PRINTSPACE << cuts->off << std::endl;
    f.close();
  }
  else if(!action.compare("src")){
    std::ofstream f("cuts_src.txt", std::ios::out | std::ios::app);
    f << PRINTSPACE << OUTPATH
      << PRINTSPACE << SRC_HIST->Integral() << std::endl;
    f.close();
  }
}

void fit_manual_bin(){
  //load data from input_data.csv
  //data must be in a csv file of form di, bi, si
  TH1::SetDefaultSumw2();
  DAT_HIST = new TH1F("DataHist", "Data", NBIN, MSWLOW, MSWHIGH);
  BKG_HIST =  new TH1F("BkgHist",  "BKG", NBIN, MSWLOW, MSWHIGH);
  SRC_HIST =  new TH1F("SrcHist" , "SRC" ,NBIN, MSWLOW, MSWHIGH);

  std::string line;
  std::vector<std::string> line_fields;
  std::ifstream input("input_data.csv");
  for(int i = 0; i < NBIN; i++){
    std::getline(input, line);
    boost::split(line_fields, line, boost::is_any_of(","));
    DAT_HIST->Sumw2(false);
    DAT_HIST->SetBinContent(i+1, std::atof(line_fields.at(0).c_str()));
    DAT_HIST->Sumw2(true);
    BKG_HIST->Sumw2(false);
    BKG_HIST->SetBinContent(i+1, std::atof(line_fields.at(1).c_str()));
    BKG_HIST->Sumw2(true);
    SRC_HIST->Sumw2(false);
    SRC_HIST->SetBinContent(i+1, std::atof(line_fields.at(2).c_str()));
    SRC_HIST->Sumw2(true);
    line_fields.clear();
  }

  //Make indices and args structs and OUTPATH
  //these are manually set; indices and OUTPATH are cosmetic
  args_t args;
  args.hist = 3;
  args.output = 11;
  args.bin_vars = 31;
  indices_t ins;
  ins.za = 1;
  ins.e = 3;
  ins.tel = 0;
  ins.az = 0;
  ins.off = 2;
  OUTPATH = "Match_only78 ZA1E2T3A0O5";

  //perform fit
  fit(ins, args, 1);
  printRawData();
  histogram_raw_data(ins);

  //clean up
  delete DAT_HIST;
  delete BKG_HIST;
  delete SRC_HIST;


}

void map_likelihood(double Pb, double Ps, std::string title_tag, indices_t ins, args_t args){
  if(!(DAT_HIST->Integral() + BKG_HIST->Integral())) return;
  int xbins = 100;
  int ybins = 100;
  //Find x limits such that xmax <= 1, xmin >= 0, and xmax - xmin = .3
  double xmax = (Pb + .15 > 1.0 ? 1.0 : Pb + .15);
  double xmin = Pb - .15 < 0.0 ? 0.0 : Pb - .15;
  if(xmax - xmin < .3){
    if(xmax == 1) xmin = .7;
    else if(xmin == 0) xmax = .3;
  }
  //Find y limits similarly
  double ymax = Ps + .15 > 1.0 ? 1.0 : Ps + .15;
  double ymin = Ps - .15 < 0.0 ? 0.0 : Ps - .15;
  if(ymax - ymin < .3){
    if(ymax == 1) ymin = .7;
    else if(ymin == 0) ymax = .3;
  }

  std::stringstream title;
  title << title_tag
        << " Likelihood ";
  if(args.bin_vars & 1) title << "ZA"<< ins.za;
  if(args.bin_vars & 2) title << "E" << ins.e;
  if(args.bin_vars & 4) title << "T" << ins.tel;
  if(args.bin_vars & 8) title << "A" << ins.az;
  if(args.bin_vars & 16) title << "O" << ins.off;

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1", 1200,800);
  TH2F *map = new TH2F("map", title.str().c_str(), xbins, xmin, xmax, ybins, ymin, ymax);
  TAxis *xax = map->GetXaxis();
  TAxis *yax = map->GetYaxis();
  xax->SetTitle("Pb");
  yax->SetTitle("Ps");
  for(int j = 1; j <= ybins; j++){
    for(int i = 1; i <= xbins; i++){
      double Pb = xax->GetBinCenter(i);
      double Ps = yax->GetBinCenter(j);
      double lnL = -src_BB(Pb,Ps);
      map->SetBinContent(i, j, lnL);
    }
  }

  gStyle->SetPalette(kSunset);
  map->Draw("SURF3");
  std::stringstream out_file;
  out_file << boost::replace_all_copy(title.str(), " ", "_") << ".png";
  c1->SaveAs(out_file.str().c_str());

  //clean up
  delete c1;
  delete map;
}
