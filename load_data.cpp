#include "load_data.h"

std::string data_excl_file = "data_src.txt";
std::string bkg_excl_file = "bkg_src.txt";

void loadsrc_csv(indices_t ins, args_t args){
  std::cout << "Loading Source" << std::endl;
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

  if(args.output & 8) print_cuts("src", 0);
}

void loadData_toy(indices_t ins, args_t args, double* alpha){
  std::string line;
  TFile *f;
  std::ifstream data_flist("data.list");
  std::ifstream bkg_flist("bkg.list");

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
  if(args.output & 8){
    print_cuts("data", &data_cuts);
    print_cuts("bkg", &bkg_cuts);
  }
  std::cout << "Data Overflow: " << data_overflow << std::endl;
  std::cout << "Bkg Overflow: " << bkg_overflow << std::endl;
  *alpha = livetime_data / livetime_bkg;
}

/*
TH1F** loadData_csv(indices_t ins, args_t args, double* alpha){
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
  TH1F* DAT_HIST = new TH1F("DataHist", "Data", NBIN, MSWLOW, MSWHIGH);
  TH1F* BKG_HIST =  new TH1F("BkgHist" , "BKG" , NBIN, MSWLOW, MSWHIGH);
  TH1F* SRC_HIST =  new TH1F("SrcHist" , "SRC" , NBIN, MSWLOW, MSWHIGH);

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

  return {DAT_HIST, BKG_HIST, SRC_HIST};
}
*/

void loadData_vegas(indices_t ins, args_t args, double* alpha){
  //Check for source exclusion files
  bool data_excl_exists = true;
  bool bkg_excl_exists = true;
  if(access(data_excl_file.c_str(), F_OK)){
    std::cerr << "No source exclusion file found at "
    << data_excl_file
    << ". This may cause problems."
    << std::endl;
    data_excl_exists = false;
  }
  if(access(bkg_excl_file.c_str(), F_OK)){
    std::cerr << "No source exclusion file found at "
    << bkg_excl_file
    << ". This may cause problems."
    << std::endl;
    bkg_excl_exists = false;
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
    if(fail_src_cuts && data_excl_exists){
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
    if(fail_src_cuts && bkg_excl_exists){
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

  if(args.output & 8){
    print_cuts("data", &data_cuts);
    print_cuts("bkg", &bkg_cuts);
  }

  std::cout << "Data Offset Overflow: " << data_overflow << std::endl;
  std::cout << "Bkg Offset Overflow: " << bkg_overflow << std::endl;

  //TODO
  *alpha = 1.0;
}

void loadData_sample(indices_t ins, args_t args, double* alpha){
  /*Load sample data that is easily and quickly repeatable for testing purposes.
   * Data hist is filled from the function 5x with ~1000 counts.
   * Bkg hist is filled from the function 5x with ~3000 counts.
   * Src hist is filled from a Gaussian with ~30000 counts.
   */

  for(int i = 1; i <= NBIN; i++){
    double dat = 5*i + 1000/NBIN + gRandom->Gaus(0, 10);
    double bkg = 5*i + 3000/NBIN + gRandom->Gaus(0, 10);
    DAT_HIST->SetBinContent(i, dat);
    BKG_HIST->SetBinContent(i, bkg);
  }
  for(int i = 0; i < 30000; i++){
    SRC_HIST->Fill(gRandom->Gaus(1, .1));
  }
  *alpha = 3/5;
}

void loadData(indices_t ins, args_t args, double *alpha){
  TH1::SetDefaultSumw2();
  DAT_HIST = new TH1F("DataHist", "Data", NBIN, MSWLOW, MSWHIGH);
  BKG_HIST = new TH1F("BkgHist", "BKG", NBIN, MSWLOW, MSWHIGH);
  SRC_HIST = new TH1F("SrcHist", "SRC", NBIN, MSWLOW, MSWHIGH);


  if(args.format == Format_t::Toy){
    loadData_toy(ins, args, alpha);
    loadsrc_csv(ins, args);
    std::cout << "Histograms loaded from Toy format." << std::endl;
  }
  else if(args.format == Format_t::Vegas){
    loadData_vegas(ins, args, alpha);
    loadsrc_csv(ins, args);
    std::cout << "Histograms loaded from Vegas format." << std::endl;
  }
  else if(args.format == Format_t::Sample){
    loadData_sample(ins, args, alpha);
    std::cout << "Histograms loaded from Sample format." << std::endl;
  }
  else{
    std::cerr << "No valid data format specified." << std::endl;
    throw 999;
  }
}