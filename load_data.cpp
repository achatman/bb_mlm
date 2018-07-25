#include "load_data.h"
#include "output.h"

std::string OUTSTR;

std::vector<std::pair<double,double>> bkg_centers;

void loadsrc_csv(indices_t ins, args_t args, TH1D* SRC_HIST){
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
}

double loadData_toy(indices_t ins, args_t args, std::string pathbase, TH1D* HIST){
  std::string line;
  TFile *f;
  std::stringstream sstream;
  sstream << pathbase << ".list";
  std::ifstream flist(sstream.str().c_str());

  //Check each point against cuts
  cuts_t cuts;
  double livetime = 0;
  int overflow = 0;

  //Construct Data Histogram
  std::cout << "Loading " << pathbase << std::endl;
  timestamp;
  while(std::getline(flist, line)){
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
      cuts.read++;
      //Telescope Cut
      if(args.bin_vars & 4){
        if(*telescopes != TBINS[ins.tel]){
          cuts.tel++;
          continue;
        }
      }

      //Energy Cut
      if(args.bin_vars & 2){
        double en = *energy * 1000;
        if(en > EBINS[ins.e+1] || en < EBINS[ins.e]){
          cuts.e++;
          continue;
        }
      }

      //ZA Cut
      if(args.bin_vars & 1){
        if(90-(*za) > ZABINS[ins.za + 1] || 90-(*za) < ZABINS[ins.za]){
          cuts.za++;
          continue;
        }
      }

      //MSW Cut
      if(*msw < MSWLOW || *msw > MSWHIGH){
        cuts.msw++;
        continue;
      }

      //AZ Cut
      if(args.bin_vars & 8){
        if(*az > AZBINS[ins.az + 1] || *az < AZBINS[ins.az]){
          cuts.az++;
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
        if(offset > 2.0) overflow++;
        if(offset > OBINS[ins.off + 1] || offset < OBINS[ins.off]){
          cuts.off++;
          continue;
        }
      }

      HIST->Fill(*msw);
      livetime += *live;
    }
    f->Close();
  }
  cuts.passed = HIST->Integral();
  std::cout << cuts.passed << " passed cuts." << std::endl;
  std::cout << cuts.read << " points read." << std::endl;
  std::cout << cuts.tel << " failed tel cut." << std::endl;
  std::cout << cuts.e << " failed energy cut." << std::endl;
  std::cout << cuts.za << " failed za cut." << std::endl;
  std::cout << cuts.msw << " failed msw cut." << std::endl;
  std::cout << cuts.az << " failed az cut." << std::endl;
  std::cout << cuts.off << " failed off cut." << std::endl;


  if(args.output & 8){
    print_cuts(pathbase, &cuts, OUTSTR);
  }
  std::cout << pathbase << " Overflow: " << overflow << std::endl;
  return livetime;
}

void cacheData_vegas(args_t args, std::string pathbase){
  std::stringstream excl_file;
  excl_file << pathbase << "_src.txt";
  //Check for source exclusion file
  bool excl_exists = true;
  if(access(excl_file.str().c_str(), F_OK)){
    std::cerr << "No source exclusion file found at "
    << excl_file.str()
    << ". This may cause problems."
    << std::endl;
    excl_exists = false;
  }
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

  std::vector<SourceCut_t> sourcecuts;
  std::string line;
  std::vector<std::string> line_fields;
  std::ifstream srclist(excl_file.str().c_str());
  while(std::getline(srclist, line)){
    boost::split(line_fields, line, boost::is_any_of(","));
    sourcecuts.emplace_back(line_fields);
    line_fields.clear();
  }

  //Setup flist
  std::stringstream infile;
  infile << pathbase << ".list";
  std::ifstream flist(infile.str().c_str());

  //Set up RA/DEC Hist
  TH2D* ra_dec = new TH2D("RA_DEC", "RA_DEC", 360, 0, 360, 360, 0, 360);

  //Set up Cache file
  std::stringstream cache_path;
  cache_path << "Cache_" << pathbase << ".csv";
  std::ofstream cache_file(cache_path.str().c_str());

  cuts_t cuts;
  timestamp;
  std::cout << "Caching " << pathbase << std::endl;
  //Cache events
  int bin_counts[6][4][2];
  for(int a = 0; a < 6; a++){
    for(int b = 0; b < 4; b++){
      for(int c = 0; c < 2; c++){
        bin_counts[a][b][c] = 0;
      }
    }
  }

  std::ofstream temp_file;
  if(pathbase == "data"){
    temp_file.open("data_check.csv");
  }

  TFile *f;
  while(std::getline(flist, line)){
    f = TFile::Open(line.c_str());
    if(!f){
      std::cout << "Cannot open " << line << std::endl;
      continue;
    }

    TTreeReader reader("SelectedEvents/CombinedEventsTree", f);
    TTreeReaderValue<VAShowerData> shower(reader, "S");

    while(reader.Next()){
      bool cont = false;
      cuts.read++;
      //Source Cut
      Double_t eventRA = shower->fDirectionRA_J2000_Rad * TMath::RadToDeg();
      Double_t eventDec = shower->fDirectionDec_J2000_Rad * TMath::RadToDeg();
      VACoordinatePair eventCoord = VACoordinatePair(eventRA,eventDec,VACoordinates::J2000,VACoordinates::Deg);
      bool fail_src_cuts;
      for(auto &it: sourcecuts){
        fail_src_cuts = it.InsideExclRadius(eventCoord);
        if(fail_src_cuts) break;
      }
      if(fail_src_cuts && excl_exists){
        cuts.src++;
        cont = true;
      }

      //Assign Cut Variables
      double msw = shower->fMSW;
      double msl = shower->fMSL;
      auto tels_used = shower->fTelUsedInReconstruction;
      int tels = count(tels_used.begin(), tels_used.end(), 1);
      double energy = shower->fEnergy_GeV;
      double za = 90.0 - (shower->fDirectionElevation_Rad * TMath::RadToDeg());
      double az = shower->fDirectionAzimuth_Rad * TMath::RadToDeg();
      VACoordinatePair shower_coords = VACoordinatePair(
        shower->fDirectionRA_J2000_Rad,
        shower->fDirectionDec_J2000_Rad,
        VACoordinates::J2000,
        VACoordinates::Rad
      );
      VACoordinatePair tracking_coords = VACoordinatePair(
        shower->fArrayTrackingRA_J2000_Rad,
        shower->fArrayTrackingDec_J2000_Rad,
        VACoordinates::J2000,
        VACoordinates::Rad
      );
      double offset = tracking_coords.angularSeparation_Deg(shower_coords);

      //MSW cut
      if(msw < MSWLOW || msw > MSWHIGH){
        cuts.msw++;
        cont = true;
      }

      //Tel cut
      if(tels != 3 && tels != 4){
        cuts.tel++;
        cont = true;
      }

      //Energy cut
      int e_bin = -1;
      for(int e = 0; e < 4; e++){
        if(energy < EBINS[e+1] && energy > EBINS[e]){
          e_bin = e;
        }
      }
      if(e_bin == -1){
        cuts.e++;
        cont = true;
      }

      //ZA cut
      int z_bin = -1;
      for(int z = 0; z < 6; z++){
        if(za < ZABINS[z+1] && za > ZABINS[z]){
          z_bin = z;
        }
      }
      if(z_bin == -1){
        cuts.za++;
        cont = true;
      }

      //AZ cut
      int a_bin = -1;
      for(int a = 0; a < 8; a++){
        if(az < AZBINS[a+1] && az > AZBINS[a]){
          a_bin = a;
        }
      }
      if(a_bin == -1){
        cuts.az++;
        cont = true;
      }

      //Offset cut
      int o_bin = -1;
      for(int o = 0; o < 8; o++){
        if(offset < OBINS[o+1] && offset > OBINS[o]){
          o_bin = o;
        }
      }
      if(o_bin == -1){
        cuts.off++;
        cont = true;
      }

      if(cont) continue;
      ra_dec->Fill(eventRA, eventDec);
      bin_counts[z_bin][e_bin][o_bin]++;

      cache_file << msw << ","
                << msl << ","
                << z_bin << ","
                << e_bin << ","
                << tels-3 << ","
                << a_bin << ","
                << o_bin << ","
                << eventRA << ","
                << eventDec << std::endl;

      temp_file << msw << ","
                << msl << ","
                << za << ","
                << energy << ","
                << tels << ","
                << az << ","
                << offset << ","
                << eventRA << ","
                << eventDec << std::endl;
    }
  }
  std::cout << cuts.passed << " passed cuts." << std::endl;
  std::cout << cuts.read << " cuts read." << std::endl;
  std::cout << cuts.src << " failed src cut." << std::endl;
  std::cout << cuts.tel << " failed tel cut." << std::endl;
  std::cout << cuts.e << " failed energy cut." << std::endl;
  std::cout << cuts.za << " failed za cut." << std::endl;
  std::cout << cuts.msw << " failed msw cut." << std::endl;
  std::cout << cuts.az << " failed az cut." << std::endl;
  std::cout << cuts.off << " failed off cut." << std::endl;

  for(int a = 0; a < 6; a++){
    for(int b = 0; b < 4; b++){
      for(int c = 0; c < 2; c++){
        std::cout << a << b << c << " " << bin_counts[a][b][c] << std::endl;
      }
    }
  }

  if(args.output & 8) print_cuts(pathbase, &cuts, "Cache");
  if(pathbase != "data"){
    bkg_centers.push_back(std::make_pair(ra_dec->GetMean(1), ra_dec->GetMean(2)));
  }
  cache_file.close();
  temp_file.close();
}

void loadData_sample(indices_t ins, args_t args, std::string pathbase, TH1D* HIST, TH2D* HIST2D = 0){
  /*Load sample data that is easily and quickly repeatable for testing purposes.
   * Data hist is filled from the function 5x with ~1000 counts.
   * Bkg hist is filled from the function 5x with ~3000 counts.
   * Src hist is filled from a Gaussian with ~30000 counts.
   */

  if(pathbase == "data"){
    for(int i = 1; i <= NBIN; i++){
      double dat = 5*i + 1000/NBIN + gRandom->Gaus(0, 10);
      HIST->SetBinContent(i, dat);
    }
  }
  if(pathbase == "bkg"){
    for(int i = 1; i <= NBIN; i++){
      double bkg = 5*i + 3000/NBIN + gRandom->Gaus(0, 10);
      HIST->SetBinContent(i, bkg);
    }
  }
  if(pathbase == "src"){
    for(int i = 0; i < 30000; i++){
      HIST->Fill(gRandom->Gaus(1, .1));
    }
  }
  if(HIST2D){
    for(int i = 1; i <= NBIN; i++){
      double nEvents = HIST->GetBinContent(i);
      for(int j = 0; j < nEvents; j++){
        double x = gRandom->Uniform(HIST->GetBinCenter(i) - .5*HIST->GetBinWidth(i), HIST->GetBinCenter(i) + .5*HIST->GetBinWidth(i));
        double y = gRandom->Gaus(1, .15);
        HIST2D->Fill(x,y);
      }
    }
  }
}

double loadData_vegas(indices_t ins, args_t args, std::string pathbase, TH1D* HIST, TH2D* HIST2D = 0){
  std::stringstream cache_path;
  cache_path << "Cache_" << pathbase << ".csv";
  std::ifstream cache_file(cache_path.str().c_str());
  if(!cache_file.good()){
    cacheData_vegas(args, pathbase);
    cache_file.open(cache_path.str().c_str());
  }

  std::cout << "Loading " << pathbase << std::endl;

  std::string line;
  std::vector<std::string> line_fields;
  cuts_t cuts;
  while(std::getline(cache_file, line)){
    cuts.read++;
    bool cont = false;
    boost::split(line_fields, line, boost::is_any_of(","));
    std::vector<double> fields;
    for(auto it : line_fields){
      fields.push_back(std::stod(it));
    }
    line_fields.clear();

    //Zenith cut
    if(args.bin_vars & 1){
      if((int)fields[2] != ins.za){
        cuts.za++;
        cont = true;
      }
    }

    //Energy cut
    if(args.bin_vars & 2){
      if((int)fields[3] != ins.e){
        cuts.e++;
        cont = true;
      }
    }

    //Tel cut
    if(args.bin_vars & 4){
      if((int)fields[4] != ins.tel){
        cuts.tel++;
        cont = true;
      }
    }

    //Azimuth cut
    if(args.bin_vars & 8){
      if((int)fields[5] != ins.az){
        cuts.az++;
        cont = true;
      }
    }

    //Offset cut
    if(args.bin_vars & 16){
      if((int)fields[6] != ins.off){
        cuts.off++;
        cont = true;
      }
    }

    if(cont) continue;
    HIST->Fill(fields[0]);
    if(HIST2D) HIST2D->Fill(fields[0], fields[1]);
  }
  cuts.passed = HIST->Integral();
  std::cout << cuts.passed << " passed cuts." << std::endl;
  std::cout << cuts.read << " cuts read." << std::endl;
  std::cout << cuts.src << " failed src cut." << std::endl;
  std::cout << cuts.tel << " failed tel cut." << std::endl;
  std::cout << cuts.e << " failed energy cut." << std::endl;
  std::cout << cuts.za << " failed za cut." << std::endl;
  std::cout << cuts.msw << " failed msw cut." << std::endl;
  std::cout << cuts.az << " failed az cut." << std::endl;
  std::cout << cuts.off << " failed off cut." << std::endl;

  if(args.output & 8) print_cuts(pathbase, &cuts, OUTSTR);
  return 1;
}

void loadData(indices_t ins, args_t args, double *alpha, hists_t hists){
  OUTSTR = hists.outpath;
  if(args.format == Format_t::Toy){
    *alpha = loadData_toy(ins, args, "data", hists.dat_hist) / loadData_toy(ins, args, "bkg", hists.bkg_hist);
    loadsrc_csv(ins, args, hists.src_hist);
    std::cout << "Histograms loaded from Toy format." << std::endl;
  }
  else if(args.format == Format_t::Vegas){
    loadData_vegas(ins, args, "data", hists.dat_hist, hists.dat_2hist);
    //This is not ideal. (Slightly better now)
    if(!access("bkg_sources.list", F_OK)){
      std::ifstream flist("bkg_sources.list");
      std::string line;
      while(std::getline(flist, line)){
        loadData_vegas(ins, args, line, hists.bkg_hist, hists.bkg_2hist);
      }
    }
    else{
      loadData_vegas(ins, args, "bkg", hists.bkg_hist, hists.bkg_2hist);
    }
    *alpha = hists.dat_hist->Integral() / hists.bkg_hist->Integral(); //TODO
    loadsrc_csv(ins, args, hists.src_hist);
    std::cout << "Histograms loaded from Vegas format." << std::endl;
  }
  else if(args.format == Format_t::Sample){
    loadData_sample(ins, args, "data", hists.dat_hist, hists.dat_2hist);
    loadData_sample(ins, args, "bkg", hists.bkg_hist, hists.bkg_2hist);
    loadData_sample(ins, args, "src", hists.src_hist);
    *alpha = 3/5;
    std::cout << "Histograms loaded from Sample format." << std::endl;
  }
  else{
    std::cerr << "No valid data format specified." << std::endl;
    throw 999;
  }

  if(bkg_centers.size()){
    for(auto &i : bkg_centers){
      for(auto &j : bkg_centers){
        if(i != j){
          VACoordinatePair icoord = VACoordinatePair(i.first,i.second,VACoordinates::J2000,VACoordinates::Deg);
          VACoordinatePair jcoord = VACoordinatePair(j.first,j.second,VACoordinates::J2000,VACoordinates::Deg);
          if(icoord.angularSeparation_Deg(jcoord) < 3.5){
            std::cerr << "Background Overlap: Two background sources are within 3.5 degrees of each other." << std::endl;
          }
        }
      }
    }
  }


}
