#include "load_data.h"
#include "output.h"

std::string OUTSTR;

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


/*  Return value is a bit vector.
 *  bit 0 - Passed cuts for msw
 *  bit 1 - Passed cuts for bdt
 *
 *  Ex. rval = 2 -> event is suitable for bdt fit but not msw.
 */
int bin_event(Cut_Params_t *params, cuts_t *cuts, args_t *args, indices_t *ins, std::vector<SourceCut_t> &source_cuts){
  bool fail = false;

  //Source Cut
  if(source_cuts.size()){
    VACoordinatePair eventCoord = VACoordinatePair(
      params->eventRa,
      params->eventDec,
      VACoordinates::J2000,
      VACoordinates::Deg);
    bool fail_src_cuts;
    for(auto &it: source_cuts){
      fail_src_cuts = it.InsideExclRadius(eventCoord);
      if(fail_src_cuts) break;
    }
    if(fail_src_cuts){
      cuts->src++;
      fail = true;
    }
  }

  //ZA cut
  int z_bin = -42;
  if(args->bin_vars & 1){
    for(int z = 0; z < 6; z++){
      if(params->zenith < ZABINS[z+1] && params->zenith > ZABINS[z]){
        z_bin = z;
      }
    }
    if(z_bin == -42){
      cuts->za++;
      fail = true;
    }
  }

  //Energy Cut
  int e_bin = -42;
  if(args->bin_vars & 2){
    for(int e = 0; e < 4; e++){
      if(params->energy < EBINS[e+1] && params->energy > EBINS[e]){
        e_bin = e;
      }
    }
    if(e_bin == -42){
      cuts->e++;
      fail = true;
    }
  }

  //Tel Cut
  if(args->bin_vars & 4){
    if(params->ntels == 2){
      cuts->tel++;
      fail = true;
    }
  }

  //AZ Cut
  int a_bin = -42;
  if(args->bin_vars & 8){
    for(int a = 0; a < 8; a++){
      if(params->azimuth < AZBINS[a+1] && params->azimuth > AZBINS[a]){
        a_bin = a;
      }
    }
    if(a_bin == -42){
      cuts->az++;
      fail = true;
    }
  }

  //Offset Cut
  int o_bin = -42;
  if(args->bin_vars & 16){
    for(int o = 0; o < 8; o++){
      if(params->offset < OBINS[o+1] && params->offset > OBINS[o]){
        o_bin = o;
      }
    }
    if(o_bin == -42){
      cuts->off++;
      fail = true;
    }
  }
  else if(params->offset > 1.75){
    cuts->off++;
    fail = true;
  }

  params->zbin = z_bin;
  params->ebin = e_bin;
  params->abin = a_bin;
  params->obin = o_bin;

  if(fail) return 0;

  int rval = 0;
  bool msw_fail = false;
  bool bdt_fail = false;

  //MSL Cut
  if(params->msl > 1.3 || params->msl < 0.05) msw_fail = true;

  //Shower Height Cut
  if(params->height < 7) msw_fail = true;

  //MSW Cut (No need to carry around over/underflow in hists)
  if(params->msw < MSWLOW || params->msw > MSWHIGH) msw_fail = true;

  //BDT Cut (See above)
  if(params->bdt < BDTLOW || params->bdt > BDTHIGH) bdt_fail = true;

  if(!msw_fail) rval += 1;
  if(!bdt_fail) rval += 2;

  return rval;
}


bool check_event(TTreeReader *reader, cuts_t *cuts, args_t *args){
  TTreeReaderValue<float> zenith(reader, "Zenith");
  TTreeReaderValue<float> energy(reader, "Energy_GeV");
  TTreeReaderValue<float> ntels(reader, "NTels");
  TTreeReaderValue<float> azimuth(reader, "Azimuth");
  TTreeReaderValue<float> offset(reader, "Offset");
  TTreeReaderValue<float> src_excl(reader, "Src_Excl");

  bool pass = true;

  //ZA Bin
  if(args->bin_vars & 1){
    if(*zenith >= ZABINS[ins->za] || *zenith < ZABINS[ins->za]){
      pass = false;
      cuts->za++;
    }
  }

  //Energy Bin
  if(args->bin_vars & 2){
    if(*energy >= EBINS[ins->e] || *energy < EBINS[ins->e]){
      pass = false;
      cuts->e++;
    }
  }

  //Tel Bin
  if(args->bin_vars & 4){
    if(*ntels != TBINS[ins->tel]){
      pass = false;
      cuts->tel++;
    }
  }

  //Azimuth Bin
  if(args->bin_vars & 8){
    if(*azimuth >= AZBINS[ins->az] || *azimuth < AZBINS[ins->az]){
      pass = false;
      cuts->az++;
    }
  }

  //Offset Bin
  if(args->bin_vars & 16){
    if(*offset >= OBINS[ins->off] || *offset < OBINS[ins->off]){
      pass = false;
      cuts->off++;
    }
  }
  //Offset Cut
  else if(*offset > MAX_OFFSET){
    pass = false;
    cuts->off++;
  }

  //Source Cut
  if(*src_excl){
    pass = false;
    cuts->src++;
  }

  return pass;
}

bool check_event_msw(TTreeReader *reader, cuts_t *cuts, args_t *args){
  TTreeReaderValue<float> msw(reader, "MSW");
  TTreeReaderValue<float> msl(reader, "MSL");
  TTreeReaderValue<float> height(reader, "ShowerHeightMax_KM");

  bool pass = true;

  //MSL Cut
  if(*msl > MAX_MSL || *msl < MIN_MSL){
    pass = false;
    cuts->msl++;
  }

  //Shower Height Cut
  if(*height < MIN_HEIGHT){
    pass = false;
    cuts->height++;
  }

  //MSW Cut
  if(*msw > MSWHIGH || *msw < MSWLOW){
    pass = false;
    cuts->par++;
  }

  return pass;
}

bool check_event_bdt(TTreeReader *reader, cuts_t *cuts, args_t *args){
  TTreeReaderValue<float> bdt(reader, "BDT");

  bool pass = true;

  //BDT Cut
  if(*bdt > BDTHIGH || *bdt < BDTLOW){
    pass = false;
    cuts->par++;
  }

  return pass;
}


double loadData_bbmlm(indices_t *ins, args_t *args, TH1D *hist, Fit_Par_t fit_par, std::string infile){
  TFile *infile = TFile::Open(infile.c_str());
  if(!infile->IsOpen()){
    std::cerr << "Cannot open " << path.str() << std::endl;
    return -1;
  }
  else{
    std::cout << "Reading " << path.str() << std::endl;
  }

  TTreeReader reader("ntuple", infile);
  TTreeReaderValue<float> msw(reader, "MSW");
  TTreeReaderValue<float> bdt(reader, "BDT");

  cuts_t cuts;
  timestamp;
  while(reader.Next()){
    cuts.read++;

    bool pass = false;
    if(fit_par == Fit_Par_t::msw){
      pass = check_event(&reader, &cuts, args)
           & check_event_msw(&reader, &cuts, args);
      if(pass) hist->Fill(*msw);
    }
    else if(fit_par == Fit_Par_t::bdt){
      pass = check_event(&reader, &cuts, args)
           & check_event_bdt(&reader, &cuts, args);
      if(pass) hist->Fill(*bdt);
    }

  }
  std::cout << hist->Integral() << " passed cuts." << std::endl;
  std::cout << cuts.read << " cuts read." << std::endl;
  std::cout << cuts.src << " failed src cut." << std::endl;
  std::cout << cuts.par << " failed parameter cut." << std::endl;
  std::cout << cuts.msl << " failed msl cut." << std::endl;
  std::cout << cuts.height << " failed height cut." << std::endl;
  std::cout << cuts.za << " failed za cut." << std::endl;
  std::cout << cuts.e << " failed energy cut." << std::endl;
  std::cout << cuts.tel << " failed tel cut." << std::endl;
  std::cout << cuts.az << " failed az cut." << std::endl;
  std::cout << cuts.off << " failed off cut." << std::endl;


  return 1;
}

void loadData(indices_t *ins, args_t *args, hists_t *hists, Fit_Par_t fit_par){
  OUTSTR = hists->outpath;
  loadData_bbmlm(ins, args, hists->dat, fit_par, "data.root");
  loadData_bbmlm(ins, args, hists->bkg, fit_par, "bkg.root");
  loadData_bbmlm(ins, args, hists->src, fit_par, "src.root");
  std::cout << "Histograms loaded." << std::endl;
}
