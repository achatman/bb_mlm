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

typedef struct Cut_Params{
  double eventRa;
  double eventDec;
  double zenith;
  double energy;
  double ntels;
  double azimuth;
  double offset;
  double msw;
  double msl;
  double height;
  double bdt;
  int zbin;
  int ebin;
  int abin;
  int obin;
} Cut_Params_t;


/*  Return value is a bit vector.
 *  bit 0 - Passed cuts for msw
 *  bit 1 - Passed cuts for bdt
 *
 *  Ex. rval = 2 -> event is suitable for bdt fit but not msw.
 */
int bin_event(Cut_Params_t *params, cuts_t *cuts, args_t *args, std::vector<SourceCut_t> &source_cuts){
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

double loadData_bbmlm(indices_t ins, args_t args, std::string pathbase, hists_t *hists){
  std::stringstream path;
  path << "Data_" << pathbase << ".root";
  TFile *infile = TFile::Open(path.str().c_str());
  if(!infile->IsOpen()){
    std::cerr << "Cannot open " << path.str() << std::endl;
    return -1;
  }
  else{
    std::cout << "Reading " << path.str() << std::endl;
  }

  TTreeReader reader("ntuple", infile);
  TTreeReaderValue<float> ra(reader, "RA");
  TTreeReaderValue<float> dec(reader, "DEC");
  TTreeReaderValue<float> zenith(reader, "Zenith");
  TTreeReaderValue<float> energy(reader, "Energy_GeV");
  TTreeReaderValue<float> ntels(reader, "NTels");
  TTreeReaderValue<float> azimuth(reader, "Azimuth");
  TTreeReaderValue<float> offset(reader, "Offset");
  TTreeReaderValue<float> msw(reader, "MSW");
  TTreeReaderValue<float> msl(reader, "MSL");
  TTreeReaderValue<float> bdt(reader, "BDT");
  TTreeReaderValue<float> height(reader, "ShowerHeightMax_KM");

  cuts_t cuts;
  timestamp;
  while(reader.Next()){
    cuts.read++;
    //Could probably just pass the TTreeReaderValue
    //Not sure if that would work with the other functions though so
    //This is how it is for now
    Cut_Params_t *params = new Cut_Params_t();
    params->eventRa = *ra;
    params->eventDec = *dec;
    params->zenith = *zenith;
    params->energy = *energy;
    params->ntels = *ntels;
    params->azimuth = *azimuth;
    params->offset = *offset;
    params->msw = *msw;
    params->msl = *msl;
    params->bdt = *bdt;
    params->height = *height;

    vector<SourceCut_t> empty;
    int r = bin_event(params, &cuts, &args, empty);
    if(r == 0) continue;

    cuts.passed++;
    if(pathbase == "data"){
      if(r&1 && hists->msw_dat) hists->msw_dat->Fill(*msw);
      if(r&2 && hists->bdt_dat) hists->bdt_dat->Fill(*bdt);
    }
    else if(pathbase == "bkg"){
      if(r&1 && hists->msw_bkg) hists->msw_bkg->Fill(*msw);
      if(r&2 && hists->bdt_bkg) hists->bdt_bkg->Fill(*bdt);
    }
    else if(pathbase == "src"){
      if(r&1 && hists->msw_src) hists->msw_src->Fill(*msw);
      if(r&2 && hists->bdt_src) hists->bdt_src->Fill(*bdt);
    }
  }
  std::cout << cuts.passed << " passed cuts." << std::endl;
  std::cout << cuts.read << " cuts read." << std::endl;
  std::cout << cuts.src << " failed src cut." << std::endl;
  std::cout << cuts.tel << " failed tel cut." << std::endl;
  std::cout << cuts.e << " failed energy cut." << std::endl;
  std::cout << cuts.za << " failed za cut." << std::endl;
  std::cout << cuts.az << " failed az cut." << std::endl;
  std::cout << cuts.off << " failed off cut." << std::endl;


  return 1;
}

void loadData(indices_t ins, args_t args, hists_t *hists){
  OUTSTR = hists->outpath;
  loadData_bbmlm(ins, args, "data", hists);
  loadData_bbmlm(ins, args, "bkg", hists);
  loadData_bbmlm(ins, args, "src", hists);
  std::cout << "Histograms loaded." << std::endl;
}
