#include "bbts.h"
#include "load_data.h"
#include "output.h"

std::string OUTPATH;
std::string LONGOUTPATH;

TH1D* DAT_HIST;
TH1D* BKG_HIST;
TH1D* SRC_HIST;

//Bin boundaries
double ZABINS[] = {0.00,14.1,25.46,32.59,37.57,42.56,47.55};
double EBINS[] = {316,630,1259,2512,5012};
int TBINS[] = {3,4};
double AZBINS[] = {0,45,90,135,180,225,270,315,360};
double OBINS[] = {0.0, 1.0, 1.4, 1.7, 2.0};


//helper function declarations
int parse_command_line(int argc, char* argv[], args_t* args);
void prepare_std_output_files(args_t args);
int optional_binning(indices_t indices, args_t args);

/**
  Calculates the negative log likelihood
  for the case with only a background
  component and barlow-beeston turned off.

  @param Pb Unweighted background fraction.
  @param F If present, F will hold the f_i histogram.
  @return Negative log likelihood.
*/
double nosrc_noBB(double Pb, TH1D* F){
  int i;
  double pb;
  double di, bi, fi;
  double dtot = 0, btot = 0, ftot = 0;
  std::ofstream f;
  double lnL = 0;
  pb = Pb * DAT_HIST->Integral() / BKG_HIST->Integral();
  for(i = 1; i <= NBIN; i++){
    di = DAT_HIST->GetBinContent(i);
    bi = BKG_HIST->GetBinContent(i);
    fi = pb * bi;
    if(di == 0) lnL -= fi;
    else        lnL += di * TMath::Log(fi) - fi;
    dtot += di; btot += bi; ftot += fi;
    if(F){
      F->Sumw2(false);
      F->SetBinContent(i, fi);
      F->Sumw2(true);
    }
  }
  return -lnL;
}

/**
  Calculates the negative log likelihood
  for the case with both background and source
  components and barlow-beeston turned off.

  @param Pb Unweighted background fraction.
  @param Ps Unweighted source fraction.
  @param F If present, F will hold the f_i histogram.
  @return Negative log likelihood.
*/
double src_noBB(double Pb, double Ps, TH1D* F){
  int i;
  double pb, ps;
  double di, bi, si, fi;
  double dtot = 0, btot = 0, ftot = 0, stot = 0;
  std::ofstream f;
  double lnL = 0;
  pb = Pb * DAT_HIST->Integral() / BKG_HIST->Integral();
  ps = Ps * DAT_HIST->Integral() / SRC_HIST->Integral();
  for(i = 1; i <= NBIN; i++){
    di = DAT_HIST->GetBinContent(i);
    bi = BKG_HIST->GetBinContent(i);
    si = SRC_HIST->GetBinContent(i);
    fi = pb * bi + ps * si;
    if(di == 0) lnL -= fi;
    else        lnL += di * TMath::Log(fi) - fi;

    dtot += di; btot += bi; ftot += fi; stot += si;
    if(F){
      F->Sumw2(false);
      F->SetBinContent(i, fi);
      F->Sumw2(true);
    }
  }
  return -lnL;
}

/**
  Calculates the negative log likelihood
  for the case with only a background
  component and barlow-beeston turned on.

  @param Pb Unweighted background fraction.
  @param F If present, F will hold the f_i histogram.
  @param B If present, B will hold the a_i histogram for background.
  @return Negative log likelihood.
*/
double nosrc_BB(double Pb, TH1D* F, TH1D* B){
  int i;
  double pb;
  double di, bi, Bi, fi;
  double Y, Z, ti;
  double dtot = 0, btot = 0, ftot = 0, Btot = 0;
  std::ofstream f;
  double lnL = 0;
  pb = Pb * DAT_HIST->Integral() / BKG_HIST->Integral();
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
  return -lnL;
}

/**
  Calculates the negative log likelihood
  for the case with both background and
  source components with barlow-beeston
  turned on.

  @param Pb Unweighted background fraction.
  @param Ps Unweighted source fraction.
  @param F If present, F will hold the f_i histogram.
  @param B If present, B will hold the a_i histogram for background.
  @return Negative log likelihood.
*/
double src_BB(double Pb, double Ps, TH1D* F, TH1D* B){
  int i;
  double pb, ps;
  double di, bi, si, Bi, fi;
  double dtot = 0, btot = 0, ftot = 0, stot = 0, Btot = 0;
  std::ofstream f;
  double X, Y, Z, ti;
  double lnL = 0;
  pb = Pb * DAT_HIST->Integral() / BKG_HIST->Integral();
  ps = Ps * DAT_HIST->Integral() / SRC_HIST->Integral();
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
  return -lnL;
}

void wrapper_nosrc_noBB(Int_t &nDim, Double_t *gout, Double_t &result, Double_t par[], Int_t flag){ result = nosrc_noBB(par[0]); }

void wrapper_src_noBB(Int_t &nDim, Double_t *gout, Double_t &result, Double_t par[], Int_t flag){ result = src_noBB(par[0], par[1], 0); }

void wrapper_nosrc_BB(Int_t &nDim, Double_t *gout, Double_t &result, Double_t par[], Int_t flag){ result = nosrc_BB(par[0]); }

void wrapper_src_BB(Int_t &nDim, Double_t *gout, Double_t &result, Double_t par[], Int_t flag){ result = src_BB(par[0], par[1], 0); }

/**
  Runs fits with and without source components
  and barlow-beeston turned on and off for a
  total of four fits. The fit results are then
  output to disk.

  @param ins Indices struct holding the current
  bin indices for each binned variable in the fit.
  @param args Arguments struct holding the
  command line arguments.
  @param fracs Array used to hold various fractions
  to be output.
  @param fit_param Which parameter is being fit on.
  Should have the value "MSW" or "BDT".
  @param hists Struct holding input histograms
*/
void fit(indices_t ins, args_t *args, Fit_Par_t fit_param, hists_t *hists){
  //Set up fitters
  TFitter* fit_nosrc_nobb = new TFitter(1);
  TFitter* fit_src_nobb   = new TFitter(2);
  TFitter* fit_nosrc_bb   = new TFitter(1);
  TFitter* fit_src_bb     = new TFitter(2);
  //Set Printout
  double p1 = args->verbosity;
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
  //Run MIGRAD
  fit_nosrc_nobb->ExecuteCommand("MIGRAD", 0, 0);
  fit_src_nobb  ->ExecuteCommand("MIGRAD", 0, 0);
  fit_nosrc_bb  ->ExecuteCommand("MIGRAD", 0, 0);
  fit_src_bb    ->ExecuteCommand("MIGRAD", 0, 0);
  //Run MINOS
  fit_nosrc_nobb->ExecuteCommand("MINOS", 0, 0);
  fit_src_nobb  ->ExecuteCommand("MINOS", 0, 0);
  fit_nosrc_bb  ->ExecuteCommand("MINOS", 0, 0);
  fit_src_bb    ->ExecuteCommand("MINOS", 0, 0);

  double fracs[4];
  fracs[0] = fit_src_nobb->GetParameter(0);
  fracs[1] = fit_src_nobb->GetParameter(1);
  fracs[2] = fit_src_bb->GetParameter(0);
  fracs[3] = fit_src_bb->GetParameter(1);
  write_to_root_file(args, hists, fracs, fit_param, 0);
  //TODO SRCEXCL
}

/**
  Loops over the binned variables and
  first loads data then runs the fit
  for each bin combination and each fit
  parameter.
*/
int main(int argc, char* argv[]){
  args_t* args = new args_t;
  if(parse_command_line(argc, argv, args)) return 1;
  init_output_root_file();
  indices_t indices;
  //TODO fix indices. Maybe make them more controllable
  //TODO Write method to create vector of indices to run. Then iterate.
  for(indices.za = 3; indices.za < 4; indices.za++){
    for(indices.e = 0; indices.e < 2; indices.e++){
      for(indices.tel = 0; indices.tel < 1; indices.tel++){
        for(indices.az = 0; indices.az < 8; indices.az++){
          for(indices.off = 0; indices.off < 4; indices.off++){
            if(optional_binning(indices, *args)) continue;
            TH1::SetDefaultSumw2();
            hists_t *hists = new hists_t;
            hists->outpath = OUTPATH;
            hists->longoutpath = LONGOUTPATH;
            Fit_Par_t fit_par;
            //MSW Fit
            if(args->fit_params & 1){
                hists->dat = new TH1D("MSWDataHist", "MSW Data", NBIN, MSWLOW, MSWHIGH);
                hists->bkg = new TH1D("MSWBkgHist", "MSW Bkg", NBIN, MSWLOW, MSWHIGH);
                hists->src = new TH1D("MSWSrcHist", "MSW Src", NBIN, MSWLOW, MSWHIGH);
                fit_par = Fit_Par_t::msw;
            }
            //BDT Fit
            if(args->fit_params & 2){
                hists->dat = new TH1D("BDTDataHist", "BDT Data", NBIN, BDTLOW, BDTHIGH);
                hists->bkg = new TH1D("BDTBkgHist", "BDT Bkg", NBIN, BDTLOW, BDTHIGH);
                hists->src = new TH1D("BDTSrcHist", "BDT Src", NBIN, BDTLOW, BDTHIGH);
                fit_par = Fit_Par_t::bdt;
            }
            loadData(indices, args, hists, fit_par);

            if(!hists->dat) throw 407;
            if(!hists->bkg) throw 408;
            if(!hists->src) throw 409;
            if(!hists->dat->Integral()) continue;
            if(!hists->bkg->Integral()) continue;
            if(!hists->src->Integral()) continue;

            //TODO Fix this. It's gross.
            DAT_HIST = hists->dat;
            BKG_HIST = hists->bkg;
            SRC_HIST = hists->src;
            OUTPATH = hists->outpath;
            LONGOUTPATH = hists->longoutpath;

            fit(indices, *args, fit_par, hists);

            //Clean up hists
            delete hists->dat;
            delete hists->bkg;
            delete hists->src;
          }
        }
      }
    }
  }
  return 0;
}

//Helpers
int parse_command_line(int argc, char* argv[], args_t* args){
  for(int i = 0; i < argc; i++){
    if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
      const char *help_text = R"(
OPTIONS:
  -b VAR, --bin-variable VAR
    Sets which variables to bin over.
    Available: zenith, energy, telescope, azimuth, offset, all.
    Default: zenith, energy, telescope.

  --fit-parameter
    Which parameter(s) to run the fit on.
    Available: msw, bdt
    Default: msw

  -h, --help
    Print help text and exit.

  -v VERBOSITY, --verbosity VERBOSITY
    Controls the verbosity of MINUIT.
    Default: -1.
      )";
      std::cout << help_text << std::endl;
      return 1;
    }
    if(!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbosity")){
      if(i < argc - 1){
        args->verbosity = atoi(argv[++i]);
      }
    }
    if(!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bin-variable")){
      if(i < argc - 1 && !strcmp(argv[i+1], "zenith") && !(args->bin_vars & 1)){
        if(args->bin_vars & 1) args->bin_vars -= 1;
        else args->bin_vars += 1;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "energy") && !(args->bin_vars & 2)){
        if(args->bin_vars & 2) args->bin_vars -= 2;
        else args->bin_vars += 2;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "telescope") && !(args->bin_vars & 4)){
        if(args->bin_vars & 4) args->bin_vars -= 4;
        else args->bin_vars += 4;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "azimuth") && !(args->bin_vars & 8)){
        if(args->bin_vars & 8) args->bin_vars -= 8;
        else args->bin_vars += 8;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "offset") && !(args->bin_vars & 16)){
        if(args->bin_vars & 16) args->bin_vars -= 16;
        else args->bin_vars += 16;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "all")){
        args->bin_vars = 31;
      }
    }
    if(!strcmp(argv[i], "--fit-parameter")){
      if(i < argc - 1 && !strcmp(argv[i+1], "msw")){
        if(args->fit_params & 1) args->fit_params -= 1;
        else args->fit_params += 1;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "bdt")){
        if(args->fit_params & 2) args->fit_params -= 2;
        else args->fit_params += 2;
      }
    }
  }
  return 0;
}

int optional_binning(indices_t indices, args_t args){
  std::stringstream path;
  std::stringstream longpath;
  if(!(args.bin_vars & 1) && indices.za != 0) return 1;
  else if(args.bin_vars & 1){
    path << "ZA" << indices.za;
    longpath << "ZA" << ZABINS[indices.za] << "-" << ZABINS[indices.za+1];
  }
  if(!(args.bin_vars & 2) && indices.e != 0) return 1;
  else if(args.bin_vars & 2){
    path << "E" << indices.e;
    longpath << "_E" << EBINS[indices.e] << "-" << EBINS[indices.e+1];
  }
  if(!(args.bin_vars & 4) && indices.tel != 0) return 1;
  else if(args.bin_vars & 4){
    path << "T" << TBINS[indices.tel];
    longpath << "_T" << TBINS[indices.tel];
  }
  if(!(args.bin_vars & 8) && indices.az != 0) return 1;
  else if(args.bin_vars & 8){
    path << "A" << indices.az;
    longpath << "_A" << AZBINS[indices.az] << "-" << AZBINS[indices.az+1];
  }
  if(!(args.bin_vars & 16) && indices.off != 0) return 1;
  else if(args.bin_vars & 16){
    path << "O" << indices.off;
    longpath << "_O" << OBINS[indices.off] << "-" << OBINS[indices.off+1];
  }
  std::cout << path.str() << std::endl;
  std::cout << longpath.str() << std::endl;
  OUTPATH = path.str();
  LONGOUTPATH = longpath.str();
  return 0;
}
