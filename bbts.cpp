#include "bbts.h"
#include "load_data.h"
#include "output.h"


//This is a global variable so that it doesn't have to
//be passed explicitly into the likelihood functions.
hists_t* HISTS;

//Bin boundaries
double ZABINS[] = {0.00,14.1,25.46,32.59,37.57,42.56,47.55};
double EBINS[] = {316,630,1259,2512,5012};
int TBINS[] = {3,4};
double AZBINS[] = {0,45,90,135,180,225,270,315,360};
double OBINS[] = {0.0, 1.0, 1.4, 1.7, 2.0};


//helper function declarations
int parse_command_line(int argc, char* argv[], args_t* args);
void get_bins(vector<indices_t> *ins_list, args_t *args);

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
  pb = Pb * HISTS->dat->Integral() / HISTS->bkg->Integral();
  for(i = 1; i <= NBIN; i++){
    di = HISTS->dat->GetBinContent(i);
    bi = HISTS->bkg->GetBinContent(i);
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
  pb = Pb * HISTS->dat->Integral() / HISTS->bkg->Integral();
  ps = Ps * HISTS->dat->Integral() / HISTS->src->Integral();
  for(i = 1; i <= NBIN; i++){
    di = HISTS->dat->GetBinContent(i);
    bi = HISTS->bkg->GetBinContent(i);
    si = HISTS->src->GetBinContent(i);
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
  pb = Pb * HISTS->dat->Integral() / HISTS->bkg->Integral();
  for(i = 1; i <= NBIN; i++){
    di = HISTS->dat->GetBinContent(i);
    bi = HISTS->bkg->GetBinContent(i);

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
  pb = Pb * HISTS->dat->Integral() / HISTS->bkg->Integral();
  ps = Ps * HISTS->dat->Integral() / HISTS->src->Integral();
  for(i = 1; i <= NBIN; i++){
    di = HISTS->dat->GetBinContent(i);
    bi = HISTS->bkg->GetBinContent(i);
    si = HISTS->src->GetBinContent(i);

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
void fit(indices_t *ins, args_t *args, Fit_Par_t fit_param){
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
  write_to_root_file(ins, HISTS, fracs, fit_param);
  //TODO SRCEXCL
}

void run_fit(indices_t *ins, args_t *args, Fit_Par_t fit_param){
  loadData(ins, HISTS, fit_param);
  //Check that everything was loaded correctly
  if(!HISTS->dat) throw 407;
  if(!HISTS->bkg) throw 408;
  if(!HISTS->src) throw 409;
  if(!HISTS->dat->Integral()) return;
  if(!HISTS->bkg->Integral()) return;
  if(!HISTS->src->Integral()) return;
  fit(ins, args, fit_param);
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
  std::vector<indices_t> *ins_list = new std::vector<indices_t>;
  get_bins(ins_list, args);
  for(indices_t &ins : *ins_list){
    TH1::SetDefaultSumw2();
    hists_t *hists = new hists_t;
    //MSW Fit
    if(args->fit_params & 1){
      hists->dat = new TH1D("MSWDataHist", "MSW Data", NBIN, MSWLOW, MSWHIGH);
      hists->bkg = new TH1D("MSWBkgHist", "MSW Bkg", NBIN, MSWLOW, MSWHIGH);
      hists->src = new TH1D("MSWSrcHist", "MSW Src", NBIN, MSWLOW, MSWHIGH);

      HISTS = hists;
      run_fit(&ins, args, Fit_Par_t::msw);
      //Clean up hists
      delete hists->dat;
      delete hists->bkg;
      delete hists->src;
    }
    //BDT Fit
    if(args->fit_params & 2){
      hists->dat = new TH1D("BDTDataHist", "BDT Data", NBIN, BDTLOW, BDTHIGH);
      hists->bkg = new TH1D("BDTBkgHist", "BDT Bkg", NBIN, BDTLOW, BDTHIGH);
      hists->src = new TH1D("BDTSrcHist", "BDT Src", NBIN, BDTLOW, BDTHIGH);

      HISTS = hists;
      run_fit(&ins, args, Fit_Par_t::bdt);
      //Clean up hists
      delete hists->dat;
      delete hists->bkg;
      delete hists->src;
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

void get_bins(vector<indices_t> *ins_list, args_t *args){
  int za_max = 6;
  int e_max = 4;
  int tel_max = 2;
  int az_max = 8;
  int off_max = 4;

  for(int se = 0; se < 2; se ++){
    for(int z = 3; z < za_max; z++){
      for(int e = -1; e < e_max; e++){
        for(int t = -1; t < tel_max; t++){
          for(int a = -1; a < az_max; a++){
            for(int o = -1; o < off_max; o++){

              if(args->bin_vars&1 && z == -1) continue;
              else if(!(args->bin_vars&1) && z != -1) continue;

              if(args->bin_vars&2 && e == -1) continue;
              else if(!(args->bin_vars&2) && e != -1) continue;

              if(args->bin_vars&4 && t == -1) continue;
              else if(!(args->bin_vars&4) && t != -1) continue;

              if(args->bin_vars&8 && a == -1) continue;
              else if(!(args->bin_vars&8) && a != -1) continue;

              if(args->bin_vars&16 && o == -1) continue;
              else if(!(args->bin_vars&16) && o != -1) continue;

              indices_t ins;
              ins.za = z;
              ins.e = e;
              ins.tel = t;
              ins.az = a;
              ins.off = o;
              ins.src_excl = se;
              ins_list->push_back(ins);
            }
          }
        }
      }
    }
  }

}

