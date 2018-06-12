#include "bbts.h"
#include "load_data.h"

std::string OUTPATH;

TH1F* DAT_HIST;
TH1F* BKG_HIST;
TH1F* SRC_HIST;

//Bin boundaries
double ZABINS[] = {25.46,32.59,37.57,42.56,47.55};
double EBINS[] = {316,630,1259,2512,5012};
int TBINS[] = {3,4};
double AZBINS[] = {0,45,90,135,180,225,270,315,360};
double OBINS[] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};


//helper function declarations
int parse_command_line(int argc, char* argv[], args_t* args);
void prepare_std_output_files(args_t args);
int optional_binning(indices_t indices, args_t args);
//output function declarations
void printRawData();
void histogram_raw_data(indices_t ins);
void histogram_fit_data(double fracs[6], indices_t ins);
void calculate_errors(double Pb, double Ps, double sigma_Pb, double sigma_Ps, indices_t ins, double alpha);
void print_cuts(std::string action, cuts_t* cuts);
void fit_manual_bin();
void map_likelihood(double Pb, double Ps, std::string title_tag, indices_t ins, args_t args);
void plot_msw_vs_msl(TH2F* dat_2hist, TH2F* bkg_2hist);


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

void fit(indices_t ins, args_t args, double alpha, double *fracs = 0){
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

  if(fracs){
    fracs[0] = fit_src_nobb->GetParameter(0);
    fracs[1] = fit_src_nobb->GetParameter(1);
    fracs[2] = fit_src_bb->GetParameter(0);
    fracs[3] = fit_src_bb->GetParameter(1);
    fracs[4] = fit_nosrc_nobb->GetParameter(0);
    fracs[5] = fit_nosrc_bb->GetParameter(0);
  }

  //Get likelihood and TS
  bool output_bins = (args.output & 1) && (DAT_HIST->Integral() + BKG_HIST->Integral());
  double lnL_nosrc_nobb = -nosrc_noBB(fit_nosrc_nobb->GetParameter(0), output_bins);
  double lnL_src_nobb = -src_noBB(fit_src_nobb->GetParameter(0), fit_src_nobb->GetParameter(1), output_bins);
  double lnL_nosrc_bb = -nosrc_BB(fit_nosrc_bb->GetParameter(0), output_bins);
  double lnL_src_bb = -src_BB(fit_src_bb->GetParameter(0), fit_src_bb->GetParameter(1), output_bins);
  double TS_nobb = 2 * (lnL_src_nobb - lnL_nosrc_nobb);
  double TS_bb = 2 * (lnL_src_bb - lnL_nosrc_bb);
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

void run(int argc, char* argv[]){
  args_t* args = new args_t;
  if(parse_command_line(argc, argv, args)) return;
  prepare_std_output_files(*args);
  int zi = 1, ei = 0, ti = 0, ai = 0, oi = 0;
  indices_t indices;
  for(indices.za = zi; indices.za < 2; indices.za++){
    for(indices.e = ei; indices.e < 4; indices.e++){
      for(indices.tel = ti; indices.tel < 2; indices.tel++){
        for(indices.az = ai; indices.az < 8; indices.az++){
          for(indices.off = oi; indices.off < 8; indices.off++){
            if(optional_binning(indices, *args)) continue;

            //Run Fit
            double alpha = 1;
            TH1::SetDefaultSumw2();
            DAT_HIST = new TH1F("DataHist", "Data", NBIN, MSWLOW, MSWHIGH);
            BKG_HIST = new TH1F("BkgHist", "BKG", NBIN, MSWLOW, MSWHIGH);
            SRC_HIST = new TH1F("SrcHist", "SRC", NBIN, MSWLOW, MSWHIGH);
            TH2F *dat_2hist, *bkg_2hist;
            if(args->graphics & 4){
              dat_2hist = new TH2F("Data_MSWvsMSL", "Data MSW vs MSL", NBIN, MSWLOW, MSWHIGH, NBIN, MSWLOW, MSWHIGH);
              bkg_2hist = new TH2F("Bkg_MSWvsMSL", "Bkg MSW vs MSL", NBIN, MSWLOW, MSWHIGH, NBIN, MSWLOW, MSWHIGH);
            }
            loadData(indices, *args, &alpha, DAT_HIST, BKG_HIST, SRC_HIST, dat_2hist, bkg_2hist);

            if(!DAT_HIST || !BKG_HIST || !SRC_HIST) throw 407;
            double fracs[6];
            fit(indices, *args, alpha, fracs);
            if(args->output & 2) printRawData();
            if(args->hist & 1) histogram_raw_data(indices);
            if(args->hist & 2) histogram_fit_data(fracs, indices);
            if(args->graphics & 4) plot_msw_vs_msl(dat_2hist, bkg_2hist);
          }
        }
      }
    }
  }
  delete args;
}

void bidirectional(int argc, char* argv[]){
  args_t* args = new args_t;
  if(parse_command_line(argc, argv, args)) return;
  prepare_std_output_files(*args);
  int zi = 1, ei = 0, ti = 0, ai = 0, oi = 0;
  indices_t indices;
  for(indices.za = zi; indices.za < 2; indices.za++){
    for(indices.e = ei; indices.e < 4; indices.e++){
      for(indices.tel = ti; indices.tel < 2; indices.tel++){
        for(indices.az = ai; indices.az < 8; indices.az++){
          for(indices.off = oi; indices.off < 8; indices.off++){
            if(optional_binning(indices, *args)) continue;
            //Run Forward Fit
            double alpha = 1;
            TH1::SetDefaultSumw2();
            DAT_HIST = new TH1F("DataHist", "Data", NBIN, MSWLOW, MSWHIGH);
            BKG_HIST = new TH1F("BkgHist", "BKG", NBIN, MSWLOW, MSWHIGH);
            SRC_HIST = new TH1F("SrcHist", "SRC", NBIN, MSWLOW, MSWHIGH);
            loadData(indices, *args, &alpha, DAT_HIST, BKG_HIST, SRC_HIST);
            if(!DAT_HIST || !BKG_HIST || !SRC_HIST) throw 407;
            double fracs_for[6];
            fit(indices, *args, alpha, fracs_for);

            //Plot Forwards
            gStyle->SetOptStat(0);
            TCanvas *c1 = new TCanvas("",OUTPATH.c_str(),1600,1600);
            c1->Divide(3,2);
            std::stringstream title;

            //Plot Forward Fit
            c1->cd(1);
            TH1F* dat_fit_forward = new TH1F("DFit_For", "Forward", NBIN, MSWLOW, MSWHIGH);
            TH1F* bkg_fit_forward = new TH1F("BFit_For", "Forward", NBIN, MSWLOW, MSWHIGH);
            dat_fit_forward->SetLineColor(4);
            bkg_fit_forward->SetLineColor(6);
            dat_fit_forward->SetStats(false);
            bkg_fit_forward->SetStats(false);
            title << OUTPATH << " " << "Fit Forward";
            dat_fit_forward->SetTitle(title.str().c_str());
            title.str("");
            src_BB(fracs_for[2], fracs_for[3], false, dat_fit_forward, bkg_fit_forward);
            bkg_fit_forward->Scale(dat_fit_forward->Integral() / bkg_fit_forward->Integral());
            TLegend *legend2 = new TLegend(0.12, 0.8, 0.3, 0.9);
            legend2->AddEntry(dat_fit_forward, "Fit Data");
            legend2->AddEntry(bkg_fit_forward, "Fit Bkg");
            dat_fit_forward->SetMinimum(0);
            dat_fit_forward->SetMaximum(std::max(dat_fit_forward->GetMaximum(), bkg_fit_forward->GetMaximum())*1.1);
            dat_fit_forward->Draw();
            dat_fit_forward->Draw("sameE0");
            bkg_fit_forward->Draw("same");
            bkg_fit_forward->Draw("sameE0");
            legend2->Draw();

            //store ts for later
            double ts1 = -2*(src_BB(fracs_for[2], fracs_for[3]) - nosrc_BB(fracs_for[5]));

            //Run Backward Fit
            TH1F* temp = DAT_HIST;
            DAT_HIST = BKG_HIST;
            BKG_HIST = temp;
            double fracs_back[6];
            fit(indices, *args, alpha, fracs_back);

            //Plot Backward Fit
            c1->cd(4);
            TH1F* dat_fit_back = new TH1F("DFit_Back", "Backward", NBIN, MSWLOW, MSWHIGH);
            TH1F* bkg_fit_back = new TH1F("BFit_Back", "Backward", NBIN, MSWLOW, MSWHIGH);
            dat_fit_back->SetLineColor(4);
            bkg_fit_back->SetLineColor(6);
            dat_fit_back->SetStats(false);
            bkg_fit_back->SetStats(false);
            title << OUTPATH << " " << "Fit Backward";
            dat_fit_back->SetTitle(title.str().c_str());
            title.str("");
            src_BB(fracs_back[2], fracs_back[3], false, dat_fit_back, bkg_fit_back);
            bkg_fit_back->Scale(dat_fit_back->Integral() / bkg_fit_back->Integral());
            TLegend *legend4 = new TLegend(0.12, 0.8, 0.3, 0.9);
            legend4->AddEntry(dat_fit_back, "Fit Data");
            legend4->AddEntry(bkg_fit_back, "Fit Bkg");
            dat_fit_back->SetMinimum(0);
            dat_fit_back->SetMaximum(std::max(dat_fit_back->GetMaximum(), bkg_fit_back->GetMaximum())*1.1);
            dat_fit_back->Draw();
            dat_fit_back->Draw("sameE0");
            bkg_fit_back->Draw("same");
            bkg_fit_back->Draw("sameE0");
            legend4->Draw();

            //Plot Raw Comp & Pulls
            c1->cd(2);
            TH1F* dat1 = new TH1F(*BKG_HIST);
            TH1F* dat2 = new TH1F(*DAT_HIST);
            dat1->SetLineColor(1);
            dat2->SetLineColor(2);
            dat1->SetTitle("Raw Comparison");
            dat2->Scale(dat1->Integral() / dat2->Integral());
            TRatioPlot_BetterError* rp_raw = new TRatioPlot_BetterError(dat1, dat2, "diffsig");
            rp_raw->SetH1DrawOpt("E0");
            rp_raw->SetH2DrawOpt("E0");
            rp_raw->Draw();
            TLegend *legend5 = new TLegend(0.12, 0.8, 0.3, 0.9);
            legend5->AddEntry(dat1, "Sample 1");
            legend5->AddEntry(dat2, "Sample 2");
            legend5->Draw();

            //Plot Fit Comp & Pulls
            c1->cd(5);
            TH1F* fit1 = new TH1F(*dat_fit_forward);
            TH1F* fit2 = new TH1F(*dat_fit_back);
            fit1->SetLineColor(1);
            fit2->SetLineColor(2);
            fit1->SetTitle("Fit Comparison");
            fit2->Scale(fit1->Integral() / fit2->Integral());
            TRatioPlot_BetterError* rp_fit = new TRatioPlot_BetterError(fit1, fit2, "diffsig");
            rp_fit->SetH1DrawOpt("E0");
            rp_fit->SetH2DrawOpt("E0");
            rp_fit->Draw();
            TLegend *legend6 = new TLegend(0.12, 0.8, 0.3, 0.9);
            legend6->AddEntry(fit1, "Sample 1");
            legend6->AddEntry(fit2, "Sample 2");
            legend6->Draw();

            //Write Forward Fit Data
            c1->cd(3);
            TPaveText *pt1 = new TPaveText(0, 0, 1, 1);
            std::stringstream line;
            pt1->AddText(.5, .95, "Forward Fit Values:");
            line << "P_b = " << fracs_for[2];
            pt1->AddText(.05, .85, line.str().c_str());
            line.str("");
            line << "P_s = " << fracs_for[3];
            pt1->AddText(.05, .8, line.str().c_str());
            line.str("");
            line << "TS = " << ts1;
            pt1->AddText(.05, .75, line.str().c_str());
            line.str("");
            line << "Data Count = " << BKG_HIST->Integral();
            pt1->AddText(.05, .7, line.str().c_str());
            line.str("");
            line << "Bkg Count = " << DAT_HIST->Integral();
            pt1->AddText(.05, .65, line.str().c_str());
            line.str("");
            pt1->AddLine(0, .5, 1, .5);
            pt1->AddText(.5, .4, "Bin Boundaries:");
            if(args->bin_vars & 1){
              line << "Zenith Angle: " << ZABINS[indices.za] << "-" << ZABINS[indices.za+1];
              pt1->AddText(.05, .3, line.str().c_str());
              line.str("");
            }
            if(args->bin_vars & 2){
              line << "Energy: " << EBINS[indices.e] << "-" << EBINS[indices.e+1];
              pt1->AddText(.05, .25, line.str().c_str());
              line.str("");
            }
            if(args->bin_vars & 4){
              line << "Telescope: " << TBINS[indices.tel];
              pt1->AddText(.05, .2, line.str().c_str());
              line.str("");
            }
            if(args->bin_vars & 8){
              line << "Azimuth: " << AZBINS[indices.az] << "-" << AZBINS[indices.az+1];
              pt1->AddText(.05, .15, line.str().c_str());
              line.str("");
            }
            if(args->bin_vars & 16){
              line << "Offset: " << OBINS[indices.off] << "-" << OBINS[indices.off+1];
              pt1->AddText(.05, .1, line.str().c_str());
              line.str("");
            }
            pt1->SetAllWith("=", "size", .05);
            pt1->SetAllWith("=", "align", 12);
            pt1->SetAllWith(": ", "size", .05);
            pt1->SetAllWith(": ", "align", 12);
            pt1->Draw();


            //Write Backward Fit Data
            c1->cd(6);
            double ts2 = -2*(src_BB(fracs_back[2], fracs_back[3]) - nosrc_BB(fracs_back[5]));
            TPaveText *pt2 = new TPaveText(0, 0, 1, 1);
            pt2->AddText(.5, .95, "Backward Fit Values:");
            line << "P_b = " << fracs_back[2];
            pt2->AddText(.05, .85, line.str().c_str());
            line.str("");
            line << "P_s = " << fracs_back[3];
            pt2->AddText(.05, .8, line.str().c_str());
            line.str("");
            line << "TS = " << ts2;
            pt2->AddText(.05, .75, line.str().c_str());
            line.str("");
            line << "Data Count = " << DAT_HIST->Integral();
            pt2->AddText(.05, .7, line.str().c_str());
            line.str("");
            line << "Bkg Count = " << BKG_HIST->Integral();
            pt2->AddText(.05, .65, line.str().c_str());
            line.str("");
            pt2->AddLine(0, .5, 1, .5);
            if(args->op_info.c_str()){
              std::ifstream infile(args->op_info);
              std::string readline;
              std::getline(infile, readline);
              pt2->AddText(.05, .4, readline.c_str())->SetTextAlign(12);
              std::getline(infile, readline);
              pt2->AddText(.05, .35, readline.c_str())->SetTextAlign(12);
              std::getline(infile, readline);
              pt2->AddText(.05, .3, readline.c_str())->SetTextAlign(12);
              std::getline(infile, readline);
              pt2->AddText(.05, .25, readline.c_str())->SetTextAlign(12);
              std::getline(infile, readline);
              pt2->AddText(.05, .2, readline.c_str())->SetTextAlign(12);
              std::getline(infile, readline);
              pt2->AddText(.05, .15, readline.c_str())->SetTextAlign(12);
              infile.close();
            }

            pt2->SetAllWith("=", "size", .05);
            pt2->SetAllWith("=", "align", 12);
            pt2->Draw();


            //Save
            title << "Bidirectional_" << OUTPATH << ".png";
            c1->SaveAs(title.str().c_str());
            c1->Clear();
            delete c1;
            delete DAT_HIST;
            delete BKG_HIST;
            delete SRC_HIST;
            delete dat_fit_forward;
            delete bkg_fit_forward;
            delete dat_fit_back;
            delete bkg_fit_back;
            delete legend2;
            delete legend4;
            delete legend5;
            delete legend6;
            delete dat1;
            delete dat2;
            delete fit1;
            delete fit2;


          }
        }
      }
    }
  }
}


int main(int argc, char** argv){
  if(!strcmp("--bidirectional", argv[1])) bidirectional(argc, argv);
  else run(argc, argv);
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

  --bidirectional
    Runs the fit normally, then swaps data and bkg and runs the fit again.

  -d FORMAT, --data-format FORMAT
    Format of the imput data.
    Available: toy, csv, vegas, sample.
    Default: toy.

  -g GRAPHICS, --graphics GRAPHICS
    Triggers output of graphics files.
    Available: none, stdlnL, bblnL, mswmsl, all.
    Default: none.

  -h, --help
    Print help text and exit.

  -hist DATA, --histogram DATA
    Triggers output of histograms.
    Available: none, raw, fit, all.
    Default: none.

  -op_info PATH
    Passes file path to be used where optional info is printed.

  -out DATA, --output DATA
    Triggers output of ascii files.
    Available: none, bins, raw, errors, cuts, all.
    Default: none.

  -v VERBOSITY, --verbosity VERBOSITY
    Controls the verbosity of MINUIT.
    Default: -1.
      )";
      std::cout << help_text << std::endl;
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
      if(i < argc - 1 && !strcmp(argv[i+1], "sample")){
        args->format = Format_t::Sample;
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
      if(i < argc - 1 && !strcmp(argv[i+1], "mswmsl") && !(args->graphics & 4)){
        args->graphics += 4;
      }
      if(i < argc - 1 && !strcmp(argv[i+1], "all")){
        args->graphics = 7;
      }
    }
    if(!strcmp(argv[i], "-op") || !strcmp(argv[i], "--op-info")){
      if(i < argc -1) args->op_info = argv[i+1];
    }
  }
  return 0;
}

void prepare_std_output_files(args_t args){
  std::ofstream f1("fitstats_bb.csv");
  std::ofstream f2("fitstats_nobb.csv");
  std::ofstream f3("summary.csv");
  if(args.bin_vars & 1)  {f1 << "ZA,"; f2 << "ZA,"; f3 << "ZA,";}
  if(args.bin_vars & 2)  {f1 << "E,"; f2 << "E,"; f3 << "E,";}
  if(args.bin_vars & 4)  {f1 << "T,"; f2 << "T,"; f3 << "T,";}
  if(args.bin_vars & 8)  {f1 << "A,"; f2 << "A,"; f3 << "A,";}
  if(args.bin_vars & 16) {f1 << "O,"; f2 << "O,"; f3 << "O,";}
  f1 << "bkgfrac_nosrc, bkgfrac, srcfrac, dataCt, bkgCt, srcCt, lnL_nosrc, lnL_src, TS" << std::endl;
  f2 << "bkgfrac_nosrc, bkgfrac, srcfrac, dataCt, bkgCt, srcCt, lnL_nosrc, lnL_src, TS" << std::endl;
  f3 << "ct_ratio,TS_noBB,TS_BB,Lima_std,Lima_bb" << std::endl;
  f1.close(); f2.close(); f3.close();
  if(args.output & 8) print_cuts("reset", 0);
}

int optional_binning(indices_t indices, args_t args){
  std::stringstream path;
  int zi = 1, ei = 0, ti = 0, ai = 0, oi = 0; //TODO
  if(!(args.bin_vars & 1)){
    if(indices.za != zi) return 1;
  }
  else path << "ZA" << indices.za;
  if(!(args.bin_vars & 2)){
    if(indices.e != ei) return 1;
  }
  else path << "E" << indices.e;
  if(!(args.bin_vars & 4)){
    if(indices.tel != ti) return 1;
  }
  else path << "T" << TBINS[indices.tel];
  if(!(args.bin_vars & 8)){
    if(indices.az != ai) return 1;
  }
  else path << "A" << indices.az;
  if(!(args.bin_vars & 16)){
    if(indices.off != oi) return 1;
  }
  else path << "O" << indices.off;
  std::cout << path.str() << std::endl;
  OUTPATH = path.str();
  return 0;
}

//Output
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

  //Set up
  TCanvas *c1 = new TCanvas("","",1200,1200);
  c1->cd();
  dathist->SetLineColor(4);
  bkghist->SetLineColor(6);

  dathist->SetStats(false);
  bkghist->SetStats(false);

  dathist->SetTitle(OUTPATH.c_str());

  bkghist->Scale(DAT_HIST->Integral() / BKG_HIST->Integral());

  TLegend *legend;
  legend = new TLegend(0.12, 0.8, 0.3, 0.9);
  legend->AddEntry(dathist, "Data");
  legend->AddEntry(bkghist, "Background Template");

  dathist->SetMinimum(0);
  dathist->SetMaximum(std::max(dathist->GetMaximum(), bkghist->GetMaximum()) * 1.1);

  //Draw
  dathist->Draw();
  dathist->Draw("sameE0");
  bkghist->Draw("same");
  bkghist->Draw("sameE0");
  legend->Draw();

  //Save and clean up
  c1->SaveAs(filepath.str().c_str());
  c1->Clear();
  delete c1;
  delete dathist;
  delete bkghist;
  delete legend;
}

void histogram_fit_data(double fracs[6], indices_t ins){
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
  TRatioPlot* rp = new TRatioPlot(F1, DAT_HIST, "diff");
  rp->Draw();

  //Write Data
  c1->cd(4);
  TPaveText *pt = new TPaveText(0, 0, 1, 1);
  std::stringstream line;
  pt->AddText(.5, .95, "Standard Fit Values:");
  line << "P_b = " << fracs[0];
  pt->AddText(.25, .8, line.str().c_str());
  line.str("");
  line << "P_s = " << fracs[1];
  pt->AddText(.75, .8, line.str().c_str());
  line.str("");
  line << "TS = " << -2 * (src_noBB(fracs[0], fracs[1]) - nosrc_noBB(fracs[4]));
  pt->AddText(.5, .7, line.str().c_str());
  pt->AddLine(0, .5, 1, .5);
  pt->AddText(.5, .45, "Barlow-Beeston Fit Values:");
  line.str("");
  line << "P_b = " << fracs[2];
  pt->AddText(.25, .3, line.str().c_str());
  line.str("");
  line << "P_s = " << fracs[3];
  pt->AddText(.75, .3, line.str().c_str());
  line.str("");
  line << "TS = " << -2 * (src_BB(fracs[2], fracs[3]) - nosrc_BB(fracs[5]));
  pt->AddText(.5, .2, line.str().c_str());
  pt->SetAllWith("=", "size", .05);
  pt->Draw();

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

  gStyle->SetPalette(kBlackBody);
  map->Draw("SURF3");
  std::stringstream out_file;
  out_file << boost::replace_all_copy(title.str(), " ", "_") << ".png";
  c1->SaveAs(out_file.str().c_str());

  //clean up
  delete c1;
  delete map;
}

void plot_msw_vs_msl(TH2F* dat_2hist, TH2F* bkg_2hist){
  dat_2hist->GetXaxis()->SetTitle("MSW");
  dat_2hist->GetYaxis()->SetTitle("MSL");
  bkg_2hist->GetXaxis()->SetTitle("MSW");
  bkg_2hist->GetYaxis()->SetTitle("MSL");
  double dat, bkg;
  dat = dat_2hist->GetMaximum();
  bkg = bkg_2hist->GetMaximum();
  dat_2hist->SetMaximum(dat > bkg ? dat : bkg);
  bkg_2hist->SetMaximum(dat > bkg ? dat : bkg);
  dat = dat_2hist->GetMinimum();
  bkg = bkg_2hist->GetMinimum();
  dat_2hist->SetMinimum(dat < bkg ? dat : bkg);
  bkg_2hist->SetMinimum(dat < bkg ? dat : bkg);

  std::stringstream title;
  title << dat_2hist->GetTitle() << " " << OUTPATH;
  dat_2hist->SetTitle(title.str().c_str());
  title.str("");
  title << bkg_2hist->GetTitle() << " " << OUTPATH;
  bkg_2hist->SetTitle(title.str().c_str());

  gStyle->SetPalette(kBlackBody);
  gStyle->SetOptStat(0);

  TCanvas c1("c1", "c1", 2400, 1200);
  c1.Divide(2,1);
  c1.cd(1);
  dat_2hist->Draw("COLZ");

  c1.cd(2);
  bkg_2hist->Draw("COLZ");

  std::stringstream out_file;
  out_file << "MSWMSL_" << OUTPATH << ".png";
  c1.SaveAs(out_file.str().c_str());
}