#include "output.h"

void printRawData(hists_t *hists, std::string fit_param){
  TH1D* dathist;
  TH1D* bkghist;
  TH1D* srchist;
  int nbins;
  if(fit_param == "MSW"){
    dathist = hists->msw_dat;
    bkghist = hists->msw_bkg;
    srchist = hists->msw_src;
    nbins = NBIN;
  }
  else if(fit_param == "BDT"){
    dathist = hists->bdt_dat;
    bkghist = hists->bdt_bkg;
    srchist = hists->bdt_src;
    nbins = NBIN;
  }

  std::stringstream path;
  path << "RawData_" << fit_param << "_" << hists->outpath << ".csv";
  std::ofstream f(path.str().c_str());
  f << "Bin: " << hists->longoutpath << std::endl;
  f << fit_param << ", di, bi, si" << std::endl;
  for(int i = 1; i <= nbins; i++){
    f << dathist->GetBinCenter(i) << ","
    << dathist->GetBinContent(i) << ","
    << bkghist->GetBinContent(i) << ","
    << srchist->GetBinContent(i) << std::endl;
  }
}

void histogram_raw_data(hists_t *hists, std::string fit_param){
  std::stringstream filepath;
  filepath << "HIST_RAW_" << fit_param << "_" << hists->outpath << ".png";

  //Copy the global histograms so that we can edit them
  TH1D *dathist, *bkghist;
  if(fit_param == "MSW"){
    dathist = new TH1D(*(hists->msw_dat));
    bkghist = new TH1D(*(hists->msw_bkg));
  }
  else if(fit_param == "BDT"){
    dathist = new TH1D(*(hists->bdt_dat));
    bkghist = new TH1D(*(hists->bdt_bkg));
  }

  //Set up
  TCanvas *c1 = new TCanvas("","",1200,1200);
  c1->cd();
  dathist->SetLineColor(4);
  bkghist->SetLineColor(6);

  dathist->SetStats(false);
  bkghist->SetStats(false);

  std::stringstream title;
  title << "Raw " << hists->longoutpath;
  dathist->SetTitle(title.str().c_str());

  bkghist->Scale(dathist->Integral() / bkghist->Integral());

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

void histogram_fit_data(double fracs[6], indices_t ins, args_t *args, hists_t *hists, std::string fit_param, TFitter *fitter){
  std::stringstream filepath;
  filepath << "HIST_FIT_" << fit_param << " " << hists->outpath << ".png";
  TH1D *dathist, *bkghist, *srchist;

  //Copy the global histograms so that we can edit them
  TH1D *F1, *B1;
  if(fit_param == "MSW"){
    dathist = new TH1D(*hists->msw_dat);
    bkghist = new TH1D(*hists->msw_bkg);
    srchist = new TH1D(*hists->msw_src);
    F1 = new TH1D("F1Hist", "MSW Fit", NBIN, MSWLOW, MSWHIGH);
    B1 = new TH1D("B1Hist", "BB Fit", NBIN, MSWLOW, MSWHIGH);
  }
  else if(fit_param == "BDT"){
    dathist = new TH1D(*hists->bdt_dat);
    bkghist = new TH1D(*hists->bdt_bkg);
    srchist = new TH1D(*hists->bdt_src);
    F1 = new TH1D("F1Hist", "BDT Fit", NBIN, BDTLOW, BDTHIGH);
    B1 = new TH1D("B1Hist", "BB Fit", NBIN, BDTLOW, BDTHIGH);
  }
  double dat_int = dathist->Integral();
  double bkg_int = bkghist->Integral();
  double src_int = srchist->Integral();
  TH1D* raw_src = new TH1D(*srchist);

  //Canvas set up
  
  TCanvas *c1 = new TCanvas("", hists->outpath.c_str(), 1600, 1600);
  c1->Divide(2,2);

  dathist->SetLineColor(4);
  bkghist->SetLineColor(6);
  srchist->SetLineColor(1);
  raw_src->SetLineColor(1);
  F1->SetLineColor(8);
  B1->SetLineColor(6);

  dathist->SetStats(false);
  bkghist->SetStats(false);
  srchist->SetStats(false);
  raw_src->SetStats(false);
  F1->SetStats(false);
  B1->SetStats(false);

  TLegend *legend0 = new TLegend(0.15, 0.6, 0.45, 0.9);
  TLegend *legend1 = new TLegend(0.15, 0.6, 0.3, 0.9);
  TLegend *legend2 = new TLegend(0.15, 0.7, 0.3, 0.9);

  //Draw Raw Data
  c1->cd(1);
  legend0->AddEntry(dathist, "Data");
  legend0->AddEntry(bkghist, "Background Template");
  legend0->AddEntry(raw_src, "Source Template");

  bkghist->Scale(dat_int / bkg_int);
  raw_src->Scale(dat_int / src_int);

  dathist->SetMinimum(0);
  dathist->SetMaximum(std::max(dathist->GetMaximum(), bkghist->GetMaximum()) * 1.1);

  dathist->Draw();
  dathist->Draw("sameE0");
  bkghist->Draw("same");
  bkghist->Draw("sameE0");
  raw_src->Draw("sameE0");
  legend0->Draw();

  //Draw BB Src
  c1->cd(2);
  src_BB(fracs[2], fracs[3], false, F1, B1);
  srchist->Scale(fracs[3]);
  srchist->Scale(dat_int / src_int);
  B1->Scale(fracs[2]);
  B1->Scale(dat_int / bkg_int);
  F1->Draw();
  F1->Draw("sameE0");
  B1->Draw("sameE0");
  srchist->Draw("sameE0");

  F1->SetMinimum(0);
  F1->SetMaximum(std::max(F1->GetMaximum(), B1->GetMaximum()) *1.1);

  legend1->AddEntry(F1, "fi");
  legend1->AddEntry(B1, "Bi");
  legend1->AddEntry(srchist, "Si");
  legend1->Draw();

  //Draw Residual
  c1->cd(3);
  TH1D* F1copy = new TH1D(*F1);
  F1copy->SetTitle("Residual");
  F1copy->SetMaximum(std::max(F1copy->GetMaximum(), dathist->GetMaximum()) *1.1);
  TRatioPlot_BetterError* rp = new TRatioPlot_BetterError(F1copy, dathist, "diffsig");
  rp->SetH1DrawOpt("E0");
  rp->SetH2DrawOpt("E0");
  rp->Draw();
  legend2->AddEntry(F1copy, "Fit");
  legend2->AddEntry(dathist, "Raw");
  legend2->Draw();

  //Get Errors from Fitters
  double a, b, c;
  double err_pb, err_ps;
  fitter->GetErrors(0, a, b, err_pb, c);
  fitter->GetErrors(1, a, b, err_ps, c);

  //Write Data
  c1->cd(4);
  TPaveText *pt = new TPaveText(0, 0, 1, 1);
  std::stringstream line;
  pt->AddText(.35, .95, "Barlow-Beeston:")->SetTextAlign(12);
  line.str("");
  line << std::scientific << std::setprecision(3);
  line << "P_b = " << fracs[2] << " #pm " << err_pb;
  pt->AddText(.1, .85, line.str().c_str())->SetTextAlign(12);
  line.str("");
  line << "P_s = " << fracs[3] << " #pm " << err_ps;
  pt->AddText(.1, .8, line.str().c_str())->SetTextAlign(12);
  line.str("");
  line << "TS = " << -2 * (src_BB(fracs[2], fracs[3]) - nosrc_BB(fracs[5]));
  pt->AddText(.1, .75, line.str().c_str())->SetTextAlign(12);
  line.str("");
  line << std::defaultfloat;
  if(args->bin_vars & 1){
    line << "ZA: " << ZABINS[ins.za] << "-" << ZABINS[ins.za+1];
    pt->AddText(.05, .65, line.str().c_str())->SetTextAlign(12);
    line.str("");
  }
  if(args->bin_vars & 2){
    line << "E: " << EBINS[ins.e] << "-" << EBINS[ins.e+1];
    pt->AddText(.45, .65, line.str().c_str())->SetTextAlign(12);
    line.str("");
  }
  if(args->bin_vars & 4){
    line << "Tel: " << TBINS[ins.tel];
    pt->AddText(.85, .65, line.str().c_str())->SetTextAlign(12);
    line.str("");
  }
  if(args->bin_vars & 8){
    line << "Az: " << AZBINS[ins.az] << "-" << AZBINS[ins.az+1];
    pt->AddText(.05, .6, line.str().c_str())->SetTextAlign(12);
    line.str("");
  }
  if(args->bin_vars & 16){
    line << "Off: " << OBINS[ins.off] << "-" << OBINS[ins.off+1];
    pt->AddText(.45, .6, line.str().c_str())->SetTextAlign(12);
    line.str("");
  }

  pt->AddLine(0, .55, 1, .55);
  pt->SetAllWith("=", "size", .05);
  pt->SetAllWith(": ", "size", .05);
  if(args->op_info.c_str()){
    std::ifstream infile(args->op_info);
    std::string readline;
    std::getline(infile, readline);
    pt->AddText(.05, .5, readline.c_str())->SetTextAlign(12);
    std::getline(infile, readline);
    pt->AddText(.05, .4, readline.c_str())->SetTextAlign(12);
    std::getline(infile, readline);
    pt->AddText(.05, .3, readline.c_str())->SetTextAlign(12);
    std::getline(infile, readline);
    pt->AddText(.05, .2, readline.c_str())->SetTextAlign(12);
    std::getline(infile, readline);
    pt->AddText(.05, .1, readline.c_str())->SetTextAlign(12);
    infile.close();
  }
  pt->Draw();

  //Save and Clean up
  c1->SaveAs(filepath.str().c_str());
  c1->Clear();
  delete c1;
  delete dathist;
  delete bkghist;
  delete srchist;
  delete F1;
  delete B1;
  delete legend0;
  delete legend1;
}

void print_errors(TFitter *fitter, std::string fit_param, std::string pathbase){
  std::stringstream fname;
  fname << "Errors_" << fit_param << ".csv";
  if(access(fname.str().c_str(), F_OK)){
    std::ofstream f(fname.str().c_str());
    f << PRINTSPACE << "Bin"
      << PRINTSPACE << "Parameter"
      << PRINTSPACE << "eplus"
      << PRINTSPACE << "eminus"
      << PRINTSPACE << "eparab"
      << PRINTSPACE << "globcc"
      << PRINTSPACE << "Parameter"
      << PRINTSPACE << "eplus"
      << PRINTSPACE << "eminus"
      << PRINTSPACE << "eparab"
      << PRINTSPACE << "globcc" << std::endl;
    f.close();
  }
  double eplus, eminus, eparab, globcc;
  fitter->GetErrors(0, eplus, eminus, eparab, globcc);

  std::ofstream f(fname.str().c_str(), std::ios::out | std::ios::app);
  f << PRINTSPACE << pathbase
    << PRINTSPACE << "bkgfrac"
    << PRINTSPACE << eplus
    << PRINTSPACE << eminus
    << PRINTSPACE << eparab
    << PRINTSPACE << globcc;

  fitter->GetErrors(1, eplus, eminus, eparab, globcc);
  f << PRINTSPACE << "srcfrac"
    << PRINTSPACE << eplus
    << PRINTSPACE << eminus
    << PRINTSPACE << eparab
    << PRINTSPACE << globcc << std::endl;
  f.close();
}

void print_cuts(std::string pathbase, cuts_t* cuts, std::string outpath){
  std::stringstream fname;
  fname << "Cuts_" << pathbase << ".txt";
  if(access(fname.str().c_str(), F_OK)){
    std::ofstream f(fname.str().c_str());
    f << PRINTSPACE << "Bin"
      << PRINTSPACE << "passed"
      << PRINTSPACE << "read"
      << PRINTSPACE << "source"
      << PRINTSPACE << "telescope"
      << PRINTSPACE << "energy"
      << PRINTSPACE << "zenith angle"
      << PRINTSPACE << "azimuth"
      << PRINTSPACE << "offset" << std::endl;
    f.close();
  }
  std::ofstream f(fname.str().c_str(), std::ios::out | std::ios::app);
  f << PRINTSPACE << outpath
    << PRINTSPACE << cuts->passed
    << PRINTSPACE << cuts->read
    << PRINTSPACE << cuts->src
    << PRINTSPACE << cuts->tel
    << PRINTSPACE << cuts->e
    << PRINTSPACE << cuts->za
    << PRINTSPACE << cuts->az
    << PRINTSPACE << cuts->off << std::endl;
  f.close();
}

void map_likelihood(double Pb, double Ps, std::string title_tag, indices_t ins, args_t args, std::string outpath, std::string longoutpath){
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
  title << "Likelihood " << title_tag << " " << longoutpath;

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1", 1200,800);
  TH2D *map = new TH2D("map", title.str().c_str(), xbins, xmin, xmax, ybins, ymin, ymax);
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
  out_file << "Likelihood_" << title_tag << "_" << outpath << ".png";
  c1->SaveAs(out_file.str().c_str());

  //clean up
  delete c1;
  delete map;
}

void plot_msw_vs_msl(hists_t *hists){
  hists->msw_msl_dat->GetXaxis()->SetTitle("MSW");
  hists->msw_msl_dat->GetYaxis()->SetTitle("MSL");
  hists->msw_msl_bkg->GetXaxis()->SetTitle("MSW");
  hists->msw_msl_bkg->GetYaxis()->SetTitle("MSL");
  double dat, bkg;
  dat = hists->msw_msl_dat->GetMaximum();
  bkg = hists->msw_msl_bkg->GetMaximum();
  hists->msw_msl_dat->SetMaximum(dat > bkg ? dat : bkg);
  hists->msw_msl_bkg->SetMaximum(dat > bkg ? dat : bkg);
  dat = hists->msw_msl_dat->GetMinimum();
  bkg = hists->msw_msl_bkg->GetMinimum();
  hists->msw_msl_dat->SetMinimum(dat < bkg ? dat : bkg);
  hists->msw_msl_bkg->SetMinimum(dat < bkg ? dat : bkg);

  std::stringstream title;
  title << hists->msw_msl_dat->GetTitle() << " " << hists->longoutpath;
  hists->msw_msl_dat->SetTitle(title.str().c_str());
  title.str("");
  title << hists->msw_msl_bkg->GetTitle() << " " << hists->longoutpath;
  hists->msw_msl_bkg->SetTitle(title.str().c_str());

  gStyle->SetPalette(kBlackBody);
  gStyle->SetOptStat(0);

  TCanvas c1("c1", "c1", 2400, 1200);
  c1.Divide(2,1);
  c1.cd(1);
  hists->msw_msl_dat->Draw("COLZ");

  c1.cd(2);
  hists->msw_msl_bkg->Draw("COLZ");

  std::stringstream out_file;
  out_file << "MSWMSL_" << hists->outpath << ".png";
  c1.SaveAs(out_file.str().c_str());
}

void init_output_root_file(){
  TFile *fout = new TFile("output.root", "NEW");
  if(!fout->IsOpen()){
    std::cerr << "Failed to init output file." << std::endl;
  }
  TDirectory *dir_src_excl = fout->mkdir("Source Excl");
  dir_src_excl->mkdir("MSW");
  dir_src_excl->mkdir("BDT");
  TDirectory *dir_all = fout->mkdir("All");
  dir_all->mkdir("MSW");
  dir_all->mkdir("BDT");
  fout->Write();
  delete fout;
}


/**
  Outputs a bin's data to the output root file.
  @param args Current arguments struct.
  @param hists Current hists struct.
  @param fracs Pointer to start of array of
  background and source fractions. The order
  should be:
    bkgfrac_src_nobb
    srcfrac_src_nobb
    bkgfrac_src_bb
    srcfrac_src_bb
  @param fit_param Which parameter is being fit on.
  Should have the value "MSW" or "BDT".
  @param srcexcl Whether sources were excluded.
*/
void write_to_root_file(args_t *args, hists_t *hists, double *fracs, std::string fit_param, bool srcexcl){
  TFile *fout = new TFile("output.root", "UPDATE");
  if (!fout->IsOpen()){
    std::cerr << "Failed to open output file." << std::endl;
  }
  TDirectory *tld;
  if(srcexcl) tld = fout->Get("Source Excl");
  else tld = fout->Get("All");
  TDirectory *paramdir = tld->Get(fit_param);
  TDirectory *maindir = paramdir->mkdir(hists->longoutpath);

  maindir->cd();
  //Write input histograms
  double parlow, parhigh;
  if(fit_param == "MSW"){
    hists->msw_dat->Write("Data");
    hists->msw_bkg->Write("Background");
    hists->msw_src->Write("Source");
    parlow = MSWLOW; parhigh = MSWHIGH;
  }
  else if(fit_param == "BDT"){
    hists->bdt_dat->Write("Data");
    hists->bdt_bkg->Write("Background");
    hists->bdt_src->Write("Source");
    parlow = BDTLOW; parhigh = BDTHIGH;
  }
  TH1D* nobb_F = new TH1D("nobb_F", "Std Fit", NBINS, parlow, parhigh);
  TH1D* bb_F = new TH1D("bb_F", "BB Fit Data", NBINS, parlow, parhigh);
  TH1D* bb_B = new TH1D("bb_B", "BB Fit Bkg", NBINS, parlow, parhigh);

  src_noBB(fracs[0], fracs[1], 0, nobb_F);
  src_BB(fracs[2], fracs[3], 0, bb_F, bb_B);

  nobb_F->Write();
  bb_F->Write();
  bb_B->Write();

  delete fout;
}
