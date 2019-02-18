#include "output.h"

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

std::string get_outpath(indices_t *ins, bool more){
    std::stringstream outpath;
    if(ins->za != -1){
      outpath << "ZA" << ZABINS[ins->za] << "-" << ZABINS[ins->za+1];
    }
    if(ins->e != -1){
      if(outpath.str().length()) outpath << "_";
      outpath << "E" << EBINS[ins->e] << "-" << EBINS[ins->e+1];
    }
    if(ins->tel != -1){
      if(outpath.str().length()) outpath << "_";
      outpath << "T" << TBINS[ins->tel];
    }
    if(ins->az != -1){
      if(outpath.str().length()) outpath << "_";
      outpath << "AZ" << AZBINS[ins->az] << "-" << AZBINS[ins->az+1];
    }
    if(ins->off != -1){
      if(outpath.str().length()) outpath << "_";
      outpath << "O" << OBINS[ins->off] << "-" << OBINS[ins->off+1];
    }
    if(more){
        outpath << "_SE" << ins->src_excl;
    }
    return outpath.str();
}



/**
  Outputs a bin's data to the output root file.
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
void write_to_root_file(indices_t *ins, hists_t *hists, double *fracs, Fit_Par_t fit_param){
  TFile *fout = new TFile("output.root", "UPDATE");
  if (!fout->IsOpen()){
    std::cerr << "Failed to open output file." << std::endl;
  }
  TDirectory *tld;
  if(ins->src_excl) tld = (TDirectory*)fout->Get("Source Excl");
  else tld = (TDirectory*)fout->Get("All");
  TDirectory *paramdir;
  if(fit_param == Fit_Par_t::msw) (TDirectory*)tld->Get("MSW");
  else if(fit_param == Fit_Par_t::bdt) (TDirectory*)tld->Get("BDT");
  TDirectory *maindir = paramdir->mkdir(get_outpath(ins).c_str());

  maindir->cd();
  //Write input histograms
  hists->dat->Write("Data");
  hists->bkg->Write("Background");
  hists->src->Write("Source");
  double parlow, parhigh;
  if(fit_param == Fit_Par_t::msw){
    parlow = MSWLOW; parhigh = MSWHIGH;
  }
  else if(fit_param == Fit_Par_t::bdt){
    parlow = BDTLOW; parhigh = BDTHIGH;
  }
  TH1D* nobb_F = new TH1D("nobb_F", "Std Fit", NBIN, parlow, parhigh);
  TH1D* bb_F = new TH1D("bb_F", "BB Fit Data", NBIN, parlow, parhigh);
  TH1D* bb_B = new TH1D("bb_B", "BB Fit Bkg", NBIN, parlow, parhigh);

  src_noBB(fracs[0], fracs[1], nobb_F);
  src_BB(fracs[2], fracs[3], bb_F, bb_B);

  nobb_F->Write();
  bb_F->Write();
  bb_B->Write();

  delete fout;
}
