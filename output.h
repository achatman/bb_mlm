#ifndef OUTPUT_H
#define OUTPUT_H
#include "bbts.h"

void printRawData(hists_t *hists, std::string fit_param);
void histogram_raw_data(hists_t *hists, std::string fit_param);
void histogram_fit_data(double fracs[6], indices_t ins, args_t *args, hists_t *hists, std::string fit_param, TFitter *fitter);
void print_errors(TFitter *fitter, std::string fit_param, std::string pathbase);
void print_cuts(std::string action, cuts_t* cuts, std::string outpath = "");
void map_likelihood(double Pb, double Ps, std::string title_tag, indices_t ins, args_t args, std::string outpath, std::string longoutpath);
void plot_msw_vs_msl(hists_t *hists);


void init_output_root_file();
void write_to_root_file(args_t *args, hists_t *hists, double *fracs, std::string fit_param, bool srcexcl);

#endif
