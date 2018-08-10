#ifndef OUTPUT_H
#define OUTPUT_H
#include "bbts.h"

void printRawData(hists_t *hists, std::string fit_param);
void histogram_raw_data(hists_t *hists, std::string fit_param);
void histogram_fit_data(double fracs[6], indices_t ins, args_t *args, hists_t hists);
void calculate_errors(double Pb_m, double Ps_m, double sigma_Pb_m, double sigma_Ps_m, indices_t ins, double alpha, hists_t hists);
void print_cuts(std::string action, cuts_t* cuts, std::string outpath = "");
void map_likelihood(double Pb, double Ps, std::string title_tag, indices_t ins, args_t args, std::string outpath, std::string longoutpath);
void plot_msw_vs_msl(hists_t *hists);

#endif
