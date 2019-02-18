#ifndef OUTPUT_H
#define OUTPUT_H
#include "bbts.h"

std::string get_outpath(indices_t *ins, bool more=false);
void init_output_root_file();
void write_to_root_file(indices_t *ins, hists_t *hists, double *fracs, Fit_Par_t fit_param);
#endif
