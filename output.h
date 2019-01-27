#ifndef OUTPUT_H
#define OUTPUT_H
#include "bbts.h"

void init_output_root_file();
void write_to_root_file(args_t *args, hists_t *hists, double *fracs, Fit_Par_t fit_param, bool srcexcl);

#endif
