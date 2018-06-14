#ifndef LOAD_DATA_H
#define LOAD_DATA_H
#include "bbts.h"

void loadData(indices_t ins, args_t args, double *alpha, TH1D* DAT_HIST, TH1D* BKG_HIST, TH1D* SRC_HIST, TH2D* DAT_2HIST=0, TH2D* BKG_2HIST=0);

#endif