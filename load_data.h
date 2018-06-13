#ifndef LOAD_DATA_H
#define LOAD_DATA_H
#include "bbts.h"

void loadData(indices_t ins, args_t args, double *alpha, TH1F* DAT_HIST, TH1F* BKG_HIST, TH1F* SRC_HIST, TH2F* DAT_2HIST=0, TH2F* BKG_2HIST=0);

#endif