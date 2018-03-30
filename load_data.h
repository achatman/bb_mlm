#ifndef LOAD_DATA_H
#define LOAD_DATA_H
#include "bbts.h"

TH1F** loadData(indices_t ins, args_t args, double* alpha);
TH1F** loadData_csv(indices_t ins, args_t args, double* alpha);
TH1F** loadData_vegas(indices_t ins, args_t args, double* alpha);
TH1F** loadData_sample(indices_t ins, args_t args, double* alpha);

#endif