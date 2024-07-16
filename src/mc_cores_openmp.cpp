
#include "common.h"

#ifndef NoRcpp

RcppExport SEXP mc_cores_openmp(void) {

#else

int mc_cores_openmp(void) {

#endif

#ifdef _OPENMP

int mc_cores_openmp=omp_get_num_threads();

#else

int mc_cores_openmp=0;

#endif

#ifndef NoRcpp

return Rcpp::wrap(mc_cores_openmp);

#else

return mc_cores_openmp;

#endif

}

