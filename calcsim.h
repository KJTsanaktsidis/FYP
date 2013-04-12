int calc_simulation(const double* Dvector,
                    const double* Rvector,
                    int nIV,
                    const double* initCond,
                    int ndt, int ndx,
                    int dt, int dx,
                    double r,
                    double* simResults);

#define SIM_SUCCESS 0
#define SIM_UNSTABLE 1
