#define SIM_UNSTABLE 1
#define SIM_SUCCESS 0

#include <stdlib.h>
#include <string.h>
#include <math.h>

int calc_simulation(const double* Dvector,
                    const double* Rvector,
                    int nIV,
                    const double* initCond,
                    int ndt, int ndx,
                    double dt, double dx,
                    double r,
                    double* simResults)
{
    int simStorageBytes = sizeof(double) * ndx;
    //this stores step i;
    double* prevStep = (double*)malloc(simStorageBytes);
    memcpy(prevStep, initCond, simStorageBytes);

    //work out our constant composition boundary conditions
    double leftBoundary = initCond[0];
    double rightBoundary = initCond[ndx-1];

    //loop in time
    for (int n = 0; n < ndt; n++)
    {
        //apply boundary conditions
        simResults[0] = leftBoundary;
        simResults[ndx-1] = rightBoundary;

        //loop over space
        for (int k = 1; k < ndx - 1; k++)
        {
            //Work out lookup indicies for D and R based on concentration
            int Ckm1Index = lround(prevStep[k-1] * (nIV - 1));
            int CkIndex = lround(prevStep[k] * (nIV - 1));
            int Ckp1Index = lround(prevStep[k+1] * (nIV - 1));

            //This tells us if our previous concentrations were outside [0,1]
            //and hence if we're about to dereference garbage
            if (Ckm1Index < 0 || CkIndex < 0 || Ckp1Index < 0 ||
                Ckm1Index >= nIV || CkIndex >= nIV || Ckp1Index >= nIV)
            {
                free(prevStep);
                return SIM_UNSTABLE;
            }

            //Indicies for calculating either central, left or right difference derivative of dDdc/dRdc;
            //Central unless we have C = 0 or C = 1
            int dLeftIndex = CkIndex == 0 ? CkIndex : CkIndex - 1;
            int dRightIndex = CkIndex == (nIV - 1) ? CkIndex : CkIndex + 1;

            //Calculate derivatives based on these indicies
            //Multiply by nIV-1 to turn indicies into units of concentration
            double dDdc = (Dvector[dRightIndex] - Dvector[dLeftIndex]) /
                        (dRightIndex - dLeftIndex) * (nIV - 1);
            double dRdc = (Rvector[dRightIndex] - Rvector[dLeftIndex]) /
                        (dRightIndex - dLeftIndex) * (nIV - 1);

            //Current values of D and R
            double Dv = Dvector[CkIndex];
            double Rv = Rvector[CkIndex];

            //Concentrations and derivatives thereof
            double C = prevStep[k];
            double dCdx = (prevStep[k+1] - prevStep[k]) / dx;
            double d2Cdx2 = (prevStep[k+1] - 2 * prevStep[k] + prevStep[k-1]) / (dx * dx);

            //the next step
            /*double simc = dt * (Dv*d2Cdx2 + dDdc*dCdx*dCdx - dDdc*dCdx*C*Rv*r - Dv*Rv*dCdx*r
                                  -Dv*dRdc*dCdx*C*r) + prevStep[k];
                                  */
            double simc = dt * (Dv*d2Cdx2 + dDdc*dCdx*dCdx - dDdc*C*Rv*r - Dv*Rv*dCdx*r
                                   -Dv*dRdc*dCdx*r) + prevStep[k];                   
            //Smooth out major fuckups in simc; for somereason these are appearing
            //Probably because dDdC is much bigger with this model
            /*if (simc > 1.0)
                simc = 1.0;
            else if (simc < 0.0)
                simc = 0.0;*/
            simResults[k] = simc;
        }

        //Now swap current step into prev step
        memcpy(prevStep, simResults, simStorageBytes);
    }

    //The result is now in memory (simResults); free and return
    free(prevStep);
    return SIM_SUCCESS;
}

