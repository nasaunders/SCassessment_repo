// simple assessment model implemented in TMB

#include <TMB.hpp>

// template for square function

template <class Type> Type square(Type x){return x*x;}

// main function that defines the objective function, reports etc.

template<class Type>
Type objective_function<Type>::operator() ()
{
  // data, fixed quantities, estimated parameters

  DATA_IVECTOR(dms);            // Number of strata and stations.
  DATA_VECTOR(areaStr);         // Surface area by strata.
  DATA_VECTOR(dCaught);         // Density caught by station.
  DATA_VECTOR(meanWght);        // Mean weight of individuals by station.
  DATA_IVECTOR(StrSt);          // Strata by station.


  PARAMETER_VECTOR(logBstr);
  PARAMETER_VECTOR(logStdBstr);

  // relevant dimensions

  int nStr = dms[1];
  int nsta = dms[2];

  // define the objective function variable and NLL components

  Type f = 0;

  // now define all the relevant process variables we need

  vector<Type> Btot(1);         // Total biomass
  vector<Type> Bdens(nsta);     // Biomass density by stratum.
  vector<Type> BstrObs(nStr);   // Biomass by stratum.
  // vector<Type> StdBstr(nStr);   // Standard deviation of biomass by stratum.

  for (int i=0;i<nsta;i++){
    Bdens(i) = dCaught(i) * meanWght(i);
  }

  Type n = 0;
  for (int st =0;st<nStr;st++){
    BstrObs(st) = 0.0;
    n = 0;
    Btot(1) = 0.0;
    for (int i=0;i<nsta;i++){
      if (StrSt[i] == (st + 1)){
        BstrObs(st) += Bdens(i) ;
        n++ ;
      }
    }
    BstrObs(st) /= double(n);
    BstrObs(st) *= areaStr(st);
    Btot(1) += BstrObs(st);
  }



  //////////////////////////////////
  // define likelihood components //
  //////////////////////////////////

  //

  Type StdBstr = exp(logStdBstr);
  Type Bstr = exp(logBstr);

  for(int st=0;st<nStr;st++)
    f -= dnorm(BstrObs(st), Bstr(st), StdBstr(st), true);

  // total objective function


  ///////////////////////
  // Reporting section //
  ///////////////////////

  // ADREPORT (things for which we want std. errors)

  ADREPORT(Btot);
  ADREPORT(BstrObs);

  // REPORT (things we just wanted predicted MLEs of)

  REPORT(f);

  return f;
}
