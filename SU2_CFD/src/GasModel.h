/*
 * GasModel_.h
 *
 *     Created on: Nov 17, 2011
 *  Last modified: Dec 21, 2015
 *         Author: Enrico Rinaldi
 *          email: enri.rinaldi@gmail.com
 */


#ifndef GASMODEL_H_
#define GASMODEL_H_


#include "../UgpWithCvCompFlow.h"
#include <deque>

#if WITH_FP
#include "fluidprop.h"
#define LEN_COMPONENTS 32


/*
 * REAL GAS - DIRECT CALLS
 *
 * Overload the thermodynamic functions contained in UgpWithCvCompFlow.h with direct calls to
 * FluidProp. See FluidProp folder for details.
 *
 * When using real gas model, the following lines must be added in Joe.in:
 * FPLIBNAME (char) name of the library. Options are: RefProp (multiparameter), StanMix (iPRSV type).  For other options see FluidProp folder.
 * FLUIDNAME (char) name of the fluid (e.g., CO2, Toluene, Water, Air, ...)
 * NCOMP (int) number of components (just 2 components implemented) (optional, dafault 1)
 * CONC1 (double) concentration of the first component (optional, default 1.0)
 */

class GasModel: virtual public UgpWithCvCompFlow {
public:
  // member variables
  string libraryName;
  string fluidName;

public:
  // constructors
  GasModel()
  {
    if (mpi_rank == 0)
      cout << "GasModel()" << endl;

    // ----------------------------------------------------------------------------------------
    // fluid prop interface
    // ----------------------------------------------------------------------------------------
    init_fluidprop();

    libraryName = getStringParam("FPLIBNAME", "RefProp");

    int  version[4];
    if (fluidprop_getversion(libraryName.c_str(), version))
      if (mpi_rank == 0)
        printf( "%s version %d.%d.%d.%d available\n", libraryName.c_str(), version[0], version[1], version[2], version[3] );
    else
    {
      if (mpi_rank == 0)
        printf( "%s library not detected\n", libraryName.c_str());
//      throw(-321);
    }

    char Comp[20][LEN_COMPONENTS];
    double Conc[20];

    // Implemented for one component only!
    fluidName = getStringParam("FLUIDNAME", "N2");
    strcpy (Comp[0], fluidName.c_str());
    Conc[0] = getDoubleParam("CONC1",1.0);
    Conc[1] = 1.0 - Conc[0];
    int ncomp = getIntParam("NCOMP", 1);

    fluidprop_setfluid(libraryName.c_str(), ncomp, Comp[0], LEN_COMPONENTS, Conc);

    R_gas   = 8.314472/fluidprop_mmol();
    rho_ref = fluidprop_density("PT",  p_ref/1.0e5, T_ref-273.15);
    GAMMA   = fluidprop_heatcapp("PT", p_ref/1.0e5, T_ref-273.15)/fluidprop_heatcapv("PT", p_ref/1.0e5, T_ref-273.15);

    if (viscMode == VISC_REFPROP)
      mu_ref = fluidprop_viscosity("PT", p_ref/1.0e5, T_ref-273.15);

    if (mpi_rank == 0)
    {
      cout << endl;
      cout << "##############################################################" << endl;
      cout << "##############################################################" << endl;
      cout << "UPDATED gas properties from FLUIDPROP" << endl;
      cout << "FLUID:               : " << fluidName << endl;
      cout << "    Pcr              : " << fluidprop_pcrit() << endl;
      cout << "    Tcr              : " << fluidprop_tcrit() << endl;
      cout << "    rhocr            : " << fluidprop_density("PT", fluidprop_pcrit(), fluidprop_tcrit()) << endl;
      cout << "-----------------------" << endl;
      cout << "    R_GAS            : " << R_gas << endl;
      cout << "    GAMMA            : " << GAMMA << endl;
      cout << "    P_REF            : " << p_ref << endl;
      cout << "    RHO_REF          : " << rho_ref << endl;
      cout << "    T_REF            : " << T_ref << endl;
      cout << "    SOS_REF          : " << fluidprop_soundspeed("PT", p_ref/1.0e5, T_ref-273.15) << endl;
      cout << "Material properties    " << endl;
      cout << "    MU_MODE          : " << viscMode << endl;
      cout << "      mu_ref         : " << mu_ref << endl;
    }

  }

  virtual ~GasModel()  {  }

public:
  /*
   * fluidprop_allprops_ calulates the complete thermodynamic state given two inputs (e.g., "Pv", "PT", "vu", etc.)
   */
  virtual void fluidprop_allprops_( const char* InputSpec, double Input1, double Input2, struct fluidstatejoe_t* OutputState, const char *props="ALL")//, double* x, double* y)
  {
//    printf("should not be here - fluidprop_allprops_ \n");
//    cout << InputSpec<< endl; throw(-1);
    fluidprop_allpropsjoe(InputSpec, Input1, Input2, OutputState);
  }

  virtual double entropyPv(const double P, const double v)
  {
    return fluidprop_entropy("Pv", P/1.0e5, v)*1.0e3;
  }

  virtual int fluidPhase(const double press, const double specVol)
  {
    double qual = fluidprop_vaporqual("Pv", press/1.0e5, specVol);

    if ((qual<1.0) && (qual>=0.0)) return(2);
    else                           return(1);
  }


  virtual double internalEnergyPv(const double press, const double specVol)
  {
//    printf("should not be here - internalEnergyPv\n"); throw(-1);

    return fluidprop_intenergy("Pv", press/1.0e5, specVol)*1.0e3;
  }

  virtual double internalEnergyPT(const double press, const double temp)
  {
//    printf("should not be here - internalEnergyPT\n"); throw(-1);

    return fluidprop_intenergy("PT", press/1.0e5, temp-273.15)*1.0e3;
  }

  /*
   * Calculates enthalpy and the derivatives of pressure with respect to density and total energy.
   * Pressure is a function of the conservative variables, i.e., P = P(rho, rhou, rhov. rhow, rhoet).
   * The ones used here are the correct expression of such derivatives, the one reported in the paper
   * from Colonna & Rebay (2004) are wrong.
   */
  virtual void calcEnthAndPressDeriv(double &hs, double &dpdr, double &dpdret,
      const double press, const double rho, const double velMagSq)
  {
//    printf("should not be here - calcEnthAndPressDeriv\n"); throw(-1);

    fluidstatejoe_t statevec;
    fluidprop_allprops_("Pv", press/1.0e5, 1.0/rho, &statevec);

    hs     = statevec.h*1000.0;
    dpdr   = statevec.alpha - (hs - press/rho)/rho*statevec.beta + 0.5*statevec.beta/rho*velMagSq;
    dpdret = statevec.beta/rho;
  }

  virtual double calcPressAndEnthVu(const double v, const double u, double &P, double &h)
  {
    P = fluidprop_pressure("vu", v, u/1.0e3)*1.0e5;
    h = fluidprop_enthalpy("vu", v, u/1.0e3)*1.0e3;
  }

  virtual void enthalpyEntropyPT(double &h, double &s, const double P, const double T)
  {
//    printf("should not be here - enthalpyEntropyPT\n"); throw(-1);

    double dummy;
    fluidstatejoe_t fp;
    fluidprop_allpropsjoe("PT", P/1.0e5, T-273.15, &fp);

    h = fp.h*1000.0;
    s = fp.s*1000.0;
  }

  virtual void calcLamMuCp(const double temp, const double rho, double &lambda, double &mul, double &cp)
  {
    switch (viscMode)
    {
    case VISC_SUTHERLAND:
      mul = mu_ref*pow(temp/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp + SL_Sref);
      lambda = GAMMA*R_gas/(GAMMA-1.0)*mul/Pr;
      cp = GAMMA / (GAMMA - 1.0) * R_gas;
      break;
    case VISC_POWERLAW:
      mul = mu_ref*pow(temp/T_ref, mu_power_law);
      lambda = GAMMA*R_gas/(GAMMA-1.0)*mul/Pr;
      cp = GAMMA / (GAMMA - 1.0) * R_gas;
      break;
    case VISC_REFPROP:
      fluidstatejoe_t state;
      fluidprop_allprops_("Tv", temp-273.15, 1.0/rho, &state);
      mul = state.eta;
      lambda = state.lambda;
      cp = state.cp*1.0e3;
      break;
    }
  }

  virtual void allPropsHS(thermoState_t &state, const double h, const double s)
  {
//    printf("should not be here - allPropsHS\n"); throw(-1);

    fluidstatejoe_t fp;
    fluidprop_allpropsjoe("hs", h*1.0e-3, s*1.0e-3, &fp);

    state.cp = fp.cp*1000.0;
    state.h = h;
    state.s = s;
    state.T = fp.T+273.15;
    state.P = fp.P*1.0e5;
    state.d = fp.d;

    switch (viscMode)
    {
    case VISC_SUTHERLAND:
      state.eta = mu_ref*pow(state.T/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(state.T + SL_Sref);
      state.lambda = GAMMA*R_gas/(GAMMA-1.0)*state.eta/Pr;
      break;
    case VISC_POWERLAW:
      state.eta = mu_ref*pow(state.T/T_ref, mu_power_law);
      state.lambda = GAMMA*R_gas/(GAMMA-1.0)*state.eta/Pr;
      break;
    case VISC_REFPROP:
      state.eta = fp.eta;
      state.lambda = fp.lambda;
      break;
    }
  }

  virtual double enthalpyPs(const double P, const double s)
  {
    return fluidprop_enthalpy("Ps", P/1.0e5, s/1.0e3)*1.0e3;
  }


  virtual void allPropsPT(thermoState_t &state, const double P, const double T)
  {
//    printf("should not be here - allPropsPT\n"); throw(-1);

    //fluidstatejoe_t fp;
    //fluidprop_allpropsjoe( "PT", P/1.0e5, T-273.15, &fp);

    fluidstatejoe_t fp;
    fluidprop_allpropsjoe( "PT", P/1.0e5, T-273.15, &fp);
//    fluidprop_allprops_( "PT", P/1.0e5, T-273.15, &fp);
    state.T = T;
    state.P = P;
    state.h = fp.h*1.0e3;
    state.cv = fp.cv*1.0e3;
    state.cp = fp.cp*1.0e3;
    state.u = fp.u*1.0e3;
    state.d = fp.d;
    state.v = fp.v;
    state.s = fp.s*1.0e3;
    state.alpha = fp.alpha;
    state.beta = fp.beta;
    state.zeta = fp.zeta;
    state.c = fp.c;

    switch (viscMode)
    {
    case VISC_SUTHERLAND:
      state.eta = mu_ref*pow(state.T/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(state.T + SL_Sref);
      state.lambda = GAMMA*R_gas/(GAMMA-1.0)*state.eta/Pr;
      break;
    case VISC_POWERLAW:
      state.eta = mu_ref*pow(state.T/T_ref, mu_power_law);
      state.lambda = GAMMA*R_gas/(GAMMA-1.0)*state.eta/Pr;
      break;
    case VISC_REFPROP:
      state.eta = fp.eta;
      state.lambda = fp.lambda;
      break;
    }
  }

  virtual void allPropsPv(thermoState_t &state, const double P, const double v)
  {
//    printf("should not be here - allPropsPv\n"); throw(-1);

    fluidstatejoe_t fp;
    fluidprop_allprops_( "Pv", P/1.0e5, v, &fp);

    state.T = fp.T+273.15;
    state.P = P;
    state.h = fp.h*1.0e3;
    state.d = 1.0/v;
    state.v = v;
    state.s = fp.s*1.0e3;
    state.u = fp.u*1.0e3;
    state.cp = fp.cp*1000.0;
    state.alpha = fp.alpha;
    state.beta = fp.beta;
    state.zeta = fp.zeta;
    state.cv = fp.cv*1.0e3;
    state.c = fp.c;


    switch (viscMode)
    {
    case VISC_SUTHERLAND:
      state.eta = mu_ref*pow(state.T/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(state.T + SL_Sref);
      state.lambda = GAMMA*R_gas/(GAMMA-1.0)*state.eta/Pr;
      break;
    case VISC_POWERLAW:
      state.eta = mu_ref*pow(state.T/T_ref, mu_power_law);
      state.lambda = GAMMA*R_gas/(GAMMA-1.0)*state.eta/Pr;
      break;
    case VISC_REFPROP:
      state.eta = fp.eta;
      state.lambda = fp.lambda;
      break;
    }
  }

  /*
   * Calculates the complete thermodynamic state after the conservative variables are advanced in time
   */
  virtual void calcStateVariables()
  {
    for (int icv=0; icv<ncv; icv++)
    {
      if (rho[icv] != rho[icv])
      {
        cout << "density nan at xcv: " << icv << ", " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
        throw(-1);
      }

      if (rho[icv] <= 0.0)
        cout << "negative density at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;

      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      double TKE = 0.0;
      if (kine != NULL)
        TKE = kine[icv];

      double vin    = 1.0/rho[icv];
      double usqhlf = 0.5*vin*vin*vecDotVec3d(rhou[icv], rhou[icv]);
      double ein    = vin*rhoE[icv] - usqhlf - TKE;

      fluidstatejoe_t fp;
      fluidprop_allprops_( "vu", vin, ein*0.001, &fp);

      double pr = fp.P*1.0e5;
      if (pr != pr)
      {
        cout << "pressure nan at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
        cout << "density: " << rho[icv] << " int en: " << ein << " int energy: " << fp.u << endl << " t: " << temp[icv] << " tt: " << fp.T << " press: " << press[icv];
        throw(-1);
      }

      if (pr <= 0.0)
      {
        cout << "negative pressure at xcv: " << pr << ", " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
        throw(-1);
      }
      else              press[icv] = pr;

      temp[icv] = fp.T + 273.15;    // temperature
      enthalpy[icv] = fp.h*1000.0;  // sensible enthalpy

      double dTdrho = -fp.zeta*vin*vin;
      double dTde   =  1.0/(fp.cv*1000.0);
      double dedrho = vin*(usqhlf - ein);
      double dedret = vin;

      dTdr[icv]  = dTdrho + dTde*dedrho;
      dTdre[icv] = dTde*dedret;

      dpdr[icv]   = fp.alpha - fp.u*1.0e3/fp.d*fp.beta + fp.beta/fp.d*usqhlf;
      dpdre[icv] = fp.beta/fp.d;
      sos[icv]   = fp.c;

      if (viscMode == VISC_REFPROP)
      {
        detadr[icv] = fp.deta_drho;
        dlamdr[icv] = fp.dlam_drho;
        detadT[icv] = fp.deta_dT;
        dlamdT[icv] = fp.dlam_dT;
      }
    }

    updateCvData(vel, REPLACE_ROTATE_DATA);
    updateCvData(press, REPLACE_DATA);
    updateCvData(temp, REPLACE_DATA);
    updateCvData(enthalpy, REPLACE_DATA);
    updateCvData(dTdr, REPLACE_DATA);
    updateCvData(dTdre, REPLACE_DATA);
    updateCvData(dpdr, REPLACE_DATA);
    updateCvData(dpdre, REPLACE_DATA);
    updateCvData(sos, REPLACE_DATA);
    updateCvData(detadr, REPLACE_DATA);
    updateCvData(detadT, REPLACE_DATA);
    updateCvData(dlamdr, REPLACE_DATA);
    updateCvData(dlamdT, REPLACE_DATA);

    if (viscMode == NO_VISC) return;

    for (int ifa = nfa_b; ifa < nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];

      double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));

      double temperature = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
      double density     = (w1*rho[icv0] + w0*rho[icv1])/(w0+w1);

      calcLamMuCp(temperature, density, lambda_fa[ifa], mul_fa[ifa], cp_fa[ifa]);
    }
  }
};


/*
 * REAL GAS - LOOK-UP TABLES
 *
 * This class overloads GasModel functions substituting direct calls to FluidProp with
 * look-up table interpolation. The table is created with FluidProp.
 *
 * When using real gas tables, the following lines must be added in Joe.in
 * TMIN (double) minimum value of temperature
 * TMAX (double) maximum value of temperature
 * RHOMIN (double) minimum value of density
 * RHOMAX (double) maximum value of density
 * SPLITAB (bool) 1 to split the table across the saturation line (optional, default 0)
 * CLUSTERING (bool) 1 to clustered nodes around the critical point (optional, default 0)
 * TABLE_DISTR (char) distribution of nodes along density. Options are: UNIFORM, LOG
 * IMAX (int) number of nodes along the density range
 * JMAX (int) number of nodes along the temperature range
 * IMAX_2PHASE (int) number of nodes along the density range in the two phase region table
 * JMAX_2PHASE (int) number of nodes along the temperature range  in the two phase region table
 *
 */
class GasModelTable: public GasModel
{
public:
  double *tab_i, *tab_j;

private:
  fluidstatejoe_t **state;

  char *intScheme;
  string tab_dist;
  bool tab_split;

  int imax, jmax;
  double Tmax, Tmin, rhomax, rhomin;

public:
  // ========================================================================== //
  //                                Constructor                                 //
  // ========================================================================== //
  GasModelTable()
  {
    if (mpi_rank==0)
      cout << "GasModelTable()" << endl;

    // register tab_i and tab_j in order to speed up cell identification
    tab_i = NULL; registerScalar(tab_i, "tab_i", CV_DATA);
    tab_j = NULL; registerScalar(tab_j, "tab_j", CV_DATA);

    Tmin = getDoubleParam("TABLE_TMIN");
    Tmax = getDoubleParam("TABLE_TMAX");
    rhomin = getDoubleParam("TABLE_RHOMIN");
    rhomax = getDoubleParam("TABLE_RHOMAX");

    imax = getIntParam("TABLE_IMAX");
    jmax = getIntParam("TABLE_JMAX");

    // ============================== //
    //        Generate table(s)       //
    // ============================== //
    int method = getIntParam("INTERPOLATION_SCHEME", 0); // 0=bilinear, 1=polynomial
    tab_split = getBoolParam("TABLE_SPLIT", 0);
    tab_dist = getStringParam("TABLE_DISTR", "UNIFORM"); // UNIFORM, LOG


    if (method==0)        intScheme="BILINEAR";
    else if (mpi_rank==0)
    {
      cout << "wrong interpolation scheme!!" << endl;
      throw(-111);
    }

    // initialize the table
    state = new fluidstatejoe_t*[imax];
    for (int i=0; i<imax; i++)
      state[i] = new fluidstatejoe_t[jmax];


    if (!readTableBIN("table.bin"))
    {
      if (mpi_rank==0) // just one core makes the table
      {
        for (int j=0; j<jmax; j++)
        {
          double dens = rhomin + (double)j/(jmax-1.0)*(rhomax-rhomin);
          if (strcmp(tab_dist.c_str(), "LOG")==0)
            dens = pow(10.0, log10(rhomin) + (double)j/(jmax-1.0)*(log10(rhomax)-log10(rhomin)));

          double Tmin_1ph = Tmin;
          if (tab_split)
            if (fabs(fluidprop_vaporqual("Td", Tmin_1ph-273.15, dens))<1.0)
            {
              // bisection to find the
              int nmax=100, it=0;
              double toll=1.0e-6, err=1.0;

              double tA = Tmax;
              double tB = Tmin;

              while (err>toll && it<nmax)
              {
                double tC = 0.5*(tA+tB);
                if (fabs(fluidprop_vaporqual("Td", tC-273.15, dens))<1.0) tB = tC;
                else                                                      tA = tC;

                err = fabs(tA-tB)/tA;
                it++;
              }
              if (it>=nmax-1) printf("bisection not converged at %lf, error: %1.5le\n", dens, err);
              Tmin_1ph = tA;
            }

          for (int i=0; i<imax; i++)
          {
            double temp = Tmin_1ph-273.15 + (double)i/(imax-1.0)*(Tmax-Tmin_1ph);
            fluidprop_allpropsjoe("Td", temp, dens, &state[i][j]);
            if (state[i][j].cp<0.0) state[i][j].cp=0.0;
            if (state[i][j].lambda<0.0) state[i][j].lambda=0.0;
          }
        }
        saveTableBIN("table.bin");
        saveTableTEC("TABLE.plt");
        cout << "Table generated." << endl;
      }
      // every core read the table
      readTableBIN("table.bin");
    }
    else
      if (mpi_rank==0)  cout << "Table read." << endl;
  }

  virtual ~GasModelTable(){}

  // Stores and uses tab_i and tab_j to speed up cell identification
  virtual void calcStateVariables()
  {
    for (int icv=0; icv<ncv; icv++)
    {
      if (rho[icv]<rhomin)
        cout << "rho < rhomin: " << rho[icv] << " at x_cv: " << x_cv[icv][0] << " " << x_cv[icv][1] << x_cv[icv][2] << endl;
      if (rho[icv]>rhomax)
        cout << "rho > rhomax: " << rho[icv] << " at x_cv: " << x_cv[icv][0] << " " << x_cv[icv][1] << x_cv[icv][2] << endl;


      if (rho[icv] != rho[icv])
      {
        cout << "density nan at xcv: " << icv << ", " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
        throw(-1);
      }

      if (rho[icv] <= 0.0)
        cout << "negative density at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;

      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      double TKE = 0.0;
      if (kine != NULL)
        TKE = kine[icv];

      double vin    = 1.0/rho[icv];
      double usqhlf = 0.5*vin*vin*vecDotVec3d(rhou[icv], rhou[icv]);
      double ein    = vin*rhoE[icv] - usqhlf - TKE;

      fluidstatejoe_t fp;
      // the index i,j that identify the cell that contains the thermodynamic state on the table are stored
      // so that they can be used at each time step as the starting point to identify the cell.
      // This speeds up the procedure.
      fp.i = (int)tab_i[icv];
      fp.j = (int)tab_j[icv];
      fluidprop_allprops_( "vu", vin, ein*0.001, &fp, "ALL");
      // store new i,j
      tab_i[icv] = (double)fp.i;
      tab_j[icv] = (double)fp.j;

      double pr = fp.P*1.0e5;
      if (pr != pr)
      {
        cout << "pressure nan at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
        cout << "density: " << rho[icv] << " int en: " << ein << " int energy: " << fp.u << endl << " t: " << temp[icv] << " tt: " << fp.T << " press: " << press[icv];
        throw(-1);
      }

      if (pr <= 0.0)
      {
        printf("%d\t\t%1.5le\t%1.5le\n", icv, vin, ein);
        cout << "negative pressure at xcv: " << pr << ", " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
        throw(-1);
      }
      else press[icv] = pr;

      temp[icv] = fp.T + 273.15;    // temperature
      if (temp[icv]<Tmin)
        cout << "temp < tmin: " << temp[icv] << " at x_cv: " << x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2] << endl;
      if (temp[icv]>Tmax)
        cout << "temp > tmax: " << temp[icv] << " at x_cv: " << x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2] << endl;


      enthalpy[icv] = fp.h*1000.0;  // sensible enthalpy

      double dTdrho = -fp.zeta*vin*vin;
      double dTde   =  1.0/(fp.cv*1000.0);
      // e = e (rho, Et) = 1/rho*(Et - 0.5*|rhou|^2/rho)
      double dedrho = vin*(usqhlf - ein);
      double dedret = vin;

      dTdr[icv]  = dTdrho + dTde*dedrho;
      dTdre[icv] = dTde*dedret;

      dpdr[icv]  = fp.alpha - fp.u*1.0e3/fp.d*fp.beta + fp.beta/fp.d*usqhlf;
      dpdre[icv] = fp.beta/fp.d;
      sos[icv]   = fp.c;

      if (viscMode == VISC_REFPROP)
      {
        detadr[icv] = fp.deta_drho;
        dlamdr[icv] = fp.dlam_drho;
        detadT[icv] = fp.deta_dT;
        dlamdT[icv] = fp.dlam_dT;
      }
    }

    updateCvData(vel, REPLACE_ROTATE_DATA);
    updateCvData(press, REPLACE_DATA);
    updateCvData(temp, REPLACE_DATA);
    updateCvData(enthalpy, REPLACE_DATA);
    updateCvData(dTdr, REPLACE_DATA);
    updateCvData(dTdre, REPLACE_DATA);
    updateCvData(dpdr, REPLACE_DATA);
    updateCvData(dpdre, REPLACE_DATA);
    updateCvData(sos, REPLACE_DATA);
    updateCvData(tab_i, REPLACE_DATA);
    updateCvData(tab_j, REPLACE_DATA);
    updateCvData(detadr, REPLACE_DATA);
    updateCvData(detadT, REPLACE_DATA);
    updateCvData(dlamdr, REPLACE_DATA);
    updateCvData(dlamdT, REPLACE_DATA);

    if (viscMode == NO_VISC) return;

    for (int ifa = nfa_b; ifa < nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];

      double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));

      double temperature = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
      double density     = (w1*rho[icv0] + w0*rho[icv1])/(w0+w1);

      calcLamMuCp(temperature, density, lambda_fa[ifa], mul_fa[ifa], cp_fa[ifa]);
    }
  }


  virtual void allPropsHS(thermoState_t &state, const double h, const double s)
  {
	// since none of the variable used to generate the table is an input, direct access
	// to the tabulated data and interpolation is not possible. We find iteratively
	// density and temperature that correspond to the input enthalpy and entropy.
	// For this reason an initial guess on rho and T is needed.
    fluidstatejoe_t fp;
    fp.d = state.d;
    fp.T = state.T-273.15;
    fluidprop_allprops_("hs", h*1.0e-3, s*1.0e-3, &fp, "ALL");

//    fluidstatejoe_t fp;
//    fluidprop_allpropsjoe("hs", h*1.0e-3, s*1.0e-3, &fp);

    state.cp = fp.cp*1000.0;
    state.h = h;
    state.s = s;
    state.T = fp.T+273.15;
    state.P = fp.P*1.0e5;
    state.d = fp.d;

    switch (viscMode)
    {
    case VISC_SUTHERLAND:
      state.eta = mu_ref*pow(state.T/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(state.T + SL_Sref);
      state.lambda = GAMMA*R_gas/(GAMMA-1.0)*state.eta/Pr;
      break;
    case VISC_POWERLAW:
      state.eta = mu_ref*pow(state.T/T_ref, mu_power_law);
      state.lambda = GAMMA*R_gas/(GAMMA-1.0)*state.eta/Pr;
      break;
    case VISC_REFPROP:
      state.eta = fp.eta;
      state.lambda = fp.lambda;
      break;
    }
  }


  virtual double entropyPv(const double P, const double v)
  {
    fluidstatejoe_t tmp_state;
    tmp_state.i = tmp_state.j = 1;
    fluidprop_allprops_("Pv", P/1.0e5, v, &tmp_state, "ENTROPY");

    return(tmp_state.s*1.0e3);
  }


  virtual double internalEnergyPv(const double press, const double specVol)
  {
    fluidstatejoe_t tmp_state;
    tmp_state.i = tmp_state.j = 1;
    fluidprop_allprops_("Pv", press/1.0e5, specVol, &tmp_state, "ENERGY");

    return(tmp_state.u*1.0e3);
  }

  virtual void calcTherCond(const double P, const double T, double &lam)
  {
    fluidstatejoe_t tmp_state;
    tmp_state.i = tmp_state.j = 1;
    fluidprop_allprops_("PT", P/1.0e5, T-273.15, &tmp_state, "THERMCOND");

    lam = tmp_state.lambda;
  }


  virtual void calcLamMuCp(const double temp, const double rho, double &lambda, double &mul, double &cp)
  {
    switch (viscMode)
    {
    case VISC_SUTHERLAND:
      mul = mu_ref*pow(temp/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp + SL_Sref);
      lambda = GAMMA*R_gas/(GAMMA-1.0)*mul/Pr;
      cp = GAMMA / (GAMMA - 1.0) * R_gas;
      break;
    case VISC_POWERLAW:
      mul = mu_ref*pow(temp/T_ref, mu_power_law);
      lambda = GAMMA*R_gas/(GAMMA-1.0)*mul/Pr;
      cp = GAMMA / (GAMMA - 1.0) * R_gas;
      break;
    case VISC_REFPROP:
      calcViscProps(temp, rho, lambda, mul, cp);
      break;
    }
  }

  void calcViscProps(const double T, const double rho, double &lam, double &eta, double &cp)
  {
    fluidstatejoe_t tmp_state;
    tmp_state.i = tmp_state.j = 1;
    fluidprop_allprops_("Tv", T-273.15, 1.0/rho, &tmp_state, "VISCPROPS");

    lam = tmp_state.lambda;
    eta = tmp_state.eta;
    cp  = tmp_state.cp*1.0e3;
  }


  virtual void calcEnthAndPressDeriv(double &hs, double &dpdr, double &dpdret,
      const double press, const double rho, const double velMagSq)
  {
    fluidstatejoe_t tmp_state;
    tmp_state.i = tmp_state.j = 1;
    fluidprop_allprops_("Pv", press/1.0e5, 1.0/rho, &tmp_state, "ENTHANDPRESSDERIV");

    hs     = tmp_state.h*1.0e3;
    double es = hs - press/rho;
    dpdr   = tmp_state.alpha - es/rho*tmp_state.beta + 0.5*tmp_state.beta/rho*velMagSq;         // dPi/drho
    dpdret = tmp_state.beta/rho;
  }

  virtual double calcPressAndEnthVu(const double v, const double u, double &P, double &h)
  {
    fluidstatejoe_t tmp_state;
    tmp_state.i = tmp_state.j = 1;
    fluidprop_allprops_("vu", v, u/1.0e3, &tmp_state, "PRESSUREANDENTHALPY");

    P = tmp_state.P*1.0e5;
    h = tmp_state.h*1.0e3;
  }

public:
  /*
   * For simplicity, only bilinear interpolation is implemented. It is shown that the accuracy is good
   * for most practical applications using a reasonable number of grid points,
   * see Rinaldi et al., J. Comput. Phys., 2014.
   */
  virtual void fluidprop_allprops_(const char* InputSpec, double Input1, double Input2, struct fluidstatejoe_t *OutputState, const char *props="ALL")
  {
    int indexSnd;
    double sndValue;

    int iTab, jTab;
    double v1, v2, fact1, fact2;

    if ((InputSpec[0]=='v') || (InputSpec[1]=='v'))
    {
      double dens, sndVal;
      if (InputSpec[0]=='v') {dens = 1.0/Input1; sndVal = Input2;}
      else                   {dens = 1.0/Input2; sndVal = Input1;}

      double ff;
      if (strcmp(tab_dist.c_str(), "LOG")==0) ff = (log10(dens)-log10(state[0][0].d))/(log10(state[0][jmax-1].d)-log10(state[0][0].d));
      else                                    ff = (dens-state[0][0].d)/(state[0][jmax-1].d-state[0][0].d);
      jTab = ceil(ff*(double)(jmax-1));

      fact1 = (dens-state[0][jTab-1].d)/(state[0][jTab].d-state[0][jTab-1].d);

      iTab = ceil((double)imax/2.0);
      if ((InputSpec[0]=='u') || (InputSpec[1]=='u'))
      {
        if (sndVal>(state[iTab][jTab-1].u+fact1*(state[iTab][jTab].u-state[iTab][jTab-1].u)))
          while((sndVal>(state[iTab][jTab-1].u+fact1*(state[iTab][jTab].u-state[iTab][jTab-1].u))) && (iTab<imax-1))   iTab++;
        else {while((sndVal<(state[iTab][jTab-1].u+fact1*(state[iTab][jTab].u-state[iTab][jTab-1].u))) && (iTab>0))   iTab--; iTab++;}

        v1 = state[iTab-1][jTab-1].u+fact1*(state[iTab-1][jTab].u-state[iTab-1][jTab-1].u);
        v2 = state[iTab  ][jTab-1].u+fact1*(state[iTab  ][jTab].u-state[iTab  ][jTab-1].u);
      }
      else if ((InputSpec[0]=='P') || (InputSpec[1]=='P'))
      {
        if (sndVal>(state[iTab][jTab-1].P+fact1*(state[iTab][jTab].P-state[iTab][jTab-1].P)))
          while((sndVal>(state[iTab][jTab-1].P+fact1*(state[iTab][jTab].P-state[iTab][jTab-1].P))) && (iTab<imax-1))   iTab++;
        else {while((sndVal<(state[iTab][jTab-1].P+fact1*(state[iTab][jTab].P-state[iTab][jTab-1].P))) && (iTab>0))   iTab--; iTab++;}

        v1 = state[iTab-1][jTab-1].P+fact1*(state[iTab-1][jTab].P-state[iTab-1][jTab-1].P);
        v2 = state[iTab  ][jTab-1].P+fact1*(state[iTab  ][jTab].P-state[iTab  ][jTab-1].P);
      }
      else if ((InputSpec[0]=='T') || (InputSpec[1]=='T'))
      {
        if (sndVal>(state[iTab][jTab-1].T+fact1*(state[iTab][jTab].T-state[iTab][jTab-1].T)))
          while((sndVal>(state[iTab][jTab-1].T+fact1*(state[iTab][jTab].T-state[iTab][jTab-1].T))) && (iTab<imax-1))   iTab++;
        else {while((sndVal<(state[iTab][jTab-1].T+fact1*(state[iTab][jTab].T-state[iTab][jTab-1].T))) && (iTab>0))   iTab--; iTab++;}

        v1 = state[iTab-1][jTab-1].T+fact1*(state[iTab-1][jTab].T-state[iTab-1][jTab-1].T);
        v2 = state[iTab  ][jTab-1].T+fact1*(state[iTab  ][jTab].T-state[iTab  ][jTab-1].T);
      }

      fact2 = (sndVal-v1)/(v2-v1);
    }
    else if (InputSpec=="PT")
    {
      cout << "not yet PT" << endl;
      throw(-111);
    }
    else if (InputSpec=="hs")
    {
      // Inverse evaluation for h and s starting from an initial guess of rho and T.
      // We use a Newton solver and the gradients are approximated at second order.
      double dens = (*OutputState).d;
      double temp = (*OutputState).T;
      double h_tg = Input1;
      double s_tg = Input2;

      double err=1.0, toll=1.0e-10;
      int    it=0, nmax=100;

      iTab = -1;

      while (it<nmax && err>toll)
      {
        it++;

        // find the cell
        double ff;
        if (strcmp(tab_dist.c_str(), "LOG")==0) ff = (log10(dens)-log10(state[0][0].d))/(log10(state[0][jmax-1].d)-log10(state[0][0].d));
        else                                    ff = (dens-state[0][0].d)/(state[0][jmax-1].d-state[0][0].d);
        jTab = ceil(ff*(double)(jmax-1));
        fact1 = (dens-state[0][jTab-1].d)/(state[0][jTab].d-state[0][jTab-1].d);

        iTab=0;
        while((temp>(state[iTab][jTab-1].T+fact1*(state[iTab][jTab].T-state[iTab][jTab-1].T))) && (iTab<imax-1))   iTab++;

        // calculate fact2
        v1 = state[iTab-1][jTab-1].T+fact1*(state[iTab-1][jTab].T-state[iTab-1][jTab-1].T);
        v2 = state[iTab  ][jTab-1].T+fact1*(state[iTab  ][jTab].T-state[iTab  ][jTab-1].T);
        fact2 = (temp-v1)/(v2-v1);

        double drho = state[0][jTab].d - state[0][jTab-1].d;
        double dT   = v2-v1;

        // calculate h* and s*
        bilinInterpolState(*OutputState, iTab, jTab, fact1, fact2, "HS");
        iTab = (*OutputState).i;

        // dhdT, dsdT
        double h1 = state[iTab-1][jTab-1].h+fact1*(state[iTab-1][jTab].h-state[iTab-1][jTab-1].h);
        double h2 = state[iTab  ][jTab-1].h+fact1*(state[iTab  ][jTab].h-state[iTab  ][jTab-1].h);
        double dhdT = (h2-h1)/dT;
        double s1 = state[iTab-1][jTab-1].s+fact1*(state[iTab-1][jTab].s-state[iTab-1][jTab-1].s);
        double s2 = state[iTab  ][jTab-1].s+fact1*(state[iTab  ][jTab].s-state[iTab  ][jTab-1].s);
        double dsdT = (s2-s1)/dT;


        // dhdrho, dsdrho
        double f1 = (temp-state[iTab-1][jTab-1].T)/(state[iTab][jTab-1].T-state[iTab-1][jTab-1].T);
        double f2 = (temp-state[iTab-1][jTab].T)/(state[iTab][jTab].T-state[iTab-1][jTab].T);
        h1 = state[iTab-1][jTab-1].h+f1*(state[iTab][jTab-1].h-state[iTab-1][jTab-1].h);
        h2 = state[iTab-1][jTab  ].h+f2*(state[iTab][jTab  ].h-state[iTab-1][jTab  ].h);
        double dhdrho = (h2-h1)/drho;
        s1 = state[iTab-1][jTab-1].s+f1*(state[iTab][jTab-1].s-state[iTab-1][jTab-1].s);
        s2 = state[iTab-1][jTab  ].s+f2*(state[iTab][jTab  ].s-state[iTab-1][jTab  ].s);
        double dsdrho = (s2-s1)/drho;

        double det = dsdrho*dhdT - dsdT*dhdrho;

        double ddens = (dhdT*((*OutputState).s - s_tg) - dsdT*((*OutputState).h - h_tg))/det;
        double dtemp = (-dhdrho*((*OutputState).s - s_tg) + dsdrho*((*OutputState).h - h_tg))/det;

        err = fabs(ddens)/dens + fabs(dtemp)/temp;

        if (fabs(ddens)<dens) dens -= 0.75*ddens;
        if (fabs(dtemp)<temp) temp -= 0.75*dtemp;
      }
      if (err>toll || it>nmax) cout << "HS not converged!" << endl;
    }

    bilinInterpolState(*OutputState, iTab, jTab, fact1, fact2, props);

  }


  void bilinInterpolState(fluidstatejoe_t &resState, const int i, const int j, const double fact1, const double fact2, const char *props)
    {
      if (props=="ALL")
      {
        resState.v      = bilinInterpol(state[i-1][j-1].v, state[i-1][j].v, state[i][j-1].v, state[i][j].v, fact1, fact2);
        resState.d      = bilinInterpol(state[i-1][j-1].d, state[i-1][j].d, state[i][j-1].d, state[i][j].d, fact1, fact2);
        resState.P      = bilinInterpol(state[i-1][j-1].P, state[i-1][j].P, state[i][j-1].P, state[i][j].P, fact1, fact2);
        resState.T      = bilinInterpol(state[i-1][j-1].T, state[i-1][j].T, state[i][j-1].T, state[i][j].T, fact1, fact2);
        resState.u      = bilinInterpol(state[i-1][j-1].u, state[i-1][j].u, state[i][j-1].u, state[i][j].u, fact1, fact2);
        resState.h      = bilinInterpol(state[i-1][j-1].h, state[i-1][j].h, state[i][j-1].h, state[i][j].h, fact1, fact2);
        resState.s      = bilinInterpol(state[i-1][j-1].s, state[i-1][j].s, state[i][j-1].s, state[i][j].s, fact1, fact2);
        resState.cp     = bilinInterpol(state[i-1][j-1].cp, state[i-1][j].cp, state[i][j-1].cp, state[i][j].cp, fact1, fact2);
        resState.cv     = bilinInterpol(state[i-1][j-1].cv, state[i-1][j].cv, state[i][j-1].cv, state[i][j].cv, fact1, fact2);
        resState.c      = bilinInterpol(state[i-1][j-1].c, state[i-1][j].c, state[i][j-1].c, state[i][j].c, fact1, fact2);
        resState.alpha  = bilinInterpol(state[i-1][j-1].alpha, state[i-1][j].alpha, state[i][j-1].alpha, state[i][j].alpha, fact1, fact2);
        resState.beta   = bilinInterpol(state[i-1][j-1].beta, state[i-1][j].beta, state[i][j-1].beta, state[i][j].beta, fact1, fact2);
        resState.eta    = bilinInterpol(state[i-1][j-1].eta, state[i-1][j].eta, state[i][j-1].eta, state[i][j].eta, fact1, fact2);
        resState.lambda = bilinInterpol(state[i-1][j-1].lambda, state[i-1][j].lambda, state[i][j-1].lambda, state[i][j].lambda, fact1, fact2);
        resState.zeta   = bilinInterpol(state[i-1][j-1].zeta, state[i-1][j].zeta, state[i][j-1].zeta, state[i][j].zeta, fact1, fact2);
        resState.deta_drho = bilinInterpol(state[i-1][j-1].deta_drho, state[i-1][j].deta_drho, state[i][j-1].deta_drho, state[i][j].deta_drho, fact1, fact2);
        resState.deta_dT   = bilinInterpol(state[i-1][j-1].deta_dT, state[i-1][j].deta_dT, state[i][j-1].deta_dT, state[i][j].deta_dT, fact1, fact2);
        resState.dlam_drho = bilinInterpol(state[i-1][j-1].dlam_drho, state[i-1][j].dlam_drho, state[i][j-1].dlam_drho, state[i][j].dlam_drho, fact1, fact2);
        resState.dlam_dT   = bilinInterpol(state[i-1][j-1].dlam_dT, state[i-1][j].dlam_dT, state[i][j-1].dlam_dT, state[i][j].dlam_dT, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="ENERGY")
      {
        resState.u      = bilinInterpol(state[i-1][j-1].u, state[i-1][j].u, state[i][j-1].u, state[i][j].u, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="ENTROPY")
      {
        resState.s      = bilinInterpol(state[i-1][j-1].s, state[i-1][j].s, state[i][j-1].s, state[i][j].s, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="THERMCOND")
      {
        resState.lambda = bilinInterpol(state[i-1][j-1].lambda, state[i-1][j].lambda, state[i][j-1].lambda, state[i][j].lambda, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="HS")
      {
        resState.h = bilinInterpol(state[i-1][j-1].h, state[i-1][j].h, state[i][j-1].h, state[i][j].h, fact1, fact2);
        resState.s = bilinInterpol(state[i-1][j-1].s, state[i-1][j].s, state[i][j-1].s, state[i][j].s, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="VISCPROPS")
      {
        resState.lambda = bilinInterpol(state[i-1][j-1].lambda, state[i-1][j].lambda, state[i][j-1].lambda, state[i][j].lambda, fact1, fact2);
        resState.eta    = bilinInterpol(state[i-1][j-1].eta, state[i-1][j].eta, state[i][j-1].eta, state[i][j].eta, fact1, fact2);
        resState.cp     = bilinInterpol(state[i-1][j-1].cp, state[i-1][j].cp, state[i][j-1].cp, state[i][j].cp, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="PRESSUREANDENTHALPY")
      {
        resState.P      = bilinInterpol(state[i-1][j-1].P, state[i-1][j].P, state[i][j-1].P, state[i][j].P, fact1, fact2);
        resState.h      = bilinInterpol(state[i-1][j-1].h, state[i-1][j].h, state[i][j-1].h, state[i][j].h, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="ENTHANDPRESSDERIV")
      {
        resState.h      = bilinInterpol(state[i-1][j-1].h, state[i-1][j].h, state[i][j-1].h, state[i][j].h, fact1, fact2);
        resState.alpha  = bilinInterpol(state[i-1][j-1].alpha, state[i-1][j].alpha, state[i][j-1].alpha, state[i][j].alpha, fact1, fact2);
        resState.beta   = bilinInterpol(state[i-1][j-1].beta, state[i-1][j].beta, state[i][j-1].beta, state[i][j].beta, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (mpi_rank==0) {cout << "wrong input props in bilinInterpolState: " << props << endl; throw(-1);}
    }

    inline double bilinInterpol(double a1, double a2, double b1, double b2, double fact1, double fact2)
    {
      double v1 = a1+fact1*(a2-a1);
      double v2 = b1+fact1*(b2-b1);
      return (v1 + fact2*(v2-v1));
    }



  // ========================================================================== //
  //                                Read and write                              //
  // ========================================================================== //
  // Read table from binary file
  int readTableBIN(const char *fileName)
  {
    FILE *fp = fopen(fileName, "rb");
    if (fp == NULL) return 0;

    int imaxR, jmaxR;
    fread(&imaxR, sizeof(int), 1, fp);
    fread(&jmaxR, sizeof(int), 1, fp);

    if ((imax != imaxR) || (jmax != jmaxR))
    {
      cout << "Tables dimension error!" << endl;
      throw(-1);
    }

    for (int i=0; i<imax; i++)
      for (int j=0; j<jmax; j++)
        fread(&state[i][j], sizeof(fluidstatejoe_t), 1, fp);

    fclose(fp);

    return 1;
  }

  void saveTableBIN(const char *fileName)
  {
    FILE *fp = fopen(fileName, "wb");

    fwrite(&imax, sizeof(int), 1, fp);
    fwrite(&jmax, sizeof(int), 1, fp);

    for (int i=0; i<imax; i++)
      for (int j=0; j<jmax; j++)
        fwrite(&state[i][j], sizeof(fluidstatejoe_t), 1, fp);

    fclose(fp);
  }


  // Save table in Tecplot format
  void saveTableTEC(const char *fileName)
  {
    FILE *fp = fopen(fileName, "wt");
    fprintf(fp, "TITLE=\"table\"\n");
    fprintf(fp, "VARIABLES=\"d\" \"T\" \"v\" \"s\" \"P\" \"u\" \"h\" \"eta\" \"lambda\" \"cp\" \"cv\" \"c\" \"alpha\" \"beta\" \"zeta\"\n");// \"deta_drho\" \"deta_dT\" \"dlam_drho\" \"dlam_dT\" \"one_cf\"\n");
    fprintf(fp, "ZONE I=%d, J=%d DATAPACKING=POINT\n", imax, jmax);

    double mmol = fluidprop_mmol();
    double R = 8.314462175 / mmol;

    for (int j=0; j<jmax; j++)
      for (int i=0; i<imax; i++)
        fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",//%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",
                state[i][j].d, state[i][j].T, state[i][j].v, state[i][j].s, state[i][j].P, state[i][j].u, state[i][j].h, state[i][j].eta, state[i][j].lambda,
                state[i][j].cp, state[i][j].cv, state[i][j].c, state[i][j].alpha, state[i][j].beta, state[i][j].zeta);//, state[i][j].deta_drho, state[i][j].deta_dT,
                //state[i][j].dlam_drho, state[i][j].dlam_dT, 1.0 - (state[i][j].P*1.0e5/state[i][j].d/R/(state[i][j].T+273.15)));
    fclose(fp);
  }


};

#else /* NO FLUIDPROP */
class GasModel {
public:
  GasModel(){}
  virtual ~GasModel(){}
};
class GasModelTable: public GasModel {
public:
  GasModelTable(){}
  virtual ~GasModelTable(){}
};
#endif

#endif /* GASMODEL_H_ */
