/*!
 * look_up_table.cpp
 * \brief Source of the main thermo-physical subroutines of the SU2 solvers.
 * \author S.Smit, R.Pecnik
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/look_up_table.hpp"

CLookUpTable::CLookUpTable(void) {
	imax = 0;
	jmax = 0;
	TableState = NULL;
	OutputState = NULL;
	table_filename = "";
	table_fluid = "";

}

CLookUpTable::CLookUpTable(CConfig *config){
	int rank = MASTER_NODE;
	#ifdef HAVE_MPI
	  int rank;
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	table_filename = config->GetTable_Name();
	table_fluid= config->GetTable_Fluid();
	table_distribution = config->GetTable_Distribution();
	table_interpolation_scheme = config->GetTable_InterpolationScheme();
	rhomin = config->GetTable_RhoMin();
	rhomax = config->GetTable_RhoMax();
	Tmin = config->GetTable_TMin();
	Tmax = config->GetTable_TMax();
	imax = config->GetTable_IMax();
	jmax = config->GetTable_JMax();

	switch (config->GetTable_InterpolationScheme()) {
	  case BILINEAR:
	    table_interpolation_scheme = "BILINEAR";
	    break;
	  case POLYNOMIAL:
	    table_interpolation_scheme = "POLYNOMIAL";
	    break;
	}
	switch (config->GetTable_Distribution()) {
	  case UNIFORM:
		table_distribution = "UNIFORM";
		break;
	  case LOGARITMIC:
		table_distribution = "LOGARITMIC";
		break;
	}

	OutputState = new FluidState;
	TableState = new FluidState*[imax];
    for (int i=0; i<imax; i++)
    	TableState[i] = new FluidState[jmax];

    if (rank == MASTER_NODE){
    	CreateTableRhoT();
    }
    else {
    	ReadTableRhoT(table_filename);
    }


}
CLookUpTable::CLookUpTable(string table_name,
						   string fluid,
						   string tab_dist,
						   int table_imax,
						   int table_jmax,
						   double rho_min,
						   double rho_max,
						   double T_min,
						   double T_max,
						   string interpolation_scheme
						   ){
	int rank = MASTER_NODE;
	#ifdef HAVE_MPI
	  int rank;
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	table_filename = table_name;
	table_fluid= fluid;
	table_distribution  = tab_dist;
	table_interpolation_scheme = interpolation_scheme;
	imax = table_imax;
	jmax = table_jmax;
	rhomin = rho_min;
	rhomax = rho_max;
	Tmin = T_min;
	Tmax = T_max;
	OutputState = new FluidState;
	TableState = new FluidState*[imax];
    for (int i=0; i<imax; i++)
    	TableState[i] = new FluidState[jmax];

    if (rank == MASTER_NODE){
    	CreateTableRhoT();
    }
    else {
    	ReadTableRhoT(table_filename);
    }

}



CLookUpTable::~CLookUpTable(void) {

}

void CLookUpTable::SetTDState_rhoe (su2double rho, su2double e) {}


void CLookUpTable::SetTDState_PT (su2double P, su2double T ){};


void CLookUpTable::SetTDState_Prho (su2double P, su2double rho ){};


void CLookUpTable::SetEnergy_Prho (su2double P, su2double rho ){};


void CLookUpTable::SetTDState_hs (su2double h, su2double s ){};


void CLookUpTable::SetTDState_rhoT (su2double rho, su2double T ){
	if (not CheckIfInterpolated_rhoT(rho,T)) InterpolateProperties("Tv",(double) T, 1.0/rho, OutputState, "ALL" );
};

void CLookUpTable::SetTDState_Ps (su2double P, su2double s ){};

//void CLookUpTable::ComputeDerivativeNRBC_Prho (su2double P, su2double rho );

// table functions

void CLookUpTable::CreateTableRhoT(void){
    for (int j=0; j<jmax; j++)
    {
      double dens = rhomin + (double)j/(jmax-1.0)*(rhomax-rhomin);
      if (strcmp(table_distribution.c_str(), "LOG")==0)
        dens = pow(10.0, log10(rhomin) + (double)j/(jmax-1.0)*(log10(rhomax)-log10(rhomin)));
      double Tmin_1ph = Tmin;
      for (int i=0; i<imax; i++)
      {
		if (fabs(CoolPropGetVaporQuality("TD",table_fluid, Tmin_1ph, dens))<1.0)
		{
		// bisection to find the
		int nmax=100, it=0;
		double toll=1.0e-6, err=1.0;

		double tA = Tmax;
		double tB = Tmin;

		while (err>toll && it<nmax)
		{
		  double tC = 0.5*(tA+tB);
		  if (fabs(CoolPropGetVaporQuality("TD", table_fluid, tC, dens))<1.0) tB = tC;
		  else                                                      tA = tC;

		  err = fabs(tA-tB)/tA;
		  it++;
		}
		if (it>=nmax-1) printf("bisection not converged at %lf, error: %1.5le\n", dens, err);
		Tmin_1ph = tA;
		}

        double temp = Tmin_1ph-273.15 + (double)i/(imax-1.0)*(Tmax-Tmin_1ph);
        CoolPropGetFluidState("TD", table_fluid, temp+273.15, dens, TableState[i][j]);
        if (TableState[i][j].cp<0.0) TableState[i][j].cp=0.0;
        if (TableState[i][j].lambda<0.0) TableState[i][j].lambda=0.0;
      }
    }
	SaveTableTEC(table_filename.c_str());
//	SaveTableBIN(table_filename.c_str());

}

//void CLookUpTable::CreateTransferTable(string InputSpec, string OutputSpec){
//
//
//}

int CLookUpTable::ReadTableBIN(const char *fileName)
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
      fread(&TableState[i][j], sizeof(FluidState), 1, fp);
  fclose(fp);
  return 1;
}

void CLookUpTable::SaveTableTEC(const char *fileName)
{
  FILE *fp = fopen(fileName, "wt");
  fprintf(fp, "TITLE=\"table\"\n");
  fprintf(fp, "VARIABLES=\"d\" \"T\" \"v\" \"s\" \"P\" \"u\" \"h\" \"eta\" \"lambda\" \"cp\" \"cv\" \"c\" \"alpha\" \"beta\" \"zeta\"\n");// \"deta_drho\" \"deta_dT\" \"dlam_drho\" \"dlam_dT\" \"one_cf\"\n");
  fprintf(fp, "ZONE I=%d, J=%d DATAPACKING=POINT\n", imax, jmax);

  for (int j=0; j<jmax; j++)
    for (int i=0; i<imax; i++)
      fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",//%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",
    		  TableState[i][j].d, TableState[i][j].T, TableState[i][j].v, TableState[i][j].s, TableState[i][j].P, TableState[i][j].u, TableState[i][j].h, TableState[i][j].eta, TableState[i][j].lambda,
    		  TableState[i][j].cp, TableState[i][j].cv, TableState[i][j].c, TableState[i][j].alpha, TableState[i][j].beta, TableState[i][j].zeta);//, state[i][j].deta_drho, state[i][j].deta_dT,
              //state[i][j].dlam_drho, state[i][j].dlam_dT, 1.0 - (state[i][j].P*1.0e5/state[i][j].d/R/(state[i][j].T+273.15)));
  fclose(fp);
}

void CLookUpTable::SaveTableBIN(const char *fileName)
{
  FILE *fp = fopen(fileName, "wb");

  fwrite(&imax, sizeof(int), 1, fp);
  fwrite(&jmax, sizeof(int), 1, fp);

  for (int i=0; i<imax; i++)
    for (int j=0; j<jmax; j++)
      fwrite(&TableState[i][j], sizeof(FluidState), 1, fp);

  fclose(fp);
}

void CLookUpTable::ReadTableRhoT(string filename){
    ReadTableBIN(filename.c_str());
}

//interpolation functions

bool CLookUpTable::CheckIfInterpolated_rhoe(su2double rho, su2double e ){
	return (OutputState->d == rho ) && (OutputState->u == e);
}

bool CLookUpTable::CheckIfInterpolated_PT(su2double P, su2double T ){
	return (OutputState->P == P ) && (OutputState->T == T);
}

bool CLookUpTable::CheckIfInterpolated_Prho(su2double P, su2double rho ){
	return (OutputState->P == P ) && (OutputState->d  == rho);
}

bool CLookUpTable::CheckIfInterpolated_hs(su2double h, su2double s ){
	return (OutputState->h  == h ) && (OutputState->s  == s);
}

bool CLookUpTable::CheckIfInterpolated_rhoT(su2double rho, su2double T ){
	return (OutputState->d  == rho ) && (OutputState->T  == T);
}

bool CLookUpTable::CheckIfInterpolated_Ps(su2double P, su2double s ){
	return (OutputState->P == P ) && (OutputState->s == s);
}

void CLookUpTable::InterpolateProperties( const char* InputSpec,
										double Input1,
										double Input2,
										struct FluidState *OutputState,
										const char *props){
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
      if (strcmp(table_distribution.c_str(), "LOG")==0) ff = (log10(dens)-log10(TableState[0][0].d))/(log10(TableState[0][jmax-1].d)-log10(TableState[0][0].d));
      else                                    ff = (dens-TableState[0][0].d)/(TableState[0][jmax-1].d-TableState[0][0].d);
      jTab = ceil(ff*(double)(jmax-1));

      fact1 = (dens-TableState[0][jTab-1].d)/(TableState[0][jTab].d-TableState[0][jTab-1].d);

      iTab = ceil((double)imax/2.0);
      if ((InputSpec[0]=='u') || (InputSpec[1]=='u'))
      {
        if (sndVal>(TableState[iTab][jTab-1].u+fact1*(TableState[iTab][jTab].u-TableState[iTab][jTab-1].u)))
          while((sndVal>(TableState[iTab][jTab-1].u+fact1*(TableState[iTab][jTab].u-TableState[iTab][jTab-1].u))) && (iTab<imax-1))   iTab++;
        else {while((sndVal<(TableState[iTab][jTab-1].u+fact1*(TableState[iTab][jTab].u-TableState[iTab][jTab-1].u))) && (iTab>0))   iTab--; iTab++;}

        v1 = TableState[iTab-1][jTab-1].u+fact1*(TableState[iTab-1][jTab].u-TableState[iTab-1][jTab-1].u);
        v2 = TableState[iTab  ][jTab-1].u+fact1*(TableState[iTab  ][jTab].u-TableState[iTab  ][jTab-1].u);
      }
      else if ((InputSpec[0]=='P') || (InputSpec[1]=='P'))
      {
        if (sndVal>(TableState[iTab][jTab-1].P+fact1*(TableState[iTab][jTab].P-TableState[iTab][jTab-1].P)))
          while((sndVal>(TableState[iTab][jTab-1].P+fact1*(TableState[iTab][jTab].P-TableState[iTab][jTab-1].P))) && (iTab<imax-1))   iTab++;
        else {while((sndVal<(TableState[iTab][jTab-1].P+fact1*(TableState[iTab][jTab].P-TableState[iTab][jTab-1].P))) && (iTab>0))   iTab--; iTab++;}

        v1 = TableState[iTab-1][jTab-1].P+fact1*(TableState[iTab-1][jTab].P-TableState[iTab-1][jTab-1].P);
        v2 = TableState[iTab  ][jTab-1].P+fact1*(TableState[iTab  ][jTab].P-TableState[iTab  ][jTab-1].P);
      }
      else if ((InputSpec[0]=='T') || (InputSpec[1]=='T'))
      {
        if (sndVal>(TableState[iTab][jTab-1].T+fact1*(TableState[iTab][jTab].T-TableState[iTab][jTab-1].T)))
          while((sndVal>(TableState[iTab][jTab-1].T+fact1*(TableState[iTab][jTab].T-TableState[iTab][jTab-1].T))) && (iTab<imax-1))   iTab++;
        else {while((sndVal<(TableState[iTab][jTab-1].T+fact1*(TableState[iTab][jTab].T-TableState[iTab][jTab-1].T))) && (iTab>0))   iTab--; iTab++;}

        v1 = TableState[iTab-1][jTab-1].T+fact1*(TableState[iTab-1][jTab].T-TableState[iTab-1][jTab-1].T);
        v2 = TableState[iTab  ][jTab-1].T+fact1*(TableState[iTab  ][jTab].T-TableState[iTab  ][jTab-1].T);
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
        if (strcmp(table_distribution.c_str(), "LOG")==0) ff = (log10(dens)-log10(TableState[0][0].d))/(log10(TableState[0][jmax-1].d)-log10(TableState[0][0].d));
        else                                    ff = (dens-TableState[0][0].d)/(TableState[0][jmax-1].d-TableState[0][0].d);
        jTab = ceil(ff*(double)(jmax-1));
        fact1 = (dens-TableState[0][jTab-1].d)/(TableState[0][jTab].d-TableState[0][jTab-1].d);

        iTab=0;
        while((temp>(TableState[iTab][jTab-1].T+fact1*(TableState[iTab][jTab].T-TableState[iTab][jTab-1].T))) && (iTab<imax-1))   iTab++;

        // calculate fact2
        v1 = TableState[iTab-1][jTab-1].T+fact1*(TableState[iTab-1][jTab].T-TableState[iTab-1][jTab-1].T);
        v2 = TableState[iTab  ][jTab-1].T+fact1*(TableState[iTab  ][jTab].T-TableState[iTab  ][jTab-1].T);
        fact2 = (temp-v1)/(v2-v1);

        double drho = TableState[0][jTab].d - TableState[0][jTab-1].d;
        double dT   = v2-v1;

        // calculate h* and s*
        bilinInterpolState(*OutputState, iTab, jTab, fact1, fact2, "HS");
        iTab = (*OutputState).i;

        // dhdT, dsdT
        double h1 = TableState[iTab-1][jTab-1].h+fact1*(TableState[iTab-1][jTab].h-TableState[iTab-1][jTab-1].h);
        double h2 = TableState[iTab  ][jTab-1].h+fact1*(TableState[iTab  ][jTab].h-TableState[iTab  ][jTab-1].h);
        double dhdT = (h2-h1)/dT;
        double s1 = TableState[iTab-1][jTab-1].s+fact1*(TableState[iTab-1][jTab].s-TableState[iTab-1][jTab-1].s);
        double s2 = TableState[iTab  ][jTab-1].s+fact1*(TableState[iTab  ][jTab].s-TableState[iTab  ][jTab-1].s);
        double dsdT = (s2-s1)/dT;


        // dhdrho, dsdrho
        double f1 = (temp-TableState[iTab-1][jTab-1].T)/(TableState[iTab][jTab-1].T-TableState[iTab-1][jTab-1].T);
        double f2 = (temp-TableState[iTab-1][jTab].T)/(TableState[iTab][jTab].T-TableState[iTab-1][jTab].T);
        h1 = TableState[iTab-1][jTab-1].h+f1*(TableState[iTab][jTab-1].h-TableState[iTab-1][jTab-1].h);
        h2 = TableState[iTab-1][jTab  ].h+f2*(TableState[iTab][jTab  ].h-TableState[iTab-1][jTab  ].h);
        double dhdrho = (h2-h1)/drho;
        s1 = TableState[iTab-1][jTab-1].s+f1*(TableState[iTab][jTab-1].s-TableState[iTab-1][jTab-1].s);
        s2 = TableState[iTab-1][jTab  ].s+f2*(TableState[iTab][jTab  ].s-TableState[iTab-1][jTab  ].s);
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

};

void CLookUpTable::bilinInterpolState(FluidState &resState,
									  const int i,
									  const int j,
							          const double fact1,
									  const double fact2,
									  const char *props){
	if (props=="ALL") {
        resState.v      = bilinInterpol(TableState[i-1][j-1].v, TableState[i-1][j].v, TableState[i][j-1].v, TableState[i][j].v, fact1, fact2);
        resState.d      = bilinInterpol(TableState[i-1][j-1].d, TableState[i-1][j].d, TableState[i][j-1].d, TableState[i][j].d, fact1, fact2);
        resState.P      = bilinInterpol(TableState[i-1][j-1].P, TableState[i-1][j].P, TableState[i][j-1].P, TableState[i][j].P, fact1, fact2);
        resState.T      = bilinInterpol(TableState[i-1][j-1].T, TableState[i-1][j].T, TableState[i][j-1].T, TableState[i][j].T, fact1, fact2);
        resState.u      = bilinInterpol(TableState[i-1][j-1].u, TableState[i-1][j].u, TableState[i][j-1].u, TableState[i][j].u, fact1, fact2);
        resState.h      = bilinInterpol(TableState[i-1][j-1].h, TableState[i-1][j].h, TableState[i][j-1].h, TableState[i][j].h, fact1, fact2);
        resState.s      = bilinInterpol(TableState[i-1][j-1].s, TableState[i-1][j].s, TableState[i][j-1].s, TableState[i][j].s, fact1, fact2);
        resState.cp     = bilinInterpol(TableState[i-1][j-1].cp, TableState[i-1][j].cp, TableState[i][j-1].cp, TableState[i][j].cp, fact1, fact2);
        resState.cv     = bilinInterpol(TableState[i-1][j-1].cv, TableState[i-1][j].cv, TableState[i][j-1].cv, TableState[i][j].cv, fact1, fact2);
        resState.c      = bilinInterpol(TableState[i-1][j-1].c, TableState[i-1][j].c, TableState[i][j-1].c, TableState[i][j].c, fact1, fact2);
        resState.alpha  = bilinInterpol(TableState[i-1][j-1].alpha, TableState[i-1][j].alpha, TableState[i][j-1].alpha, TableState[i][j].alpha, fact1, fact2);
        resState.beta   = bilinInterpol(TableState[i-1][j-1].beta, TableState[i-1][j].beta, TableState[i][j-1].beta, TableState[i][j].beta, fact1, fact2);
        resState.eta    = bilinInterpol(TableState[i-1][j-1].eta, TableState[i-1][j].eta, TableState[i][j-1].eta, TableState[i][j].eta, fact1, fact2);
        resState.lambda = bilinInterpol(TableState[i-1][j-1].lambda, TableState[i-1][j].lambda, TableState[i][j-1].lambda, TableState[i][j].lambda, fact1, fact2);
        resState.zeta   = bilinInterpol(TableState[i-1][j-1].zeta, TableState[i-1][j].zeta, TableState[i][j-1].zeta, TableState[i][j].zeta, fact1, fact2);
        resState.deta_drho = bilinInterpol(TableState[i-1][j-1].deta_drho, TableState[i-1][j].deta_drho, TableState[i][j-1].deta_drho, TableState[i][j].deta_drho, fact1, fact2);
        resState.deta_dT   = bilinInterpol(TableState[i-1][j-1].deta_dT, TableState[i-1][j].deta_dT, TableState[i][j-1].deta_dT, TableState[i][j].deta_dT, fact1, fact2);
        resState.dlam_drho = bilinInterpol(TableState[i-1][j-1].dlam_drho, TableState[i-1][j].dlam_drho, TableState[i][j-1].dlam_drho, TableState[i][j].dlam_drho, fact1, fact2);
        resState.dlam_dT   = bilinInterpol(TableState[i-1][j-1].dlam_dT, TableState[i-1][j].dlam_dT, TableState[i][j-1].dlam_dT, TableState[i][j].dlam_dT, fact1, fact2);
        resState.dPdrho_e  = bilinInterpol(TableState[i-1][j-1].dPdrho_e, TableState[i-1][j].dPdrho_e, TableState[i][j-1].dPdrho_e, TableState[i][j].dPdrho_e, fact1, fact2);
        resState.dPde_rho  = bilinInterpol(TableState[i-1][j-1].dPde_rho, TableState[i-1][j].dPde_rho, TableState[i][j-1].dPde_rho, TableState[i][j].dPde_rho, fact1, fact2);
        resState.dTdrho_e  = bilinInterpol(TableState[i-1][j-1].dTdrho_e, TableState[i-1][j].dTdrho_e, TableState[i][j-1].dTdrho_e, TableState[i][j].dTdrho_e, fact1, fact2);
        resState.dTde_rho  = bilinInterpol(TableState[i-1][j-1].dTde_rho, TableState[i-1][j].dTde_rho, TableState[i][j-1].dTde_rho, TableState[i][j].dTde_rho, fact1, fact2);
        resState.dhdrho_P  = bilinInterpol(TableState[i-1][j-1].dhdrho_P, TableState[i-1][j].dhdrho_P, TableState[i][j-1].dhdrho_P, TableState[i][j].dhdrho_P, fact1, fact2);
        resState.dhdP_rho  = bilinInterpol(TableState[i-1][j-1].dhdP_rho, TableState[i-1][j].dhdP_rho, TableState[i][j-1].dhdP_rho, TableState[i][j].dhdP_rho, fact1, fact2);
        resState.dsdrho_P  = bilinInterpol(TableState[i-1][j-1].dsdrho_P, TableState[i-1][j].dsdrho_P, TableState[i][j-1].dsdrho_P, TableState[i][j].dsdrho_P, fact1, fact2);
        resState.dsdP_rho  = bilinInterpol(TableState[i-1][j-1].dsdP_rho, TableState[i-1][j].dsdP_rho, TableState[i][j-1].dsdP_rho, TableState[i][j].dsdP_rho, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="ENERGY")
      {
        resState.u      = bilinInterpol(TableState[i-1][j-1].u, TableState[i-1][j].u, TableState[i][j-1].u, TableState[i][j].u, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="ENTROPY")
      {
        resState.s      = bilinInterpol(TableState[i-1][j-1].s, TableState[i-1][j].s, TableState[i][j-1].s, TableState[i][j].s, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="THERMCOND")
      {
        resState.lambda = bilinInterpol(TableState[i-1][j-1].lambda, TableState[i-1][j].lambda, TableState[i][j-1].lambda, TableState[i][j].lambda, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="HS")
      {
        resState.h = bilinInterpol(TableState[i-1][j-1].h, TableState[i-1][j].h, TableState[i][j-1].h, TableState[i][j].h, fact1, fact2);
        resState.s = bilinInterpol(TableState[i-1][j-1].s, TableState[i-1][j].s, TableState[i][j-1].s, TableState[i][j].s, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="VISCPROPS")
      {
        resState.lambda = bilinInterpol(TableState[i-1][j-1].lambda, TableState[i-1][j].lambda, TableState[i][j-1].lambda, TableState[i][j].lambda, fact1, fact2);
        resState.eta    = bilinInterpol(TableState[i-1][j-1].eta, TableState[i-1][j].eta, TableState[i][j-1].eta, TableState[i][j].eta, fact1, fact2);
        resState.cp     = bilinInterpol(TableState[i-1][j-1].cp, TableState[i-1][j].cp, TableState[i][j-1].cp, TableState[i][j].cp, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="PRESSUREANDENTHALPY")
      {
        resState.P      = bilinInterpol(TableState[i-1][j-1].P, TableState[i-1][j].P, TableState[i][j-1].P, TableState[i][j].P, fact1, fact2);
        resState.h      = bilinInterpol(TableState[i-1][j-1].h, TableState[i-1][j].h, TableState[i][j-1].h, TableState[i][j].h, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
      else if (props=="ENTHANDPRESSDERIV")
      {
        resState.h      = bilinInterpol(TableState[i-1][j-1].h, TableState[i-1][j].h, TableState[i][j-1].h, TableState[i][j].h, fact1, fact2);
        resState.alpha  = bilinInterpol(TableState[i-1][j-1].alpha, TableState[i-1][j].alpha, TableState[i][j-1].alpha, TableState[i][j].alpha, fact1, fact2);
        resState.beta   = bilinInterpol(TableState[i-1][j-1].beta, TableState[i-1][j].beta, TableState[i][j-1].beta, TableState[i][j].beta, fact1, fact2);
        resState.i      = (double)i;
        resState.j      = (double)j;
      }
//      else if (mpi_rank==0) {cout << "wrong input props in bilinInterpolState: " << props << endl; throw(-1);}
};

double CLookUpTable::bilinInterpol( double a1,
									       double a2,
										   double b1,
										   double b2,
										   double fact1,
										   double fact2){
  double v1 = a1+fact1*(a2-a1);
  double v2 = b1+fact1*(b2-b1);
  return (v1 + fact2*(v2-v1));
};

// coolprop functions

void CLookUpTable::CoolPropGetFluidState(string InputSpec, string fluid_name, double input1, double input2, struct FluidState &state){
    state.P=CoolProp::PropsSI("P", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.T=CoolProp::PropsSI("T", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name)-273.15	;
    state.d=CoolProp::PropsSI("D", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.v=1.0/state.d;
    state.h=CoolProp::PropsSI("H", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.s=CoolProp::PropsSI("S", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.u=CoolProp::PropsSI("U", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.cv=CoolProp::PropsSI("CVMASS", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.cp=CoolProp::PropsSI("CPMASS", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.c=CoolProp::PropsSI("A", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.eta=CoolProp::PropsSI("V", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.lambda=CoolProp::PropsSI("L", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.phase=CoolProp::PropsSI("PHASE", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.dlam_dT=CoolPropCalcDerivative("L", "T", "D", InputSpec, input1, input2, fluid_name);
    state.dlam_drho=CoolPropCalcDerivative("L", "D", "T", InputSpec, input1, input2, fluid_name);
    state.deta_dT=CoolPropCalcDerivative("V", "T", "D", InputSpec, input1, input2, fluid_name);
    state.deta_drho=CoolPropCalcDerivative("V", "D", "T", InputSpec, input1, input2, fluid_name);
    state.dPdrho_e=CoolProp::PropsSI("d(P)/d(D)|U", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.dPde_rho=CoolProp::PropsSI("d(P)/d(U)|D", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.dTdrho_e=CoolProp::PropsSI("d(T)/d(D)|U", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.dTde_rho=CoolProp::PropsSI("d(T)/d(U)|D", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.dhdrho_P=CoolProp::PropsSI("d(H)/d(D)|P", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.dhdP_rho=CoolProp::PropsSI("d(H)/d(P)|D", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.dsdrho_P=CoolProp::PropsSI("d(S)/d(D)|P", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
    state.dsdP_rho=CoolProp::PropsSI("d(S)/d(P)|D", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
}

double CLookUpTable::CoolPropCalcDerivative(string property, string withrespect_to, string at_constant, string InputSpec, double input1, double input2, string fluid_name){

    double fph, fmh;
    double eps=0.0001;
    if (string(1,InputSpec[0])==withrespect_to) {
        fph = CoolProp::PropsSI(property, string(1,InputSpec[0]), input1+.5*eps, string(1,InputSpec[1]), input2, fluid_name);
        fmh = CoolProp::PropsSI(property, string(1,InputSpec[0]), input1-.5*eps, string(1,InputSpec[1]), input2, fluid_name);
    }
    else if (string(1,InputSpec[1])==withrespect_to) {
        fph = CoolProp::PropsSI(property, string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2+.5*eps, fluid_name);
        fmh = CoolProp::PropsSI(property, string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2-.5*eps, fluid_name);
    }
    return (fph-fmh)/eps;
}

double CLookUpTable::CoolPropGetVaporQuality(string InputSpec, string fluid_name, double input1, double input2){

	return CoolProp::PropsSI("Q", string(1,InputSpec[0]), input1, string(1,InputSpec[1]), input2, fluid_name);
}


