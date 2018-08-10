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

}

CLookUpTable::CLookUpTable(string table_name,string table_dist, int column1_size, int column2_size){
	imax = column1_size;
	jmax = column2_size;
	tab_dist = table_dist;
	OutputState = new FluidState;
	ReadTableRhoT(table_name);
}

CLookUpTable::~CLookUpTable(void) {

}

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


void CLookUpTable::ReadTableRhoT(string filename){
	TableState = new FluidState*[imax];
    for (int i=0; i<imax; i++)
    	TableState[i] = new FluidState[jmax];
    ReadTableBIN(filename.c_str());
}


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

void CLookUpTable::SetTDState_rhoe (su2double rho, su2double e) {}


void CLookUpTable::SetTDState_PT (su2double P, su2double T ){};


void CLookUpTable::SetTDState_Prho (su2double P, su2double rho ){};


void CLookUpTable::SetEnergy_Prho (su2double P, su2double rho ){};


void CLookUpTable::SetTDState_hs (su2double h, su2double s ){};


void CLookUpTable::SetTDState_rhoT (su2double rho, su2double T ){
	if (not CheckIfInterpolated_rhoT(rho,T)) fluidprop_allprops_("Tv",(double) T, (double)1.0/rho, OutputState, "ALL" );
};


void CLookUpTable::SetTDState_Ps (su2double P, su2double s ){};

//void CLookUpTable::ComputeDerivativeNRBC_Prho (su2double P, su2double rho );

void CLookUpTable::fluidprop_allprops_( const char* InputSpec,
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
      if (strcmp(tab_dist.c_str(), "LOG")==0) ff = (log10(dens)-log10(TableState[0][0].d))/(log10(TableState[0][jmax-1].d)-log10(TableState[0][0].d));
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
        if (strcmp(tab_dist.c_str(), "LOG")==0) ff = (log10(dens)-log10(TableState[0][0].d))/(log10(TableState[0][jmax-1].d)-log10(TableState[0][0].d));
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
