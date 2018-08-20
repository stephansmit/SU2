/*!
 * \file look_up_table.hpp
 * \brief Headers of the main thermodynamic subroutines of the SU2 solvers.
 * \author S. Smit
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

#pragma once

#include "../../Common/include/mpi_structure.hpp"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>
#include "CoolProp.h"

#define LEN_COMPONENTS 32

#include "stdio.h"
#include "math.h"

#include "../../Common/include/config_structure.hpp"

using namespace std;


class CLookUpTable {
protected:

int imax, jmax;
string table_filename;
string table_fluid;
string table_distribution;
string table_interpolation_scheme;
bool create_table;
double rhomin, rhomax;
double Tmin, Tmax;

struct FluidState {
        double P;
        double T;
        double v;
        double d;
        double h;
        double s;
        double u;
        double cv;
        double cp;
        double c;
        double alpha;
        double beta;
        double zeta;
        double eta;
        double lambda;
        int i;
        int j;
        int phase;
        double dlam_dT;
        double dlam_drho;
        double deta_dT;
        double deta_drho;
        double dPdrho_e;
        double dPde_rho;
        double dTdrho_e;
        double dTde_rho;
        double dhdrho_P;
        double dhdP_rho;
        double dsdrho_P;
        double dsdP_rho;
};
struct TransferFluidState {
        double P;
        double d;
        double h;
        double s;
};
struct TransferFluidState2 {
        double P;
        double T;
        double d;
};

void CreateTableRhoT(void);
void CreateTablePs(void);
void CreateTablePT(void);


void ReadTableRhoT(string table_name);
void ReadTablePs(string table_name);
void ReadTablePT(string table_name);



// Table writing an reading

int ReadTableBIN(const char *fileName);
int ReadTransferTableBIN(const char *fileName);
int ReadTransferTable2BIN(const char *fileName);



void SaveTableBIN(const char *fileName);
void SaveTableTEC(const char *fileName);


void SaveTransferTableBIN(const char *fileName);
void SaveTransferTableTEC(const char *fileName);
void SaveTransferTable2BIN(const char *fileName);
void SaveTransferTable2TEC(const char *fileName);

// CoolProp functions

double CoolPropCalcDerivative( string property,
						string withrespect_to,
						string at_constant,
						string InputSpec,
						double input1,
						double input2,
						string fluid_name);

void CoolPropGetFluidState( string InputSpec,
							string fluid_name,
							double input1,
							double input2,
							struct FluidState &state);

void CoolPropGetTransferFluidStatePS( string InputSpec,
									string fluid_name,
									double input1,
									double input2,
									struct TransferFluidState &state);

void CoolPropGetTransferFluidStatePT( string InputSpec,
									string fluid_name,
									double input1,
									double input2,
									struct TransferFluidState2 &state);
double CoolPropGetVaporQuality( string InputSpec,
							string fluid_name,
							double input1,
							double input2);
// Interpolation functions


void InterpolateProperties(const char* InputSpec,
		double Input1,
		double Input2,
		struct FluidState *state, const char *props="ALL");

void bilinInterpolState(FluidState &resState,
		  const int i,
		  const int j,
          const double fact1,
		  const double fact2,
		  const char *props);
double bilinInterpolTransferState(
		  const int i,
		  const int j,
		  const double fact1,
		  const double fact2
		);

double bilinInterpolTransferState2(
		  const int i,
		  const int j,
		  const double fact1,
		  const double fact2
		);

double bilinInterpol(double a1, double a2, double b1, double b2, double fact1, double fact2);


public:
		FluidState **TableState;
		TransferFluidState **TransferTableState;
		TransferFluidState2 **TransferTableState2;
		FluidState *OutputState;

	     /*!
		 * \brief Constructor of the class.
		 */
		CLookUpTable(void);

		/*!
		 * \brief Constructor of the class.
		 */
		CLookUpTable(string table_name,
					 string fluid,
					 string tab_dist,
					 int table_imax,
					 int table_jmax,
					 double rho_min,
					 double rho_max,
					 double T_min,
					 double T_max,
					 string interpolation_scheme,
					 bool createTable);


		/*!
		 * \brief Constructor of the class.
		 */
		CLookUpTable(CConfig *config);


		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CLookUpTable(void);


		/*!
		 * \brief Get fluid pressure.
		 */
		su2double GetPressure ();

		/*!
		 * \brief Get fluid temperature.
		 */
		su2double GetTemperature ();

		/*!
		 * \brief Get fluid entropy.
		 */
		su2double GetEntropy ();

		/*!
		 * \brief Get fluid internal energy.
		 */
		su2double GetStaticEnergy ();

		/*!
		 * \brief Get fluid density.
		 */
		su2double GetDensity ();

		/*!
		 * \brief Get fluid speed of sound.
		 */
		su2double GetSoundSpeed ();

		/*!
		 * \brief Get fluid speed of sound squared.
		 */
		su2double GetSoundSpeed2 ();

		/*!
		 * \brief Get fluid specific heat at constant pressure.
		 */
		su2double GetCp ();

		/*!
		 * \brief Get fluid dynamic viscosity
		 */

		su2double GetLaminarViscosity ();

		/*!
		 * \brief Get fluid thermal conductivity
		 */

		su2double GetThermalConductivity ();

		/*!
		 * \brief Get fluid pressure partial derivative.
		 */
		su2double GetdPdrho_e ();

		/*!
		 * \brief Get fluid pressure partial derivative.
		 */
		su2double GetdPde_rho ();

		/*!
		 * \brief Get fluid temperature partial derivative.
		 */
		su2double GetdTdrho_e ();

		/*!
		 * \brief Get fluid temperature partial derivative.
		 */
		su2double GetdTde_rho ();

		/*!
		 * \brief Get fluid pressure partial derivative.
		 */
		su2double Getdhdrho_P ();

		/*!
		 * \brief Get fluid pressure partial derivative.
		 */
		su2double GetdhdP_rho ();

		/*!
		 * \brief Get fluid temperature partial derivative.
		 */
		su2double Getdsdrho_P ();

		/*!
		 * \brief Get fluid temperature partial derivative.
		 */
		su2double GetdsdP_rho ();

		/*!
		 * \brief Get fluid dynamic viscosity partial derivative.
		 */
		su2double Getdmudrho_T ();

		/*!
		 * \brief Get fluid dynamic viscosity partial derivative.
		 */
		su2double GetdmudT_rho ();

		/*!
		 * \brief Get fluid thermal conductivity partial derivative.
		 */
		su2double Getdktdrho_T ();

		/*!
		 * \brief Get fluid thermal conductivity partial derivative.
		 */
		su2double GetdktdT_rho ();


		void SetTDState_rhoe (su2double rho, su2double e );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("PT").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (T).
		 */

		void SetTDState_PT (su2double P, su2double T );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 */

		void SetTDState_Prho (su2double P, su2double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 *
		 */

		void SetEnergy_Prho (su2double P, su2double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("hs").
		 * \param[in] th1 - first thermodynamic variable (h).
		 * \param[in] th2 - second thermodynamic variable (s).
		 *
		 */
		void SetTDState_hs (su2double h, su2double s );


		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (T).
		 *
		 */
		void SetTDState_rhoT (su2double rho, su2double T );


		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (s).
		 */

		void SetTDState_Ps (su2double P, su2double s );


		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 *
		 */
		void ComputeDerivativeNRBC_Prho (su2double P, su2double rho );


		bool CheckIfInterpolated_rhoe(su2double rho, su2double e );

		bool CheckIfInterpolated_PT(su2double P, su2double T );

		bool CheckIfInterpolated_Prho(su2double P, su2double rho );

		bool CheckIfInterpolated_hs(su2double h, su2double s );

		bool CheckIfInterpolated_rhoT(su2double rho, su2double T );

		bool CheckIfInterpolated_Ps(su2double P, su2double s );


};

#include "look_up_table.inl"
