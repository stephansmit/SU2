/*!
 * transport_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Smit, R.Pecnik
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

#include "../include/transport_model.hpp"

/*-------------------------------------------------*/
/*--------------- Viscosity Model -----------------*/
/*-------------------------------------------------*/

CViscosityLUT::CViscosityLUT(void) : CViscosityModel() {

}

CViscosityLUT::CViscosityLUT(CLookUpTable* lookuptable) : CViscosityModel() {
	LookUpTable = lookuptable;
}

CViscosityLUT::~CViscosityLUT(void) { }


void CViscosityLUT::SetViscosity(su2double T, su2double rho) {
	if (LookUpTable->CheckIfInterpolated_rhoT(rho, T)){
		Mu = LookUpTable->GetLaminarViscosity();
	}
}

void CViscosityLUT::SetDerViscosity(su2double T, su2double rho) {
	if (LookUpTable->CheckIfInterpolated_rhoT(rho, T)){
	  dmudrho_T = LookUpTable->Getdmudrho_T();
	  dmudT_rho = LookUpTable->GetdmudT_rho();
	}
}

/*-------------------------------------------------*/
/*---------- Thermal Conductivity Models ----------*/
/*-------------------------------------------------*/


CConductivityLUT::CConductivityLUT(void) : CConductivityModel() {

}

CConductivityLUT::CConductivityLUT(CLookUpTable *lookuptable) : CConductivityModel() {
	LookUpTable = lookuptable;
}

void CConductivityLUT::SetConductivity(su2double T, su2double rho) {
	if (LookUpTable->CheckIfInterpolated_rhoT(rho, T)){
		Kt = LookUpTable->GetThermalConductivity();
	}
}

void CConductivityLUT::SetDerConductivity(su2double T, su2double rho) {
	if (LookUpTable->CheckIfInterpolated_rhoT(rho, T)){
		  dktdrho_T = LookUpTable->Getdktdrho_T();
		  dktdT_rho = LookUpTable->GetdktdT_rho();
	}
}

CConductivityLUT::~CConductivityLUT(void) { }


