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
#include "CoolProp.h"
/*-------------------------------------------------*/
/*--------------- Viscosity Model -----------------*/
/*-------------------------------------------------*/

CViscosityCP::CViscosityCP(void) : CViscosityModel() {

}

CViscosityCP::CViscosityCP(string fluid_name) : CViscosityModel() {
	FluidName  = fluid_name;
}

CViscosityCP::~CViscosityCP(void) { }


void CViscosityCP::SetViscosity(su2double T, su2double rho) {
	Mu = CoolProp::PropsSI("V", "T", T, "D", rho, FluidName );
}

void CViscosityCP::SetDerViscosity(su2double T, su2double rho) {
	dmudrho_T = 0;
	dmudT_rho = 0;
}

/*-------------------------------------------------*/
/*---------- Thermal Conductivity Models ----------*/
/*-------------------------------------------------*/


CConductivityCP::CConductivityCP(void) : CConductivityModel() {

}

CConductivityCP::CConductivityCP(string fluid_name) : CConductivityModel() {
	FluidName = fluid_name;
}

void CConductivityCP::SetConductivity(su2double T, su2double rho) {
	Kt = CoolProp::PropsSI("L", "T", T, "D", rho, FluidName );
}

void CConductivityCP::SetDerConductivity(su2double T, su2double rho) {
	dktdrho_T = 0;
	dktdT_rho = 0;

}

CConductivityCP::~CConductivityCP(void) { }


