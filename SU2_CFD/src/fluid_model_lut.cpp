/*!
 * fluid_model_lut.cpp
 * \brief Source of the look-up-table model.
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

#include "../include/fluid_model.hpp"

CLUTFluidModel::CLUTFluidModel() : CFluidModel() {

  Gamma = 0.0;
  Gamma_Minus_One = 0.0;
  Gas_Constant = 0.0;
  Cp = 0.0;

}


//CLUTFluidModel::CLUTFluidModel(string table_name,string fluid, string tab_dist, int table_imax, int table_jmax, string interpolation_scheme) : CFluidModel() {
//  Gamma = 0;
//  Gamma_Minus_One = Gamma - 1.0;
//  Gas_Constant = 0.0;
//  Cp = Gamma/Gamma_Minus_One*Gas_Constant;
//  LookUpTable = new CLookUpTable(table_name, fluid, tab_dist, table_imax, table_jmax, interpolation_scheme);
//
//}
CLUTFluidModel::CLUTFluidModel(CConfig *config) : CFluidModel() {
  Gamma = 0;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = 0.0;
  Cp = Gamma/Gamma_Minus_One*Gas_Constant;
  LookUpTable = new CLookUpTable(config);

}


CLUTFluidModel::~CLUTFluidModel(void) {

}

void CLUTFluidModel::SetTDState_rhoe (su2double rho, su2double e ) {
  Density = rho;
  StaticEnergy = e;
  Pressure = Gamma_Minus_One*Density*StaticEnergy;
  Temperature = Gamma_Minus_One*StaticEnergy/Gas_Constant;
  SoundSpeed2 = Gamma*Pressure/Density;
  Entropy = (1.0/Gamma_Minus_One*log(Temperature) + log(1.0/Density))*Gas_Constant;
  dPdrho_e = Gamma_Minus_One*StaticEnergy;
  dPde_rho = Gamma_Minus_One*Density;
  dTdrho_e = 0.0;
  dTde_rho = Gamma_Minus_One/Gas_Constant;
}

void CLUTFluidModel::SetTDState_PT (su2double P, su2double T ) {


}

void CLUTFluidModel::SetTDState_Prho (su2double P, su2double rho ) {


}

void CLUTFluidModel::SetEnergy_Prho (su2double P, su2double rho ) {

}

void CLUTFluidModel::SetTDState_hs (su2double h, su2double s ) {



}

void CLUTFluidModel::SetTDState_Ps (su2double P, su2double s ) {



}

void CLUTFluidModel::SetTDState_rhoT (su2double rho, su2double T ) {



}

void CLUTFluidModel::ComputeDerivativeNRBC_Prho(su2double P, su2double rho ){



}









