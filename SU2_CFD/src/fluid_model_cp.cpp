/*!
 * fluid_model_pig.cpp
 * \brief Source of the ideal gas model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
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

CCPFluidModel::CCPFluidModel() : CFluidModel() {

  Gamma = 0.0;
  Gamma_Minus_One = 0.0;
  Gas_Constant = 0.0;
  Cp = 0.0;
}


CCPFluidModel::CCPFluidModel(string fluid_name ) : CFluidModel() {
  FluidName = fluid_name;
}


CCPFluidModel::~CCPFluidModel(void) {

}

void CCPFluidModel::SetTDState_rhoe (su2double rho, su2double e ) {
	StaticEnergy = 	CoolProp::PropsSI("U","D",rho,"U",e,FluidName);
	Entropy = 		CoolProp::PropsSI("S","D",rho,"U",e,FluidName);
	Density = 		CoolProp::PropsSI("D","D",rho,"U",e,FluidName);
	Pressure = 		CoolProp::PropsSI("P","D",rho,"U",e,FluidName);
	SoundSpeed2 = 	CoolProp::PropsSI("A","D",rho,"U",e,FluidName);
	SoundSpeed2 *= SoundSpeed2;
	Temperature = 	CoolProp::PropsSI("T","D",rho,"U",e,FluidName);
	dPdrho_e = CoolProp::PropsSI("d(P)/d(D)|U","D",rho,"U",e,FluidName);
	dPde_rho = CoolProp::PropsSI("d(P)/d(U)|D","D",rho,"U",e,FluidName);
	dTdrho_e = CoolProp::PropsSI("d(T)/d(D)|U","D",rho,"U",e,FluidName);
	dTde_rho = CoolProp::PropsSI("d(T)/d(U)|D","D",rho,"U",e,FluidName);
}

void CCPFluidModel::SetTDState_PT (su2double P, su2double T ) {
	StaticEnergy = 	CoolProp::PropsSI("U","P",P,"T",T,FluidName);
	Entropy = 		CoolProp::PropsSI("S","P",P,"T",T,FluidName);
	Density = 		CoolProp::PropsSI("D","P",P,"T",T,FluidName);
	Pressure = 		CoolProp::PropsSI("P","P",P,"T",T,FluidName);
	SoundSpeed2 = 	CoolProp::PropsSI("A","P",P,"T",T,FluidName);
	SoundSpeed2 *= SoundSpeed2;
	Temperature = 	CoolProp::PropsSI("T","P",P,"T",T,FluidName);
	dPdrho_e = CoolProp::PropsSI("d(P)/d(D)|U","P",P,"T",T,FluidName);
	dPde_rho = CoolProp::PropsSI("d(P)/d(U)|D","P",P,"T",T,FluidName);
	dTdrho_e = CoolProp::PropsSI("d(T)/d(D)|U","P",P,"T",T,FluidName);
	dTde_rho = CoolProp::PropsSI("d(T)/d(U)|D","P",P,"T",T,FluidName);
}

void CCPFluidModel::SetTDState_Prho (su2double P, su2double rho ) {
	StaticEnergy = 	CoolProp::PropsSI("U","P",P,"D",rho,FluidName);
	Entropy = 		CoolProp::PropsSI("S","P",P,"D",rho,FluidName);
	Density = 		CoolProp::PropsSI("D","P",P,"D",rho,FluidName);
	Pressure = 		CoolProp::PropsSI("P","P",P,"D",rho,FluidName);
	SoundSpeed2 = 	CoolProp::PropsSI("A","P",P,"D",rho,FluidName);
	SoundSpeed2 *= SoundSpeed2;
	Temperature = 	CoolProp::PropsSI("T","P",P,"D",rho,FluidName);
	dPdrho_e = CoolProp::PropsSI("d(P)/d(D)|U","P",P,"D",rho,FluidName);
	dPde_rho = CoolProp::PropsSI("d(P)/d(U)|D","P",P,"D",rho,FluidName);
	dTdrho_e = CoolProp::PropsSI("d(T)/d(D)|U","P",P,"D",rho,FluidName);
	dTde_rho = CoolProp::PropsSI("d(T)/d(U)|D","P",P,"D",rho,FluidName);
}

void CCPFluidModel::SetEnergy_Prho (su2double P, su2double rho ) {
	StaticEnergy = 	CoolProp::PropsSI("U","P",P,"D",rho,FluidName);
}

void CCPFluidModel::SetTDState_hs (su2double h, su2double s ) {
	StaticEnergy = 	CoolProp::PropsSI("U","H",h,"S",s,FluidName);
	Entropy = 		CoolProp::PropsSI("S","H",h,"S",s,FluidName);
	Density = 		CoolProp::PropsSI("D","H",h,"S",s,FluidName);
	Pressure =		CoolProp::PropsSI("P","H",h,"S",s,FluidName);
	SoundSpeed2 = 	CoolProp::PropsSI("A","H",h,"S",s,FluidName);
	SoundSpeed2 *= SoundSpeed2;
	Temperature = 	CoolProp::PropsSI("T","H",h,"S",s,FluidName);
	dPdrho_e = CoolProp::PropsSI("d(P)/d(D)|U","H",h,"S",s,FluidName);
	dPde_rho = CoolProp::PropsSI("d(P)/d(U)|D","H",h,"S",s,FluidName);
	dTdrho_e = CoolProp::PropsSI("d(T)/d(D)|U","H",h,"S",s,FluidName);
	dTde_rho = CoolProp::PropsSI("d(T)/d(U)|D","H",h,"S",s,FluidName);
}

void CCPFluidModel::SetTDState_Ps (su2double P, su2double s ) {
	StaticEnergy = 	CoolProp::PropsSI("U","P",P,"S",s,FluidName);
	Entropy = 		CoolProp::PropsSI("S","P",P,"S",s,FluidName);
	Density = 		CoolProp::PropsSI("D","P",P,"S",s,FluidName);
	Pressure = 		CoolProp::PropsSI("P","P",P,"S",s,FluidName);
	SoundSpeed2 = 	CoolProp::PropsSI("A","P",P,"S",s,FluidName);
	SoundSpeed2 *= SoundSpeed2;
	Temperature = 	CoolProp::PropsSI("T","P",P,"S",s,FluidName);
	dPdrho_e = CoolProp::PropsSI("d(P)/d(D)|U","P",P,"S",s,FluidName);
	dPde_rho = CoolProp::PropsSI("d(P)/d(U)|D","P",P,"S",s,FluidName);
	dTdrho_e = CoolProp::PropsSI("d(T)/d(D)|U","P",P,"S",s,FluidName);
	dTde_rho = CoolProp::PropsSI("d(T)/d(U)|D","P",P,"S",s,FluidName);
}

void CCPFluidModel::SetTDState_rhoT (su2double rho, su2double T ) {
	StaticEnergy = 	CoolProp::PropsSI("U","D",rho,"T",T,FluidName);
	Entropy = 		CoolProp::PropsSI("S","D",rho,"T",T,FluidName);
	Density = 		CoolProp::PropsSI("D","D",rho,"T",T,FluidName);
	Pressure = 		CoolProp::PropsSI("P","D",rho,"T",T,FluidName);
	SoundSpeed2 = 	CoolProp::PropsSI("A","D",rho,"T",T,FluidName);
	SoundSpeed2 *= SoundSpeed2;
	Temperature = 	CoolProp::PropsSI("T","D",rho,"T",T,FluidName);
	dPdrho_e = CoolProp::PropsSI("d(P)/d(D)|U","D",rho,"T",T,FluidName);
	dPde_rho = CoolProp::PropsSI("d(P)/d(U)|D","D",rho,"T",T,FluidName);
	dTdrho_e = CoolProp::PropsSI("d(T)/d(D)|U","D",rho,"T",T,FluidName);
	dTde_rho = CoolProp::PropsSI("d(T)/d(U)|D","D",rho,"T",T,FluidName);
}

void CCPFluidModel::ComputeDerivativeNRBC_Prho(su2double P, su2double rho ){
	dhdrho_P = CoolProp::PropsSI("d(H)/d(D)|P","P",P,"D",rho,FluidName);
	dhdP_rho = CoolProp::PropsSI("d(H)/d(P)|D","P",P,"D",rho,FluidName);
	dsdrho_P = CoolProp::PropsSI("d(S)/d(D)|P","P",P,"D",rho,FluidName);
	dsdP_rho = CoolProp::PropsSI("d(S)/d(P)|D","P",P,"D",rho,FluidName);
}









