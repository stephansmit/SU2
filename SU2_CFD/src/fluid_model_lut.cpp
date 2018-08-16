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


CLUTFluidModel::CLUTFluidModel(string table_name,
							   string fluid,
							   string tab_dist,
							   int table_imax,
							   int table_jmax,
							   double rho_min,
							   double rho_max,
							   double T_min,
							   double T_max,
							   string interpolation_scheme,
							   bool createTable) : CFluidModel() {
  Gamma = 0;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = 0.0;
  Cp = Gamma/Gamma_Minus_One*Gas_Constant;
  LookUpTable = new CLookUpTable(
		  table_name,
		  fluid,
		  tab_dist,
		  table_imax,
		  table_jmax,
		  rho_min,
		  rho_max,
		  T_min,
		  T_max,
		  interpolation_scheme,
		  createTable);
}


CLUTFluidModel::CLUTFluidModel(CConfig *config) : CFluidModel() {
  Gamma = 0;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = 0.0;
  Cp = Gamma/Gamma_Minus_One*Gas_Constant;
  LookUpTable = new CLookUpTable(config);
//  this->SetThermalConductivityModel(config);
//  this->SetLaminarViscosityModel(config);
}


CLUTFluidModel::~CLUTFluidModel(void) {

}

void CLUTFluidModel::SetTDState_rhoe (su2double rho, su2double e ) {
//	cout << "SetTDState_rhoe" << endl;
	LookUpTable->SetTDState_rhoe(rho,e);
	CopyStateLUT();
}

void CLUTFluidModel::SetTDState_PT (su2double P, su2double T ) {
//	cout << "SetTDState_PT" << endl;
	LookUpTable->SetTDState_PT(P,T);
	CopyStateLUT();
}

void CLUTFluidModel::SetTDState_Prho (su2double P, su2double rho ) {
//	cout << "SetTDState_Prho" << endl;
	LookUpTable->SetTDState_Prho(P,rho);
	CopyStateLUT();
}

void CLUTFluidModel::SetEnergy_Prho (su2double P, su2double rho ) {
//	cout << "SetEnergy_Prho" << endl;
	LookUpTable->SetEnergy_Prho(P,rho);
	StaticEnergy = LookUpTable->GetStaticEnergy();
}

void CLUTFluidModel::SetTDState_hs (su2double h, su2double s ) {
//	cout << "SetTDState_hs" << endl;
	LookUpTable->SetTDState_hs(h,s);
	CopyStateLUT();
}

void CLUTFluidModel::SetTDState_Ps (su2double P, su2double s ) {
//	cout << "SetTDState_Ps" << endl;
	LookUpTable->SetTDState_Ps(P,s);
	CopyStateLUT();
}

void CLUTFluidModel::SetTDState_rhoT (su2double rho, su2double T ) {
//	cout << "SetTDState_rhoT" << endl;
	LookUpTable->SetTDState_rhoT(rho,T);

	CopyStateLUT();
}

void CLUTFluidModel::ComputeDerivativeNRBC_Prho(su2double P, su2double rho ){
	//	cout << "SetTDState_rhoT" << endl;
	LookUpTable->SetTDState_Prho(P, rho);
	CopyStateLUT();
}


void CLUTFluidModel::CopyStateLUT(){
//	if (isnan(LookUpTable->GetStaticEnergy())) cout << "GetStaticEnergy" <<endl;
//	if (isnan(LookUpTable->GetEntropy())) cout << "GetEntropy" <<endl;
//	if (isnan(LookUpTable->GetDensity())) cout << "GetDensity" <<endl;
//	if (isnan(LookUpTable->GetPressure())) cout << "GetPressure" <<endl;
//	if (isnan(LookUpTable->GetSoundSpeed2())) cout << "GetSoundSpeed2" <<endl;
//	if (isnan(LookUpTable->GetTemperature())) cout << "GetTemperature" <<endl;
//	if (isnan(LookUpTable->GetdPdrho_e())) cout << "GetdPdrho_e" <<endl;
//	if (isnan(LookUpTable->GetdPde_rho())) cout << "GetdPde_rho" <<endl;
//	if (isnan(LookUpTable->GetdTdrho_e())) cout << "GetdTdrho_e" <<endl;
//	if (isnan(LookUpTable->GetdTde_rho())) cout << "GetdTde_rho" <<endl;
//	if (isnan(LookUpTable->Getdhdrho_P())) cout << "Getdhdrho_P" <<endl;
//	if (isnan(LookUpTable->GetdhdP_rho())) cout << "GetdhdP_rho" <<endl;
//	if (isnan(LookUpTable->Getdsdrho_P())) cout << "Getdsdrho_P" <<endl;
//	if (isnan(LookUpTable->GetdsdP_rho())) cout << "GetdsdP_rho" <<endl;
//	if (isnan(LookUpTable->GetCp())) cout << "GetCp" <<endl;
//	if (isnan(LookUpTable->GetLaminarViscosity())) cout << "GetLaminarViscosity" <<endl;
//	if (isnan(LookUpTable->Getdmudrho_T())) cout << "Getdmudrho_T" <<endl;
//	if (isnan(LookUpTable->GetdmudT_rho())) cout << "GetdmudT_rho" <<endl;
//	if (isnan(LookUpTable->GetThermalConductivity())) cout << "GetThermalConductivity" <<endl;
//	if (isnan(LookUpTable->Getdktdrho_T())) cout << "Getdktdrho_T" <<endl;
//	if (isnan(LookUpTable->GetdktdT_rho())) cout << "GetdktdT_rho" <<endl;

	StaticEnergy = LookUpTable->GetStaticEnergy();
    Entropy = LookUpTable->GetEntropy();
    Density = LookUpTable->GetDensity();
    Pressure = LookUpTable->GetPressure();
    SoundSpeed2 = LookUpTable->GetSoundSpeed2();
    Temperature = LookUpTable->GetTemperature();
    dPdrho_e = LookUpTable->GetdPdrho_e();
    dPde_rho = LookUpTable->GetdPde_rho();
    dTdrho_e = LookUpTable->GetdTdrho_e();
    dTde_rho = LookUpTable->GetdTde_rho();
    dhdrho_P = LookUpTable->Getdhdrho_P();
    dhdP_rho = LookUpTable->GetdhdP_rho();
    dsdrho_P = LookUpTable->Getdsdrho_P();
    dsdP_rho = LookUpTable->GetdsdP_rho();
    Cp = LookUpTable->GetCp();
    Mu = LookUpTable->GetLaminarViscosity();
    dmudrho_T = LookUpTable->Getdmudrho_T();
    dmudT_rho = LookUpTable->GetdmudT_rho();
    Kt = LookUpTable->GetThermalConductivity();
    dktdrho_T = LookUpTable->Getdktdrho_T();
    dktdT_rho = LookUpTable->GetdktdT_rho();


//	dhdrho_P= LookUpTable->
//	dhdP_rho= 1.0/dPde_rho +1.0/rho;
//	dPds_rho= rho*rho*(SoundSpeed2 - dPdrho_T)/dPdT_rho;
//	dsdP_rho= 1.0/dPds_rho;
//	dsdrho_P= -SoundSpeed2/dPds_rho;
}






