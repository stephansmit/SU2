/*!
 * \file look_up_table.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
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

inline su2double CLookUpTable::GetPressure () { return OutputState->P; }
inline su2double CLookUpTable::GetSoundSpeed () { return OutputState->c; }
inline su2double CLookUpTable::GetSoundSpeed2 () { return OutputState->c*OutputState->c; }
inline su2double CLookUpTable::GetDensity () { return OutputState->d; }
inline su2double CLookUpTable::GetEntropy () { return OutputState->s; }
inline su2double CLookUpTable::GetStaticEnergy () { return OutputState->u; }
inline su2double CLookUpTable::GetTemperature () { return OutputState->T; }
inline su2double CLookUpTable::GetCp () { return OutputState->cp; }
inline su2double CLookUpTable::GetdPdrho_e () { return OutputState->dPdrho_e; }
inline su2double CLookUpTable::GetdPde_rho () { return OutputState->dPde_rho; }
inline su2double CLookUpTable::GetdTdrho_e () { return OutputState->dTdrho_e; }
inline su2double CLookUpTable::GetdTde_rho () { return OutputState->dTde_rho; }
inline su2double CLookUpTable::Getdhdrho_P () {return OutputState->dhdrho_P;}
inline su2double CLookUpTable::GetdhdP_rho () {return OutputState->dhdP_rho;}
inline su2double CLookUpTable::Getdsdrho_P () {return OutputState->dsdrho_P;}
inline su2double CLookUpTable::GetdsdP_rho () {return OutputState->dsdP_rho;}
inline su2double CLookUpTable::GetLaminarViscosity (){ return OutputState->eta;} 
inline su2double CLookUpTable::Getdmudrho_T () { return OutputState->deta_drho;}
inline su2double CLookUpTable::GetdmudT_rho () { return OutputState->deta_dT;}
inline su2double CLookUpTable::GetThermalConductivity () {return OutputState->lambda;}
inline su2double CLookUpTable::Getdktdrho_T () { return OutputState->dlam_drho;}
inline su2double CLookUpTable::GetdktdT_rho () { return OutputState->dlam_dT;}

