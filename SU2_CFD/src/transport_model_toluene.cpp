/*!
 * transport_model_toluene.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Smit
 * \version 6.0.1 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

CViscosityToluene::CViscosityToluene(void) : CViscosityModel() {
  alpha_0 = 0.401080;  alpha_1 = -0.476409;  alpha_2 = 0.0;  alpha_3 = 0.069442;
  b_0 =-19.572881; b_1 = 219.73999;  b_2 =-1015.3226;  b_3 =  2471.0125; b_4 = -3375.1717;
  b_5 = 2491.6597; b_6 = -787.26086; b_7 = 14.085455;  b_8 =-0.34664158;
  c_0 = 19.919216; c_1 = -2.6557905; c_2 =-135.904211; c_3 = -7.9962719; c_4 = -11.014795; c_5 = -10.113817;
  M = 92.13842;
  sigma = 0.524;
  EnergyScale = 472.0;
  N_A = 6.0221409e+23;
  T_c=591.75;
  rho_c=291.987;
}

CViscosityToluene::~CViscosityToluene(void) { }

su2double CViscosityToluene::Mu_0(su2double T) {
  su2double Mu_0;
  Mu_0 = 0.021357*(pow((M*T),0.5)/(pow(sigma,2.0)*Sstar_mu(T)));
  return Mu_0;
}
su2double CViscosityToluene::dMu_0dT_rho(su2double T) {
  su2double dMu_0dT_rho;
  dMu_0dT_rho = ((0.021357*0.5)*M*Sstar_mu(T)-0.021357*M*T*dSstar_mudT_rho(T))/(pow(sigma,2.0)*pow(M*T,0.5)*pow(Sstar_mu(T),2.0));
  return dMu_0dT_rho;
}

su2double CViscosityToluene::Sstar_mu(su2double T){
  su2double Sstar_mu,Tstar;
  Tstar = T/EnergyScale;
  Sstar_mu = exp(alpha_0*pow(log(Tstar),0.0)+alpha_1*pow(log(Tstar),1.0)+alpha_2*pow(log(Tstar),2.0)+alpha_3*pow(log(Tstar),3.0));
  return Sstar_mu;
}

su2double CViscosityToluene::dSstar_mudT_rho(su2double T){
  su2double dSstar_mudT_rho,Tstar;
  Tstar = T/EnergyScale;
  dSstar_mudT_rho = ( (alpha_1/T) + (2*alpha_2*log(Tstar)/T) + (3*alpha_3*pow(log(Tstar),2.0)/T) )*Sstar_mu(T);
  return dSstar_mudT_rho;
}

su2double CViscosityToluene::Mu_1(su2double T) {
  su2double Mu_1;
  Mu_1=Mu_0(T)*B_mu(T);
  return Mu_1;
}

su2double CViscosityToluene::dMu_1dT_rho(su2double T) {
  su2double dMu_1dT_rho;
  dMu_1dT_rho=dMu_0dT_rho(T)*B_mu(T)+Mu_0(T)*dB_mudT_rho(T);
  return dMu_1dT_rho;
}

su2double CViscosityToluene::B_mu(su2double T) {
  su2double B_mu;
  B_mu = Bstar_mu(T)*((N_A*pow(sigma*1e-9,3.0))/(M*1e-3));
  return B_mu;
}

su2double CViscosityToluene::dB_mudT_rho(su2double T) {
  su2double dB_mudT_rho;
  dB_mudT_rho = dBstar_mudT_rho(T)*((N_A*pow(sigma*1e-9,3.0))/(M*1e-3));
  return dB_mudT_rho;
}
su2double CViscosityToluene::Bstar_mu(su2double T) {
  su2double Bstar_mu, Tstar;
  Tstar = T/EnergyScale;
  Bstar_mu =      b_0*pow(Tstar, (-0.25*0.0))+
		  b_1*pow(Tstar, (-0.25*1.0))+
		  b_2*pow(Tstar, (-0.25*2.0))+
		  b_3*pow(Tstar, (-0.25*3.0))+
		  b_4*pow(Tstar, (-0.25*4.0))+
		  b_5*pow(Tstar, (-0.25*5.0))+
		  b_6*pow(Tstar, (-0.25*6.0))+
		  b_7*pow(Tstar, -2.5)+
		  b_8*pow(Tstar, -5.5);
  return Bstar_mu;
}

su2double CViscosityToluene::dBstar_mudT_rho(su2double T) {
  su2double dBstar_mudT_rho, Tstar;
  Tstar = T/EnergyScale;
  dBstar_mudT_rho =
	          (-0.25*1.0*b_1/EnergyScale)*pow(Tstar, (-0.25*1.0-1.0))+
	          (-0.25*2.0*b_2/EnergyScale)*pow(Tstar, (-0.25*2.0-1.0))+
	          (-0.25*3.0*b_3/EnergyScale)*pow(Tstar, (-0.25*3.0-1.0))+
	          (-0.25*4.0*b_4/EnergyScale)*pow(Tstar, (-0.25*4.0-1.0))+
	          (-0.25*5.0*b_5/EnergyScale)*pow(Tstar, (-0.25*5.0-1.0))+
	          (-0.25*6.0*b_6/EnergyScale)*pow(Tstar, (-0.25*6.0-1.0))+
		  (-2.5*b_7/EnergyScale)*pow(Tstar, (-2.5-1.0))+
		  (-5.5*b_8/EnergyScale)*pow(Tstar, (-5.5-1.0));
  return dBstar_mudT_rho;
}

su2double CViscosityToluene::DeltaMu(su2double T, su2double rho) {
  su2double DeltaMu, T_r, rho_r;
  T_r=T/T_c;
  rho_r=rho/rho_c;
  DeltaMu = (pow(rho_r,2.0/3.0)*pow(T_r,0.5))*
		  	  	  (((c_0*rho_r+c_1*pow(rho_r,4.0))/T_r)+
				  ((c_2*pow(rho_r,3.0))/(pow(rho_r,2.0)+c_3+c_4*T_r))+
				  (c_5*rho_r));
  return DeltaMu;

}

su2double CViscosityToluene::dDeltaMudT_rho(su2double T, su2double rho){
  su2double dDeltaMudT_rho, term1, term2;
  su2double T_r, rho_r;
  T_r=T/T_c;
  rho_r=rho/rho_c;
  term1 = (0.5/T_c)*(pow(rho_r,2.0/3.0)*pow(T_r,-0.5))*
		  	  	  (
		  	  			  ((c_0*rho_r+c_1*pow(rho_r,4.0))/T_r)+
		  	  			  ((c_2*pow(rho_r,3.0))/(pow(rho_r,2.0)+c_3+c_4*T_r))+
		  	  			  (c_5*rho_r)
		  	  	  );

  term2=(pow(rho_r,2.0/3.0)*pow(T_r,0.5))*
		  (
				  -(T_c*(c_1*pow(rho_r,4.0)+c_0*rho_r)/pow(T,2.0))
				  -(c_2*c_4*pow(rho_r,3.0)/(T_c*pow((c_4*T_r+c_3+pow(rho_r,2.0)),2.0)))
		  );

  dDeltaMudT_rho=term1+term2;
  return dDeltaMudT_rho;
}

su2double CViscosityToluene::dDeltaMudrho_T(su2double T, su2double rho){
  su2double dDeltaMudrho_T, term1, term2;
  su2double T_r, rho_r;
  T_r=T/T_c;
  rho_r=rho/rho_c;
  term1 = (2.0/3.0)*(1/rho_c)*(pow(rho_r,-1.0/3.0)*pow(T_r,0.5))*
		  	  	  (
		  	  	      ((c_0*rho_r+c_1*pow(rho_r,4.0))/T_r)+
				      ((c_2*pow(rho_r,3.0))/(pow(rho_r,2.0)+c_3+c_4*T_r))+
				      (c_5*rho_r)
				  );
  term2 = pow(rho_r,2.0/3.0)*pow(T_r,0.5)*
		  	  	  (
					  (((c_0/rho_c)+(4.0/rho_c)*c_1*pow(rho_r,3.0))/T_r)-
					  (2*c_2*pow(rho,4.0)/(pow(rho_c,5.0)*pow(c_4*T_r+c_3+pow(rho_r,2.0),2.0)))+
					  (((3.0/rho_c)*c_2*pow(rho_r,2.0))/(pow(rho_r,2.0)+c_3+c_4*T_r))+
					  (c_5/rho_c)
		  	  	  );
  dDeltaMudrho_T = term1+term2;
  return dDeltaMudrho_T;
}
////(pow(rho/r,2.0/3.0)*pow(T/t,0.5))*(((c_0*(rho/r)+c_1*pow((rho/r),4.0))/(T/t))+((c_2*pow((rho/r),3.0))/(pow((rho/r),2.0)+c_3+c_4*(T/t)))+(c_5*(rho/r)))

void CViscosityToluene::SetViscosity(su2double T, su2double rho) {
  Mu = (Mu_0(T) + Mu_1(T)*rho + DeltaMu(T,rho))*1e-6;
}

void CViscosityToluene::SetDerViscosity(su2double T, su2double rho) {
  dmudrho_T = (Mu_1(T)+ dDeltaMudrho_T(T, rho))*1e-6;
  dmudT_rho = (dMu_0dT_rho(T)+dMu_1dT_rho(T)*rho+dDeltaMudT_rho(T,rho))*1e-6;
}

/*-------------------------------------------------*/
/*---------- Thermal Conductivity Model- ----------*/
/*-------------------------------------------------*/

CConductivityToluene::CConductivityToluene(void) : CConductivityModel() {
  T_c=591.75;
  rho_c=291.992;
  B_11 = -5.18530e-2, B_21 =  5.17449e-2;
  B_12 =  1.33846e-1, B_22 = -1.21902e-1;
  B_13 = -1.20446e-1, B_23 =  1.37748e-1;
  B_14 =  5.30211e-2, B_24 = -7.32792e-2;
  B_15 = -1.00604e-2, B_25 =  1.72914e-2;
  B_16 =  6.33457e-4, B_26 = -1.38585e-3;
  C_1=0.20e-3, C_2=4.50e-2, C_3=0.090;
}
su2double CConductivityToluene::Kt_0(su2double T) {
  su2double Kt_0;
  Kt_0 =  5.8808-6.1693e-2*T   +
		  3.4151e-4*pow(T,2.0) -
  	  	  3.0420e-7*pow(T,3.0) +
  	  	  1.2868e-10*pow(T,4.0)-
  	  	  2.1303e-14*pow(T,5.0);
  return Kt_0*1e-3;
}

su2double CConductivityToluene::dKt_0dT_rho(su2double T) {
  su2double dKt_0dT_rho;
  dKt_0dT_rho =  -6.1693e-2        			+
				  2.0*3.4151e-4*T			-
				  3.0*3.0420e-7*pow(T,2.0)	+
				  4.0*1.2868e-10*pow(T,3.0)	-
				  5.0*2.1303e-14*pow(T,4.0);
  return dKt_0dT_rho*1e-3;
}

su2double CConductivityToluene::DeltaKt(su2double T,su2double rho) {
  su2double DeltaKt;
  DeltaKt = (B_11+B_21*(T/T_c))*pow(rho/rho_c,1.0)+
		  	(B_12+B_22*(T/T_c))*pow(rho/rho_c,2.0)+
		  	(B_13+B_23*(T/T_c))*pow(rho/rho_c,3.0)+
		  	(B_14+B_24*(T/T_c))*pow(rho/rho_c,4.0)+
		  	(B_15+B_25*(T/T_c))*pow(rho/rho_c,5.0)+
		  	(B_16+B_26*(T/T_c))*pow(rho/rho_c,6.0);
  return DeltaKt;
}

su2double CConductivityToluene::dDeltaKtdrho_T(su2double T,su2double rho) {
  su2double dDeltaKtdrho_T;
  dDeltaKtdrho_T =  (B_11+B_21*(T/T_c))*(1.0/rho_c)+
					(B_12+B_22*(T/T_c))*(2.0/rho_c)*(rho/rho_c)+
					(B_13+B_23*(T/T_c))*(3.0/rho_c)*pow(rho/rho_c,2.0)+
					(B_14+B_24*(T/T_c))*(4.0/rho_c)*pow(rho/rho_c,3.0)+
					(B_15+B_25*(T/T_c))*(5.0/rho_c)*pow(rho/rho_c,4.0)+
					(B_16+B_26*(T/T_c))*(6.0/rho_c)*pow(rho/rho_c,5.0);
  return dDeltaKtdrho_T;
}

su2double CConductivityToluene::dDeltaKtdT_rho(su2double T,su2double rho) {
   su2double dDeltaKtdT_rho;
   dDeltaKtdT_rho = (B_21/T_c)*pow(rho/rho_c,1.0)+
					(B_22/T_c)*pow(rho/rho_c,2.0)+
					(B_23/T_c)*pow(rho/rho_c,3.0)+
					(B_24/T_c)*pow(rho/rho_c,4.0)+
					(B_25/T_c)*pow(rho/rho_c,5.0)+
					(B_26/T_c)*pow(rho/rho_c,6.0);
   return dDeltaKtdT_rho;
}
su2double CConductivityToluene::DeltaKt_c(su2double T,su2double rho) {
  su2double DeltaKt_c;
  su2double DeltaRho_c, DeltaT_c;
  DeltaRho_c = (rho/rho_c)-1;
  DeltaT_c = (T/T_c)-1;
  DeltaKt_c = (C_1/(C_2+abs(DeltaT_c)))*exp(-pow(C_3*DeltaRho_c, 2.0));
  return DeltaKt_c;
}

void CConductivityToluene::SetConductivity(su2double T, su2double rho, su2double mu, su2double cp) {
  Kt = Kt_0(T)+DeltaKt(T,rho);//+DeltaKt_c(rho,T);
}

void CConductivityToluene::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {
  dktdrho_T = dDeltaKtdrho_T(T,rho);
  dktdT_rho = dKt_0dT_rho(T)+dDeltaKtdT_rho(T,rho);
}

CConductivityToluene::~CConductivityToluene(void) { }

