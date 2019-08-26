/*!
 * \file CIncTGVSolution.cpp
 * \brief Implementations of the member functions of CIncTGVSolution.
 * \author T. Economon, E. van der Weide
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "../../../include/toolboxes/MMS/CIncTGVSolution.hpp"

CIncTGVSolution::CIncTGVSolution(void) : CVerificationSolution() { }

CIncTGVSolution::CIncTGVSolution(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 unsigned short val_iMesh,
                                 CConfig*       config)
: CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
    
  /*--- Disable this solution for now, as it has not been tested. ---*/
  
  SU2_MPI::Error("CIncTGVSolution not yet fully implemented/tested.",
                 CURRENT_FUNCTION);
  
  /*--- Write a message that the solution is initialized for the
   Taylor-Green vortex test case. ---*/
  
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the incompressible Taylor-Green vortex case!!!" << endl;
    cout << endl << flush;
  }
  
  /*--- Store TGV specific parameters here. ---*/
  
  tgvLength    = 1.0;
  tgvVelocity  = 1.0;
  tgvDensity   = config->GetDensity_FreeStreamND();
  tgvViscosity = config->GetViscosity_FreeStreamND();
  
  /*--- We keep a copy of the freestream temperature just to be safe
   when we set the solution, even though this is an isothermal case. ---*/
  
  Temperature = config->GetTemperature_FreeStreamND();
  
  /*--- Perform some sanity and error checks for this solution here. ---*/
  
  if((config->GetUnsteady_Simulation() != TIME_STEPPING) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_1ST) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(Kind_Solver != INC_EULER && Kind_Solver != INC_NAVIER_STOKES && Kind_Solver != INC_RANS )
    SU2_MPI::Error("Incompressible flow equations must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(Kind_Solver != INC_NAVIER_STOKES)
    SU2_MPI::Error("Navier Stokes equations must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_FluidModel() != CONSTANT_DENSITY)
    SU2_MPI::Error("Constant density fluid model must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
    SU2_MPI::Error("Constant viscosity must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(config->GetEnergy_Equation())
    SU2_MPI::Error("Energy equation must be disabled (isothermal) for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(nDim != 2)
    SU2_MPI::Error("2D calculation required for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
}

CIncTGVSolution::~CIncTGVSolution(void) { }

void CIncTGVSolution::GetBCState(const su2double *val_coords,
                                 const su2double val_t,
                                 su2double       *val_solution) {
  
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CIncTGVSolution::GetSolution(const su2double *val_coords,
                                  const su2double val_t,
                                  su2double       *val_solution) {
  
  /* The exact solution is set for the incompressible Taylor-Green
   vortex case. This is the classic solution from the original work
   of Taylor and Green for the specific 2D situation where the
   exact solution can be derived for an incompressible flow. */
  
  /* Store the termporal term more easily (Taylor expansion). */
  su2double F = 1.0 - 2.0*(tgvViscosity/tgvDensity)*val_t;
  
  /* Compute the primitive variables. */
  su2double u =  tgvVelocity * F * (sin(val_coords[0]/tgvLength)*
                                    cos(val_coords[1]/tgvLength));
  su2double v = -tgvVelocity * F * (cos(val_coords[0]/tgvLength)*
                                    sin(val_coords[1]/tgvLength));
  
  su2double B = (cos(2.0*val_coords[0]/tgvLength) +
                 cos(2.0*val_coords[1]/tgvLength));
  su2double p = -(tgvDensity/4.0)*B*F*F;
  
  /* Compute the conservative variables. Note that both 2D and 3D
   cases are treated correctly. */
  val_solution[0]      = p;
  val_solution[1]      = u;
  val_solution[2]      = v;
  val_solution[nVar-1] = Temperature;
}
