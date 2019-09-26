/*!
 * \file CDiscAdjMeshSolver.cpp
 * \brief Main subroutines for solving the discrete adjoint mesh problem.
 * \author Ruben Sanchez
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


#include "../../include/solvers/CDiscAdjMeshSolver.hpp"
#include "../../include/variables/CDiscAdjMeshVariable.hpp"
#include "../../include/variables/CDiscAdjMeshBoundVariable.hpp"

CDiscAdjMeshSolver::CDiscAdjMeshSolver(void) : CSolver (){

  KindDirect_Solver = 0;

  direct_solver = NULL;

}

CDiscAdjMeshSolver::CDiscAdjMeshSolver(CGeometry *geometry, CConfig *config)  : CSolver(){

  KindDirect_Solver = 0;

}

CDiscAdjMeshSolver::CDiscAdjMeshSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver)  : CSolver(){

  unsigned short iVar, iMarker, iDim;
  unsigned long iPoint;
  long iVertex;
  bool isVertex;

  nVar = geometry->GetnDim();
  nDim = geometry->GetnDim();

  /*-- Store some information about direct solver ---*/
  this->direct_solver = direct_solver;

  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliary vectors related to the residual for problems with a BGS strategy---*/

  if (config->GetMultizone_Residual()){

    Residual_BGS      = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]      = 1.0;
    Residual_Max_BGS  = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
    }

  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 1e-16;

  /*--- Initialize the node structure ---*/
  node       = new CVariable*[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- In principle, the node is not at the boundary ---*/
    isVertex = false;
    /*--- Looping over all markers ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      /*--- If the marker is flagged as moving, retrieve the node vertex ---*/
      if (config->GetMarker_All_Deform_Mesh(iMarker) == YES) iVertex = geometry->node[iPoint]->GetVertex(iMarker);
      else iVertex = -1;

      if (iVertex != -1){isVertex = true; break;}
    }

    /*--- The MeshBound variable includes the displacements at the boundaries ---*/
    if (isVertex) node[iPoint] = new CDiscAdjMeshBoundVariable(Solution, nDim, config);
    else          node[iPoint] = new CDiscAdjMeshVariable(Solution, nDim, config);

  }


}

CDiscAdjMeshSolver::~CDiscAdjMeshSolver(void){

}


void CDiscAdjMeshSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config_container, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output){


}

void CDiscAdjMeshSolver::SetRecording(CGeometry* geometry, CConfig *config){


  unsigned long iPoint;
  /*--- Reset the solution to the initial (converged) solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->SetBound_Disp(node[iPoint]->GetBoundDisp_Direct());
  }

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void CDiscAdjMeshSolver::RegisterSolution(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, nPoint = geometry->GetnPoint();
  bool input = true;

  /*--- Register reference mesh coordinates ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->Register_MeshCoord(input);
  }

}

void CDiscAdjMeshSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset){

  /*--- Register boundary displacements as input ---*/

  unsigned long iPoint, nPoint = geometry->GetnPoint();
  bool input = true;

  /*--- Register reference mesh coordinates ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->Register_BoundDisp(input);
  }

}

void CDiscAdjMeshSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  /*--- Extract the sensitivities of the mesh coordinates ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution from the original mesh coordinates ---*/

    direct_solver->node[iPoint]->GetAdjoint_MeshCoord(Solution);

    /*--- Store the adjoint solution (the container is reused) ---*/

    node[iPoint]->SetSolution(Solution);

  }

}

void CDiscAdjMeshSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  /*--- Extract the sensitivities of the boundary displacements ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution of the boundary displacements ---*/

    direct_solver->node[iPoint]->GetAdjoint_BoundDisp(Solution);

    /*--- Store the sensitivities of the boundary displacements ---*/

    node[iPoint]->SetBoundDisp_Sens(Solution);

  }

}

void CDiscAdjMeshSolver::SetSensitivity(CGeometry *geometry, CSolver **solver, CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;
  su2double Sensitivity, eps;
  bool time_stepping = (config->GetUnsteady_Simulation() != STEADY);

  /*--- Extract the sensitivities ---*/
  ExtractAdjoint_Solution(geometry, config);

  /*--- Extract the adjoint variables: sensitivities of the boundary displacements ---*/
  ExtractAdjoint_Variables(geometry, config);

  /*--- Store the sensitivities in the flow adjoint container ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    for (iDim = 0; iDim < nDim; iDim++) {

      /*--- The sensitivity was extracted using ExtractAdjoint_Solution ---*/
      Sensitivity = node[iPoint]->GetSolution(iDim);

      /*--- If sharp edge, set the sensitivity to 0 on that region ---*/
      if (config->GetSens_Remove_Sharp()) {
        eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
        if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetAdjSharp_LimiterCoeff()*eps )
          Sensitivity = 0.0;
      }

      /*--- Store the sensitivities ---*/
      if (!time_stepping) {
        solver[ADJFLOW_SOL]->node[iPoint]->SetSensitivity(iDim, Sensitivity);
      } else {
        solver[ADJFLOW_SOL]->node[iPoint]->SetSensitivity(iDim, solver[ADJFLOW_SOL]->node[iPoint]->GetSensitivity(iDim) + Sensitivity);
      }
    }
  }
  solver[ADJFLOW_SOL]->SetSurface_Sensitivity(geometry, config);

}

void CDiscAdjMeshSolver::ComputeResidual_Multizone(CGeometry *geometry, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++){
      SetRes_BGS(iVar,0.0);
      SetRes_Max_BGS(iVar,0.0,0);
  }

  /*--- Set the residuals ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    /*--- Only for the boundary vertices ---*/
    if (node[iPoint]->Get_isVertex()){
      for (iVar = 0; iVar < nVar; iVar++){
          /*--- Compute only for the sensitivities of the boundary displacements ---*/
          residual = node[iPoint]->GetBoundDisp_Sens(iVar) - node[iPoint]->Get_BGSSolution_k(iVar);
          AddRes_BGS(iVar,residual*residual);
          AddRes_Max_BGS(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
    }
  }

  SetResidual_BGS(geometry, config);

}


void CDiscAdjMeshSolver::UpdateSolution_BGS(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  /*--- To nPoint: The solution must be communicated beforehand ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    if (node[iPoint]->Get_isVertex()) node[iPoint]->Set_BGSSolution_k();
  }

}


void CDiscAdjMeshSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

}

