/*!
 * \file CMeshSolver.hpp
 * \brief Declaration and inlines of the class to compute the 
 *        the discrete adjoint of the linear-elastic mesh solver.
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

#pragma once

#include "../solver_structure.hpp"

/*!
 * \class CDiscAdjMeshSolver
 * \brief Main class for defining the discrete adjoint solver for mesh deformation problems.
 * \ingroup Discrete_Adjoint
 * \author R. Sanchez
 */
class CDiscAdjMeshSolver : public CSolver {
private:
  unsigned short KindDirect_Solver;
  CSolver *direct_solver;

public:

  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjMeshSolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CDiscAdjMeshSolver(CGeometry *geometry, CConfig *config);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Initialize the discrete adjoint solver with the corresponding direct solver.
   * \param[in] Kind_Solver - The kind of direct solver.
   */
  CDiscAdjMeshSolver(CGeometry *geometry, CConfig *config, CSolver* solver);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjMeshSolver(void);

  /*!
   * \brief Performs the preprocessing of the AD-based mesh adjoint solver.
   *        Registers all necessary variables on the tape. Called while tape is active.
   * \param[in] geometry_container - The geometry container holding all grid levels.
   * \param[in] config_container - The particular config.
   */
  void RegisterSolution(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Sets the adjoint values of the input variables of the flow (+turb.) iteration
   *        after tape has been evaluated.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Extract and set the geometrical sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The solver container holding all terms of the solution.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSensitivity(CGeometry *geometry, CSolver **solver, CConfig *config);

  /*!
   * \brief Set the value of the max residual and RMS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store the BGS solution in the previous subiteration in the corresponding vector.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void UpdateSolution_BGS(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Prepare the solver for a new recording.
   * \param[in] kind_recording - Kind of AD recording.
   */
  void SetRecording(CGeometry *geometry, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void RegisterVariables(CGeometry *geometry, CConfig *config, bool reset = false);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the dual-time derivatives.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);

};
