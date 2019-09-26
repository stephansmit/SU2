/*!
 * \file CMeshVariable.hpp
 * \brief Declaration and inlines of the class
 *        to define the variables of the mesh movement.
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

#include "CVariable.hpp"

class CMeshVariable : public CVariable {
protected:

  unsigned short nDim;

  su2double WallDistance;   /*!< \brief Store the wall distance in reference coordinates. */

  su2double *Mesh_Coord;           /*!< \brief Store the reference coordinates of the mesh. */


public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_coor - Values of the coordinates (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CMeshVariable(const su2double *val_coor, unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMeshVariable(void);

  /*!
   * \brief Get the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the original coordinate iDim.
   */
  inline su2double GetMesh_Coord(unsigned short iDim) const final { return Mesh_Coord[iDim]; }

  /*!
   * \brief Get the undeformed coordinates.
   * \return Pointer to the reference coordinates.
   */
  inline su2double *GetMesh_Coord() final { return Mesh_Coord; }

  /*!
   * \brief Set the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \param[in] val_coord - Value of Mesh_Coord[nDim]
   */
  inline void SetMesh_Coord(unsigned short iDim, const su2double val_coord) final { Mesh_Coord[iDim] = val_coord;}

  /*!
   * \brief Move Displacement into Displacement_Old.
   */
  inline void SetSolution_Old(void){
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Solution_Old[iDim] = Solution[iDim];
  }

  /*!
   * \brief Get the value of the wall distance in reference coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the wall distance in reference coordinates.
   */
  inline su2double GetWallDistance(void) const final { return WallDistance; }

  /*!
   * \brief Set the value of the wall distance in reference coordinates.
   * \param[in] val_dist - Value of wall distance.
   */
  inline void SetWallDistance(const su2double val_dist) final { WallDistance = val_dist; }

  /*!
   * \brief Determine whether the node is a moving vertex.
   * \return False. The node is not at the boundary.
   */
  inline virtual bool Get_isVertex(void) const override { return false; }

  /*!
   * \brief Register the reference coordinates of the mesh.
   * \param[in] input - Defines whether we are registering the variable as input or as output.
   */
  inline void Register_MeshCoord(bool input) final {
    if (input) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++)
        AD::RegisterInput(Mesh_Coord[iVar]);
    }
    else { for (unsigned short iVar = 0; iVar < nVar; iVar++)
        AD::RegisterOutput(Mesh_Coord[iVar]);
    }
  }

  /*!
   * \brief Recover the value of the adjoint of the mesh coordinates.
   */
  inline void GetAdjoint_MeshCoord(su2double *adj_mesh) const final{
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        adj_mesh[iVar] = SU2_TYPE::GetDerivative(Mesh_Coord[iVar]);
    }
  }

};
