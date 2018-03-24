/*!
 * \file wall_model.cpp
 * \brief Function for the wall model functions for hom large eddy simulations.
 * \author E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
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

#include "../include/wall_model.hpp"

CWallModel::CWallModel(void){
  thickness = 0.0;
}

CWallModel::~CWallModel(void){}

void CWallModel::ComputeWallShear(const unsigned short val_marker){}

void CWallModel::ComputeWallHeatFlux(const unsigned short val_marker){}

void CWallModel::Initialize(CBoundaryFEM * boundary, CConfig *config, CGeometry *geometry){}

void CWallModel::SetUpExchange(CBoundaryFEM * boundary, CConfig *config, CGeometry *geometry){}

void CWallModel::SetPoints(CSurfaceElementFEM * thisElem, vector<su2double> exchangeCoords){}

CWallModel1DEQ::CWallModel1DEQ(void) : CWallModel(){
  expansionRatio = 0.0;
  numPoints = 0.0;
}

CWallModel1DEQ::~CWallModel1DEQ(void){}

void CWallModel1DEQ::ComputeWallShear(const unsigned short val_marker){}

void CWallModel1DEQ::ComputeWallHeatFlux(const unsigned short val_marker){}

void CWallModel1DEQ::Initialize(CBoundaryFEM * boundary, CConfig *config, CGeometry *geometry){

  /*--- Set up the wall model for this boundary marker ---*/

  /*--- Get the marker name for this boundary ---*/
  std::string markerName = boundary->markerTag;

  /*--- Get double data from the configuration class object. This should include
   * the wall model thickness and the expansion ratio of the points ---*/
  su2double * tempDoubleData = config->GetWallFunction_DoubleInfo(markerName);
  thickness = tempDoubleData[0];
  expansionRatio = tempDoubleData[1];

  initCondExists = false;

  // Output values for debugging purposes
  std::cout << "Thickness = " << thickness << std::endl;
  std::cout << "Expansion Ratio = " << expansionRatio << std::endl;

  /*--- Get the integer data from the configuration class object. This includes
   * the number of points to be used in the wall model ---*/
  unsigned short int * tempIntData = config->GetWallFunction_IntInfo(markerName);
  numPoints = tempIntData[0];

  // Output the number of points for debugging purposes
  std::cout << "Number of WM Points = " << numPoints << std::endl;
  std::cout << "*******************" << std::endl;

  /*--- Initialize some things on the boundary, itself ---*/
  boundary->wallModelBoundary = true;
  boundary->nWallModelPoints = numPoints;
  boundary->wallModelExpansionRatio = expansionRatio;
  boundary->wallModelThickness = thickness;

  this->SetUpExchange(boundary, config, geometry);

}

void CWallModel1DEQ::SetUpExchange(CBoundaryFEM * boundary, CConfig *config, CGeometry *geometry){

  /*--- This method will set up the exchange locations in each boundary cell. ---*/

  /*--- Create the DGGeometry pointer by casting the geometry as a CMeshFEM_DG ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  // Get number of dimensions
  unsigned short nDim = DGGeometry->GetnDim();

  // Get array of volume elements
  CVolumeElementFEM * volElements = DGGeometry->GetVolElem();

  /*--- Loop over the boundary elements ---*/
  vector<CSurfaceElementFEM> surfElem = boundary->surfElem;
  for(unsigned long iSurfElem = 0; iSurfElem < surfElem.size(); ++iSurfElem) {

    // Get a pointer to current surface element
    CSurfaceElementFEM * curElem = &(surfElem[iSurfElem]);

    // Get the number of degrees of freedom of this face and the corresponding element
    const unsigned short nDofsFace = curElem->DOFsSolFace.size();
    const unsigned short nDofsElement = curElem->DOFsSolElement.size();
    const unsigned long volID = curElem->volElemID;

    // Set the size of the exchange point ID vector and initialize to zero
    curElem->exchangePointIDs.assign(nDofsFace,0);

    // Get the current face's corresponding volume element
    CVolumeElementFEM * curVolume = &(volElements[volID]);

    // Get the current corresponding volume element's polynomial degree
    unsigned short nPoly = curVolume->nPolySol;

    // Step through the dofs of the surface element
    for(unsigned short iDof = 0; iDof < nDofsFace; ++iDof){

      // Set the initial difference to be the specified thickness
      su2double difference = this->thickness;

      // For each surface DOF, there will be nPoly solution DOFs "above" it in the wall-normal direction
      for(unsigned short i = 0; i < nPoly; ++i){
        // We want to look at the nPoly solution DOFs above the current surface DOF. These solution DOFs will be
        // at the current index plus (nPoly+1)^2, 2*(nPoly+1)^2, 3*(nPoly+1)^2, etc in the vector.

        // Calculate the index of the ith solution DOF above the current surface DOF
        unsigned short thisSolDOFIndex = iDof + (i+1)*((nPoly+1)*(nPoly+1));

        // Get the wall distance of this solution DOF
        su2double thisWallDist = curVolume->wallDistanceSolDOFs[thisSolDOFIndex];

        // Compare the wall distance of this DOF with that specified by the wall model for this wall
        // If the difference is smaller than the current difference, this point is closer to the specified
        // wall thickness
        if( fabs(this->thickness - thisWallDist) < difference ){
          // Set the difference between the specified wall model thickness and this solution DOF to be the new difference
          difference = fabs(this->thickness - thisWallDist);

          // Set the current surface DOF's exchange point ID to that of the current solution DOF
          curElem->exchangePointIDs[iDof] = thisSolDOFIndex;
        }
      }
    }
    su2double exchangeHeight = curVolume->wallDistanceSolDOFs[curElem->exchangePointIDs[0]];
    this->SetPoints(curElem,exchangeHeight);
  }
}

void CWallModel1DEQ::SetPoints(CSurfaceElementFEM * thisElem, su2double exchangeHeight){
  /*--- Allocate memory/vector of wall model points for this element. This assumes that there is only
   * one vector of wall model points needed per cell. This assumption is based on the volume cell being
   * a regular hex or triangular prism, where the distances from each wall integration point to its
   * corresponding exchange locations is the same throughout this boundary cell. This may not be a good
   * assumption in general cases, but will work for plane channel flow for now. ---*/

  thisElem->wallModelPoints.resize(this->numPoints);

  /*--- for now, we are just using the y-coordinate ---*/
  /*--- Use a geometric expansion to get wall model point coords between the wall and the exchange location ---*/
  /*--- Find the first cell thickness ---*/
  su2double firstCellThickness = exchangeHeight * (1 - this->expansionRatio) / (1 - std::pow(this->expansionRatio,this->numPoints-1));

  // Output for debugging
  if(thisElem->boundElemIDGlobal == 0)
  {
    std::cout << "For boundElemIDGlobal = " << thisElem->boundElemIDGlobal << std::endl;
    std::cout << "Specified WM thickness = " << this->thickness << std::endl;
    std::cout << "Actual exchange wall distance = " << exchangeHeight << std::endl;
  }

  // Set first wall model point to be at the wall (y=0)
  thisElem->wallModelPoints[0] = 0.0;
  for(unsigned short iPoint = 1; iPoint < numPoints; ++iPoint){
    thisElem->wallModelPoints[iPoint] = thisElem->wallModelPoints[iPoint-1] + firstCellThickness * std::pow(this->expansionRatio,iPoint-1);
  }
  if(thisElem->boundElemIDGlobal == 0){
    for(unsigned short i = 0; i < numPoints; ++i){
      std::cout << "y[" << i << "] = " << thisElem->wallModelPoints[i] << std::endl;
    }
  }
}

void CWallModel1DEQ::SolveCoupledSystem(CSurfaceElementFEM * thisElem, su2double u_bc, su2double T_bc, su2double P_bc, su2double * wallShear, su2double * heatFlux){
  /*--- This function solves the coupled system of Tri-diagonal systems ---*/

  /*--- First, set up vectors ---*/
  std::vector<su2double> y(numPoints,0.0);
  std::vector<su2double> yPlus(numPoints,0.0);
  std::vector<su2double> u(numPoints,u_bc);
  std::vector<su2double> T(numPoints,T_bc);
  std::vector<su2double> mu(numPoints,0.0);
  std::vector<su2double> muTurb(numPoints,0.0);
  std::vector<su2double> energyFlux(numPoints,0.0);
  su2double tauWall;

  bool converged = false;

  while( converged == false ){
    /*--- Set initial condition if it doesn't already exist ---*/
    if( initCondExists == false ){
      for(unsigned short i = 0; i<numPoints; i++){
        /*--- Set the viscosity according to Sutherland's law ---*/
        su2double C_1 = 1.458e-6;
        su2double S = 110.4;
        mu[i] = C_1 * std::pow(T[i],1.5) / (T[i] + S);

        /*--- Set the yPlus vector based on the molecular viscosity ---*/
        su2double R = 287.058;
        su2double rho = P_bc / (R * T[i]); // Pressure is assumed constant through the WM
        su2double u_tau = std::sqrt(tauWall/rho);

        /*--- Set the turbulent viscosity ---*/
        su2double kappa = 0.41;
        su2double A = 17;
        muTurb[i] = kappa * y[i] * std::sqrt(rho*tauWallInit) * std::pow((1 - std::exp(-1*yPlus[i]/A)),2);
      }
    }


  }






}

void CWallModel1DEQ::GetExchangeValues(unsigned short intPointID, su2double * u_exchange, su2double * T_exchange,  su2double * P_exchange){

}
